using BinDeps

@BinDeps.setup

blasvendor = Base.BLAS.vendor()

libs = Dict(:direct=>"libscsdir", :indirect=>"libscsindir")

if (is_apple() ? (blasvendor == :openblas64) : false)
    aliases = ["libscsdir64"]
else
    aliases = ["libscsdir"]
end

direct = library_dependency("direct", aliases=[libs[:direct]])
indirect = library_dependency("indirect", aliases=[libs[:indirect]])

if is_apple()
    using Homebrew
    provides(Homebrew.HB, "scs", scs, os = :Darwin)
end

version = "2.0.2"
win_version = "2.0.2" # The windows binaries are not consistent with this version yet.

provides(Sources, URI("https://github.com/cvxgrp/scs/archive/v$version.tar.gz"),
    [direct, indirect], os=:Unix, unpacked_dir="scs-$version")

# Windows binaries built in Cygwin as follows:
# CFLAGS="-DDLONG -DCOPYAMATRIX -DUSE_LAPACK -DCTRLC=1 -DBLAS64 -DBLASSUFFIX=_64_" LDFLAGS="-L$HOME/julia/usr/bin -lopenblas64_" make CC=x86_64-w64-mingw32-gcc out/libscsdir.dll
# mv out bin64
# make clean
# CFLAGS="-DDLONG -DCOPYAMATRIX -DUSE_LAPACK -DCTRLC=1" LDFLAGS="-L$HOME/julia32/usr/bin -lopenblas" make CC=i686-w64-mingw32-gcc out/libscsdir.dll
# mv out bin32
provides(Binaries, URI("https://cache.julialang.org/https://bintray.com/artifact/download/tkelman/generic/scs-$win_version-r2.7z"),
    [scs], unpacked_dir="bin$(Sys.WORD_SIZE)", os = :Windows,
    SHA="62bb4feeb7d2cd3db595f05b86a20fc93cfdef23311e2e898e18168189072d02")

prefix = joinpath(BinDeps.depsdir(direct), "usr")
srcdir = joinpath(BinDeps.depsdir(direct), "src", "scs-$version/")

libnames = Dict(k => v*".$(Libdl.dlext)" for (k,v) in libs)

ldflags = ""
if is_apple()
    ldflags = "$ldflags -undefined suppress -flat_namespace"
end
cflags = "-DCOPYAMATRIX -DDLONG -DUSE_LAPACK -DCTRLC=1"
if blasvendor == :openblas64
    cflags = "$cflags -DBLAS64 -DBLASSUFFIX=_64_"
end
if blasvendor == :mkl
    if Base.USE_BLAS64
        cflags = "$cflags -DMKL_ILP64 -DBLAS64"
        ldflags = "$ldflags -lmkl_intel_ilp64"
    else
        ldflags = "$ldflags -lmkl_intel"
    end
    cflags = "$cflags -fopenmp"
    ldflags = "$ldflags -lmkl_gnu_thread -lmkl_rt -lmkl_core -lgomp"
end

ENV2 = copy(ENV)
ENV2["LDFLAGS"] = ldflags
ENV2["CFLAGS"] = cflags

provides(SimpleBuild,
    (@build_steps begin
        GetSources(direct)
        CreateDirectory(joinpath(prefix, "lib"))

        FileRule(joinpath(prefix, "lib", libnames[:direct]), @build_steps begin
            ChangeDirectory(srcdir)
            setenv(`make BLASLDFLAGS= out/$(libnames[:direct])`, ENV2)
            `mv out/$(libnames[:direct]) $prefix/lib`
        end)

        FileRule(joinpath(prefix, "lib", libnames[:indirect]), @build_steps begin
            ChangeDirectory(srcdir)
            setenv(`make BLASLDFLAGS= out/$(libnames[:indirect])`, ENV2)
            `mv out/$(libnames[:indirect]) $prefix/lib`
        end)

    end), [direct, indirect], os=:Unix)

@BinDeps.install Dict(:direct => :direct, :indirect => :indirect)
