module SCS

using Libdl
using Requires

if haskey(ENV, "JULIA_SCS_LIBRARY_PATH") || VERSION < v"1.3"
    if isfile(joinpath(dirname(@__FILE__), "..", "deps", "deps.jl"))
        include(joinpath(dirname(@__FILE__), "..", "deps", "deps.jl"))
    else
        error("SCS not properly installed. Please run Pkg.build(\"SCS\") and restart julia")
    end

    function __init__()
        vnum = VersionNumber(SCS_version())
        depsdir = realpath(joinpath(dirname(@__FILE__),"..","deps"))
        if vnum < v"2.1.0"
            error("Current SCS version installed is $vnum, but we require version 2.1.*")
        end
    end
else
    import SCS_jll
    const indirect = SCS_jll.libscsindir
    const direct = SCS_jll.libscsdir

    function __init__()
        @require CUDA_jll="e9e359dc-d701-5aa8-82ae-09bbf812ea83" include("c_wrapper_gpu.jl")
    end
end

include("types.jl")

const available_solvers = [DirectSolver, IndirectSolver]

include("c_wrapper.jl")
include("MPB_wrapper.jl")
include("MOI_wrapper.jl")

if Base.VERSION >= v"1.4.2"
    include("precompile.jl")
    _precompile_()
end

end # module
