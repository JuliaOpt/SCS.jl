module SCS

using Libdl
using Requires

if haskey(ENV, "JULIA_SCS_LIBRARY_PATH") || VERSION < v"1.3"
    haskey(ENV, "JULIA_SCS_LIBRARY_PATH") && @info "using local SCS libraries at $(ENV["JULIA_SCS_LIBRARY_PATH"])"
    if isfile(joinpath(dirname(@__FILE__), "..", "deps", "deps.jl"))
        include(joinpath(dirname(@__FILE__), "..", "deps", "deps.jl"))
    else
        error("SCS not properly installed. Please run Pkg.build(\"SCS\") and restart julia")
    end

    function __init__()
        vnum = VersionNumber(SCS_version())
        if vnum < v"3.0.0"
            error("Current SCS version installed is $vnum, but we require version 3.0.*")
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

end # module
