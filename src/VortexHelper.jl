module VortexHelper

include("run_vortex.jl")
export run_vortex,
    set_test_runs,
    vortex_output_directory

include("read_vortex.jl")
export read_particle_data

end # module
