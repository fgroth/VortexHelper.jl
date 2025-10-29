using GadgetIO
using Formatting

global test_runs = "/home/moon/fgroth/phd/test_runs/test_collection/test_runs/"
"""
    set_test_runs(dir::String="/home/moon/fgroth/phd/test_runs/test_collection/test_runs/")

Set the global `test_runs` variable that is used to locate the output of simulations.
"""
function set_test_runs(dir::String="/home/moon/fgroth/phd/test_runs/test_collection/test_runs/")
    global test_runs = dir
end

"""
    run_vortex(cluster::String, method::String; 
               start_snap_num::Int64=100,end_snap_num::Int64=145, scale::Number=3, snaps_todo=nothing, 
               vortex_directory::String="vortex-p",
               limit_resources::Bool=true,
               slurm_submission::Bool=false, keep_tmp_dir::Bool=false,
               # adjust the vortex parameters
               filtering_length::Real=0, smooth_filtering::Bool=false,
               cells_per_direction::Int64=128,
               n_levels::Int64=9, n_particles_refinement::Int64=8, minimum_refinement_patchsize::Int64=-1,
               weighting::String="")

Run Vortex for given snapshots of `cluster` and `method`.

`filtering_length<0` means multi-scale filtering, typical values (corresponding to negative of maximum filtering legnth) are -1000.0. `filtering_length>0` means constant size filtering, and `filtering_length=0` un-filtered.

If `weighting` is an empty string, the default is assumed (no further suffix to executable). Otherwise, it is added as suffix to the executable and the output directory.
"""
function run_vortex(cluster::String, method::String;
                    start_snap_num::Int64=100,end_snap_num::Int64=145, scale::Number=3, snaps_todo=nothing,
                    vortex_directory::String="vortex-p",
                    limit_resources::Bool=true,
                    slurm_submission::Bool=false, keep_tmp_dir::Bool=false,
                    # adjust the vortex parameters
                    filtering_length::Real=0, smooth_filtering::Bool=false,
                    cells_per_direction::Int64=128,
                    n_levels::Int64=9, n_particles_refinement::Int64=8, minimum_refinement_patchsize::Int64=-1,
                    weighting::String="") # empty string uses default weighting
    
    last_snapnum=zeros(Int64,1)
    for snapnum in end_snap_num:-1:start_snap_num
        if isdir(joinpath(test_runs, "out_"*cluster*"_"*method, "snapdir_"*sprintf1("%03d",snapnum)))
            last_snapnum[1] = snapnum
            break
        else
            continue
        end
    end

    println("last snap num is ",last_snapnum[1])
    
    this_dir = pwd()

    # create temporary directory to run vortex.
    tmp_dir = mktempdir("./", cleanup=(!slurm_submission&&!keep_tmp_dir)) # make sure the path is not deleted if running in slurm mode, as we might close the julia session.
    vortex_exec = if occursin("mfm",method)
        "vortex_mfm"
    else
        "vortex_sph"
    end
    vortex_exec = if filtering_length != 0
        vortex_exec * "_filtered"
    else
        vortex_exec * "_unfiltered"
    end
    if weighting != ""
        vortex_exec = vortex_exec*"_"*weighting
    end
    cp(joinpath(vortex_directory, "src"), joinpath(tmp_dir, "src"))
    cd(joinpath(tmp_dir, "src"))
    println("working in temporary directory ",tmp_dir)

    # this is where vortex takes the input data from
    symlink(joinpath(test_runs, "out_"*cluster*"_"*method),"./simulation")

    vortex_output = vortex_output_directory(cluster, method, filtering_length=filtering_length, weighting=weighting)
    try
        mkdir(vortex_output)
    catch
        println(vortex_output*" already exists")
        # directory already exists
    end

    snaps_todo = if snaps_todo == nothing
        start_snap_num:last_snapnum[1]
    else
        snaps_todo
    end

    # prepare the executable to run vortex
    run_sh = open("run.sh","w")
    write(run_sh, "#!/bin/bash\n")
    if slurm_submission
        write(run_sh, "#SBATCH -J vortex_"*cluster*"_"*method*" # name of the job\n")
        write(run_sh, "#SBATCH -o ../../%x.%j.out            # output log file with name <job_name>.<job_id>.out\n")
        write(run_sh, "#SBATCH -e ../../%x.%j.err            # error log file with name <job_name>.<job_id>.err\n")
        write(run_sh, "#SBATCH -D ./                         # output directory\n")
        write(run_sh, "#SBATCH --nodes=1                     # number of nodes\n")
        write(run_sh, "#SBATCH --ntasks-per-node=1           # number of MPI ranks per node\n")
        write(run_sh, "#SBATCH --cpus-per-task=80            # number of OpenMP threads per MPI rank\n")
        write(run_sh, "#SBATCH --time=3-00:00:00             # time limit of the run\n")
        write(run_sh, "#SBATCH --mem=200GB                   # maximum mermory to be used by the job\n\n")
    end
    write(run_sh, "\n")
    if limit_resources
        write(run_sh, "ulimit -s 200000000\n")
        write(run_sh, "ulimit -v unlimited\n")
        write(run_sh, "ulimit -c 0\n")
    else
        write(run_sh, "ulimit -s unlimited\n")
        write(run_sh, "ulimit -v unlimited\n")
        write(run_sh, "ulimit -c 0\n")
    end
    write(run_sh, "echo running in\n")
    write(run_sh, "echo | pwd\n")
    write(run_sh, "export OMP_NUM_THREADS=32\n")
    write(run_sh, "export OMP_STACKSIZE=4000m\n")
    write(run_sh, "export OMP_PROC_BIND=true\n")
    write(run_sh, "\n")
    write(run_sh,"./"*vortex_exec*"\n")
    if slurm_submission
        write(run_sh, "mv -f output_files "*joinpath(vortex_output, sprintf1("%03d",snaps_todo[1]))*"\n")
        write(run_sh, "cd "*this_dir*"\n")
        if !keep_tmp_dir
            write(run_sh, "rm -rf "*tmp_dir*"\n")
        end
    end
    close(run_sh)
    chmod("run.sh",0o700)
    
    if length(snaps_todo) > 1 && slurm_submission
        error("slurm_submission does only work for single snapshot at the moment.")
    end
    
    for i_snap in snaps_todo
        println("running ",i_snap)
        snap = joinpath(test_runs, "out_"*cluster*"_"*method, "snapdir_"*sprintf1("%03d",i_snap), "snap_"*sprintf1("%03d",i_snap))
        sub = joinpath(test_runs * "out_"*cluster*"_"*method, "groups_"*sprintf1("%03d",i_snap), "sub_"*sprintf1("%03d",i_snap))

        halo_positions = read_subfind(sub, "GPOS")
        
        first_halo_position = halo_positions[:,1]
        println("first halo position ",first_halo_position)
        
        # prepare the vortex parameter file
        par_name = "./vortex.dat"
        this_par = open(par_name,"w")

        println(this_par, "***********************************************************************")
        println(this_par, "*                  VORTEX-GADGET PARAMETERS FILE                      *")
        println(this_par, "***********************************************************************")
        println(this_par, "*       General parameters block                                      *")
        println(this_par, "***********************************************************************")
        println(this_par, "Files: first, last, every, num files per snapshot -------------------->")
        # adjust snap number
        n_snap = sum(startswith.(readdir(joinpath(test_runs, "out_"*cluster*"_"*method, "snapdir_"*sprintf1("%03d",i_snap))), "snap_"*sprintf1("%03d",i_snap)))
        println(this_par, sprintf1("%d",i_snap)*","*sprintf1("%d",i_snap)*",1,"*sprintf1("%d",n_snap)*"")
        println(this_par, "Cells per direction (NX,NY,NZ) --------------------------------------->")
        println(this_par, sprintf1("%d",cells_per_direction)*","*sprintf1("%d",cells_per_direction)*","*sprintf1("%d",cells_per_direction)*"")
        println(this_par, "Max box sidelength (in input length units) --------------------------->")
        # adjust size
        size = 25e3
        println(this_par, sprintf1("%f",2*size)*"")
        println(this_par, "Domain to keep particles (in input length units; x1,x2,y1,y2,z1,z2) -->")
        println(this_par, sprintf1("%f",first_halo_position[1]-size)*","*sprintf1("%f",first_halo_position[1]+size)*","*
            sprintf1("%f",first_halo_position[2]-size)*","*sprintf1("%f",first_halo_position[2]+size)*","*
            sprintf1("%f",first_halo_position[3]-size)*","*sprintf1("%f",first_halo_position[3]+size)*"")
        println(this_par, "***********************************************************************")
        println(this_par, "*       Output customisation (0=no, 1=yes)                            *")
        println(this_par, "***********************************************************************")
        println(this_par, "Gridded data: kernel length, density (mutually exclusive), velocity -->")
        println(this_par, "0,1,1")
        println(this_par, "Gridded results: vcomp, vsol, scalar_pot, vector_pot, div(v), curl(v)->")
        println(this_par, "1,1,1,1,1,1")
        println(this_par, "Particle results: interpolation error, particle-wise results --------->")
        println(this_par, "1,1")
        println(this_par, "Filter: gridded Mach/ABVC, shocked cells, filtering length, vturb ---->")
        if filtering_length != 0
            println(this_par, "1,1,1,1")
        else
            println(this_par, "0,0,0,0")
        end
        println(this_par, "***********************************************************************")
        println(this_par, "*       Mesh creation parameters                                      *")
        println(this_par, "***********************************************************************")
        println(this_par, "Number of levels ----------------------------------------------------->")
        println(this_par, sprintf1("%d",n_levels)*"")
        println(this_par, "Number of particles for a cell to be refinable ----------------------->")
        println(this_par, sprintf1("%d",n_particles_refinement)*"")
        println(this_par, "Minimum size of a refinement patch to be accepted -------------------->")
        println(this_par, sprintf1("%d", minimum_refinement_patchsize))
        println(this_par, "Cells not to be refined from the border (base grid) ------------------>")
        println(this_par, "2")
        println(this_par, "***********************************************************************")
        println(this_par, "*       Velocity interpolation parameters                             *")
        println(this_par, "***********************************************************************")
        println(this_par, "Number of neighbours for interpolation ------------------------------->")
        if occursin("mfm",method)
            println(this_par, "32")
        else
            println(this_par, "295")
        end
        println(this_par, "***********************************************************************")
        println(this_par, "*       Poisson solver                                                *")
        println(this_par, "***********************************************************************")
        println(this_par, "SOR presion parameter, SOR max iter, border for AMR patches ---------->")
        println(this_par, "1e-9,1000,2")
        println(this_par, "***********************************************************************")
        println(this_par, "*       Turbulent filter                                              *")
        println(this_par, "***********************************************************************")
        println(this_par, "Apply filter (1: multiscale filter; 2: fix-scale filter) ------------->")
        if filtering_length == 0
            filtering = 0
        elseif filtering_length < 0
            filtering = 1
        else # filtering_length > 0
            filtering = 2
        end
        println(this_par, sprintf1("%d", filtering)*"")
        println(this_par, "Filtering parameters: tolerance, growing step, max. num. of its. ----->")
        println(this_par, "0.1,1.05,200")
        println(this_par, "Maximum (for multiscale) or fix filt. length (input length units) ---->")
        println(this_par, sprintf1("%g",abs(filtering_length)))
        println(this_par, "Smooth filtering length before applying the filter (0=no, 1=yes) ----->")
        println(this_par, smooth_filtering ? "1" : "0")
        println(this_par, "***********************************************************************")
        println(this_par, "*       On-the-fly shock detection (for multifiltering)               *")
        println(this_par, "***********************************************************************")
        println(this_par, "Threshold on velocity divergence (negative, input units) ------------->")
        println(this_par, "-1.25")
        println(this_par, "Threshold on artificial bulk viscosity constant ---------------------->")
        println(this_par, "1.")
        println(this_par, "Use particle's MACH field (0=no, 1=yes), Mach threshold -------------->")
        println(this_par, "1,2.0")

        close(this_par)

        # this directory is required by vortex, all outputfiles are stored there
        mkdir("output_files")
        # run vortex using the executable created above
        if slurm_submission
            run(`sbatch run.sh`)
        else
            run(`./run.sh`)
            # move the output files to a standardized directory
            mv("output_files",joinpath(vortex_output, sprintf1("%03d",i_snap)),force=true)
        end
    end
    
    # change back to original directory and clean up
    cd(this_dir)
    if !slurm_submission && !keep_tmp_dir
        rm(tmp_dir,force=true, recursive=true)
    end 
end

"""
    vortex_output_directory(cluster::AbstractString, method::AbstractString;
                            filtering_length::Real=-1, weighting::AbstractString="")

Return desired location of directory containing vortex output.
"""
function vortex_output_directory(cluster::AbstractString, method::AbstractString;
                                 filtering_length::Real=-1, weighting::AbstractString="")
    if filtering_length < 0
        prefix = "multi-filtered_"
    elseif filtering_length > 0
        prefix = "filtered-"*sprintf1("%g",filtering_length)*"_"
    else # filtering_length == 0
        prefix = "unfiltered_"
    end
    if weighting == ""
        suffix=""
    else
        suffix="_"*weighting
    end
    
    return joinpath(test_runs, "vortex_analysis", prefix*cluster*"_"*method*suffix)
end
