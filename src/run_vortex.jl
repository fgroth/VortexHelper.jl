using GadgetIO
using Formatting

global test_runs = "/home/moon/fgroth/phd/test_runs/test_collection/test_runs/"
"""
    set_test_runs(dir::String="/home/moon/fgroth/phd/test_runs/test_collection/test_runs/")

Set the global `testruns` variable that is used to locate the output of simulations.
"""
function set_test_runs(dir::String="/home/moon/fgroth/phd/test_runs/test_collection/test_runs/")
    global test_runs = dir
end

"""
    run_vortex(cluster::String, method::String; 
               start_snap_num=100,end_snap_num=145, scale=3, snaps_todo=nothing, 
               vortex_directory::String="vortex-p", # use "vortex-GADGET" for old code version
               # adjust the vortex parameters
               filtering::Bool=false,
               n_snap::Int64=4,
               cells_per_direction::Int64=128,
               n_levels::Int64=9, n_particles_refinement::Int64=8)

Run Vortex for given snapshots of `cluster` and `method`.
"""
function run_vortex(cluster::String, method::String;
                    start_snap_num=100,end_snap_num=145, scale=3, snaps_todo=nothing,
                    vortex_directory::String="vortex-p", # use "vortex-GADGET" for old code version
                    # adjust the vortex parameters
                    filtering::Bool=false,
                    n_snap::Int64=4,
                    cells_per_direction::Int64=128,
                    n_levels::Int64=9, n_particles_refinement::Int64=8)
    
    last_snapnum=zeros(Int64,1)
    for snapnum in end_snap_num:-1:start_snap_num
        if isdir(test_runs*"/out_"*cluster*"_"*method*"/"*"snapdir_"*sprintf1("%03d",snapnum))
            last_snapnum[1] = snapnum
            break
        else
            continue
        end
    end

    println("last snap num is ",last_snapnum[1])
    
    this_dir = pwd()

    # create temporary directory to run vortex.
    tmp_dir = mktempdir("./")
    vortex_exec = if occursin("mfm",method)
        "vortex_mfm"
    else
        "vortex_sph"
    end
    vortex_exec = if filtering
        vortex_exec * "_filtered"
    else
        vortex_exec * "_unfiltered"
    end
    cp(vortex_directory*"/src/", tmp_dir*"/src")
    cd(tmp_dir*"/src/")
    
    symlink(test_runs*"/out_"*cluster*"_"*method*"/","./simulation")

    prefix = if filtering
        "filtered_"
    else
        ""
    end

    try
        mkdir(test_runs*"/vortex_analysis/"*prefix*cluster*"_"*method)
    catch
        println(test_runs*"/vortex_analysis/"*prefix*cluster*"_"*method*" already exists")
        # directory already exists
    end

    run_sh = open("run.sh","w")
    write(run_sh, "#!/bin/bash\n")
    write(run_sh, "\n")
    write(run_sh, "ulimit -s 128000000\n")
    write(run_sh, "ulimit -v 500000000\n")
    write(run_sh, "ulimit -c 0\n")
    write(run_sh, "export OMP_NUM_THREADS=16\n")
    write(run_sh, "export OMP_STACKSIZE=4000m\n")
    write(run_sh, "export OMP_PROC_BIND=true\n")
    write(run_sh, "\n")
    write(run_sh,"./"*vortex_exec*"\n")
    close(run_sh)
    chmod("run.sh",0o700)
    
    snaps_todo = if snaps_todo == nothing
        start_snap_num:last_snapnum[1]
    else
        snaps_todo
    end
    for i_snap in snaps_todo
        println("running ",i_snap)
        snap = test_runs * "/out_"*cluster*"_"*method*"/snapdir_"*sprintf1("%03d",i_snap)*"/snap_"*sprintf1("%03d",i_snap)
        sub = test_runs * "/out_"*cluster*"_"*method*"/groups_"*sprintf1("%03d",i_snap)*"/sub_"*sprintf1("%03d",i_snap)

        halo_positions = read_subfind(sub, "GPOS")
        halo_radii = try
            read_subfind(sub, "RTOP")
        catch
            read_subfind(sub, "RVIR")
        end
        
        first_halo_position = halo_positions[:,1]
        println("first halo position ",first_halo_position)
        
        first_halo_radius = halo_radii[1]
        println("first halo radius ",first_halo_radius)

        par_name = "./vortex.dat"
        this_par = open(par_name,"w")

        write(this_par, "***********************************************************************
*                  VORTEX-GADGET PARAMETERS FILE                      *
***********************************************************************
*       General parameters block                                      *
***********************************************************************
Files: first, last, every, num files per snapshot -------------------->\n")
        # adjust snap number
        write(this_par, sprintf1("%d",i_snap)*","*sprintf1("%d",i_snap)*",1,"*sprintf1("%d",n_snap)*"\n")
        write(this_par, "Cells per direction (NX,NY,NZ) --------------------------------------->\n")
        write(this_par, sprintf("%g",cells_per_direction)*","*sprintf("%g",cells_per_direction)*","*sprintf("%g",cells_per_direction)*"\n")
        write(this_par, "Max box sidelength (in input length units) --------------------------->\n")
        # adjust size
        #size = scale*first_halo_radius
        size = 25e3
        write(this_par, sprintf1("%g",2*size)*"\n")
        write(this_par, "Domain to keep particles (in input length units; x1,x2,y1,y2,z1,z2) -->\n")
        write(this_par, sprintf1("%g",first_halo_position[1]-size)*","*sprintf1("%g",first_halo_position[1]+size)*","*
            sprintf1("%g",first_halo_position[2]-size)*","*sprintf1("%g",first_halo_position[2]+size)*","*
            sprintf1("%g",first_halo_position[3]-size)*","*sprintf1("%g",first_halo_position[3]+size)*"\n")
        write(this_par, "!***********************************************************************\n!*       Output customisation (0=no, 1=yes)                            *\n!***********************************************************************\n")
        write(this_par, "!Gridded data: kernel length, density (mutually exclusive), velocity -->\n0,1,1\n")
        write(this_par, "Gridded results: vcomp, vsol, scalar_pot, vector_pot, div(v), curl(v)->\n1,1,1,1,1,1\n")
        write(this_par, "Particle results: interpolation error, particle-wise results --------->\n1,1\n")
        write(this_par, "Filter: gridded Mach/ABVC, shocked cells, filtering length, vturb ---->\n")
        if filtering
            write(this_par, "1,1,1,1\n")
        else
            write(this_par, "0,0,0,0\n")
        end
        write(this_par, "***********************************************************************
*       Mesh creation parameters                                      *
***********************************************************************
Number of levels ----------------------------------------------------->\n")
        write(this_par, sprintf1("%d",n_levels)*"\n")
        write(this_par, "Number of particles for a cell to be refinable ----------------------->\n")
        write(this_par, sprintf1("%g",n_particles_refinement)*"\n")
        write(this_par, "Minimum size of a refinement patch to be accepted -------------------->
6
Cells not to be refined from the border (base grid) ------------------>
2
***********************************************************************
*       Velocity interpolation parameters                             *
***********************************************************************
Number of neighbours for interpolation ------------------------------->\n")
        if occursin("mfm",method)
            write(this_par, "32\n")
        else
            write(this_par, "295\n")
        end
        write(this_par,"***********************************************************************
*       Poisson solver                                                *
***********************************************************************
SOR presion parameter, SOR max iter, border for AMR patches ---------->
1e-9,1000,2
***********************************************************************
*       Multifiltering                                                *
***********************************************************************
Multiscale filter: apply filter -------------------------------------->\n")
        if filtering
            write(this_par, "1,1,1,1\n")
        else
            write(this_par, "0,0,0,0\n")
        end
write(this_par, "Filtering parameters: tolerance, growing step, max. num. of its. ----->
0.1,1.05,200
***********************************************************************
*       On-the-fly shock detection (for multifiltering)               *
***********************************************************************
Threshold on velocity divergence (negative, input units) ------------->
-1.25
Threshold on artificial bulk viscosity constant ---------------------->
1.
Use particle's MACH field (0=no, 1=yes), Mach threshold -------------->
1,2.0\n")

        # adjust position

        close(this_par)

        mkdir("output_files/")
        run(`./run.sh`)
        mv("output_files/",test_runs*"/vortex_analysis/"*prefix*cluster*"_"*method*"/"*sprintf1("%03d",i_snap),force=true)
    end
    cd(this_dir)
    rm(tmp_dir,force=true, recursive=true)
    
end
