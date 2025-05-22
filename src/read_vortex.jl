using FortranFiles

"""
    read_particle_data(file::String)

Read particle data from unformatted fortran output of vortex.
Return `Dict` containing `"pos"`, `"vel_orig"`, `"mass"`, `"vel"`, `"velcomp"`, `"velrot"` vectors.
"""
function read_particle_data(file::String)

    f = FortranFile(file) 
    n = read(f,Int32)

    data=Dict()
    
    # read positions
    data["pos"] = zeros(Float32,3,n)
    data["pos"][1,:] = read(f,data["pos"][1,:])
    data["pos"][2,:] = read(f,data["pos"][2,:])
    data["pos"][3,:] = read(f,data["pos"][3,:])

    # read velocities
    data["vel_orig"] = zeros(Float32,3,n)
    data["vel_orig"][1,:] = read(f,data["vel_orig"][1,:])
    data["vel_orig"][2,:] = read(f,data["vel_orig"][2,:])
    data["vel_orig"][3,:] = read(f,data["vel_orig"][3,:])
    
    # read mass
    data["mass"] = zeros(Float32,n)
    data["mass"] = read(f,data["mass"])
    
    # read velocities
    data["vel"] = zeros(Float32,3,n)
    data["vel"][1,:] = read(f,data["vel"][1,:])
    data["vel"][2,:] = read(f,data["vel"][2,:])
    data["vel"][3,:] = read(f,data["vel"][3,:])
    
    # read compressive velocities
    data["velcomp"] = zeros(Float32,3,n)
    data["velcomp"][1,:] = read(f,data["velcomp"][1,:])
    data["velcomp"][2,:] = read(f,data["velcomp"][2,:])
    data["velcomp"][3,:] = read(f,data["velcomp"][3,:])

    # read solenoidal velocities
    data["velrot"] = zeros(Float32,3,n)    
    data["velrot"][1,:] = read(f,data["velrot"][1,:])
    data["velrot"][2,:] = read(f,data["velrot"][2,:])
    data["velrot"][3,:] = read(f,data["velrot"][3,:])

    close(f)
    
    return data
end
