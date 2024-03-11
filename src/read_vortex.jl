using FortranFiles

"""
    read_particle_data(file::String)

Read particle data from unformatted fortran output of vortex.
Return `pos`, `mass`, `vel`, `velcomp`, `velrot` vectors.
"""
function read_particle_data(file::String)

    f = FortranFile(file) 
    n = read(f,Int32)

    # read positions
    pos = zeros(Float32,3,n)
    pos[1,:] = read(f,pos[1,:])
    pos[2,:] = read(f,pos[2,:])
    pos[3,:] = read(f,pos[3,:])

    # read mass
    mass = zeros(Float32,n)
    mass = read(f,mass)

    # read velocities
    vel = zeros(Float32,3,n)
    vel[1,:] = read(f,vel[1,:])
    vel[2,:] = read(f,vel[2,:])
    vel[3,:] = read(f,vel[3,:])
    
    # read compressive velocities
    velcomp = zeros(Float32,3,n)
    velcomp[1,:] = read(f,velcomp[1,:])
    velcomp[2,:] = read(f,velcomp[2,:])
    velcomp[3,:] = read(f,velcomp[3,:])

    # read solenoidal velocities
    velrot = zeros(Float32,3,n)    
    velrot[1,:] = read(f,velrot[1,:])
    velrot[2,:] = read(f,velrot[2,:])
    velrot[3,:] = read(f,velrot[3,:])

    return pos, mass, vel, velcomp, velrot

end
