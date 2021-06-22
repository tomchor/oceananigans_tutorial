using Oceananigans
using Oceananigans.Units
using CUDA: has_cuda

Nx, Ny, Nz = 256, 3200, 128

if has_cuda()
    arch = GPU()
else
    arch = CPU()
    Nx = Int(Nx/4)
    Ny = Int(Ny/4)
    Nz = Int(Nz/4)
end

topology = (Periodic, Bounded, Bounded)
grid = RegularRectilinearGrid(size=(Nx, Ny, Nz),
                              x=(0, 200),
                              y=(0, 2000),
                              z=(-100, 0),
                              topology=topology)
println("\n", grid, "\n")


model = IncompressibleModel(architecture = arch,
                grid = grid,
                closure=nothing,
                )
println("\n", model, "\n")



start_time = 1e-9*time_ns()
using Oceanostics: SingleLineProgressMessenger
simulation = Simulation(model, Δt=10seconds,
                        stop_time=10hours,
                        wall_time_limit=23.5hours,
                        iteration_interval=5,
                        progress=SingleLineProgressMessenger(LES=false, initial_wall_time_seconds=start_time),
                        stop_iteration=Inf,)

println("\n", simulation,"\n",)
@info "---> Starting run!\n"
run!(simulation, pickup=false)

