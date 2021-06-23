using Oceananigans

# Create grid
grid = RegularRectilinearGrid(size=(128, 128), extent=(2π, 2π), 
                              topology=(Periodic, Periodic, Flat))

# Create model
model = IncompressibleModel(architecture=CPU(),
                            timestepper = :RungeKutta3,
                              advection = UpwindBiasedFifthOrder(),
                                   grid = grid,
                               buoyancy = nothing,
                                tracers = nothing,
                                closure = IsotropicDiffusivity(ν=1e-5)
                           )

# Set IC
using Statistics
u₀ = rand(size(model.grid)...)
u₀ .-= mean(u₀)
set!(model, u=u₀, v=u₀)


# Create Simulation
progress(sim) = @info "Iteration: $(sim.model.clock.iteration), time: $(round(Int, sim.model.clock.time))"
simulation = Simulation(model, Δt=0.2, stop_time=50, iteration_interval=100, progress=progress)


# Create Diagnostics
u, v, w = model.velocities

ω = ∂x(v) - ∂y(u)
ω_field = ComputedField(ω)

s = sqrt(u^2 + v^2)
s_field = ComputedField(s)


# Set-up output
simulation.output_writers[:fields] = NetCDFOutputWriter(model, (u=u, v=v, ω=ω_field, s=s_field),
                                                      schedule = TimeInterval(2),
                                                      filepath = "two_dimensional_turbulence.nc",
                                                      mode = "c")

# Run!
run!(simulation)


