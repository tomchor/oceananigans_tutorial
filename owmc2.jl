using Random
using Printf

using Oceananigans
using Oceananigans.Units



# Create grid
arch = CPU()
Lz = 32
Nx = Ny = Nz = 32
S = 0.8 # Stretching factor
hyperbolically_spaced_nodes(k) = -Lz-Lz*(tanh(S * ( (-(k-34) - 1) / Nz - 1)) / tanh(S))
grid = VerticallyStretchedRectilinearGrid(size = (Nx, Ny, Nz), 
                                          architecture = arch,
                                          x = (0,64),
                                          y = (0,64),
                                          halo = (3, 3, 3),
                                          z_faces = hyperbolically_spaced_nodes)
println(grid, "\n\n")



# Create prelim variables for Model (e.g. BC and Buoyancy model)
buoyancy = Buoyancy(model=BuoyancyTracer())
Qᵇ = 1.5e-7 # m/s² m/s
dbdz = 1.6e-5
b_sfc_flux = FluxBoundaryCondition(Qᵇ)
b_bot_grad = GradientBoundaryCondition(dbdz)
b_bcs = TracerBoundaryConditions(grid, top=b_sfc_flux, bottom=b_bot_grad)

u₁₀ = 10     # m s⁻¹, average wind velocity 10 meters above the ocean
cᴰ = 2.5e-3  # dimensionless drag coefficient
ρₒ = 1026 # kg m⁻³, average density at the surface of the world ocean
ρₐ = 1.225   # kg m⁻³, average density of air at sea-level
Qᵘ = - ρₐ / ρₒ * cᴰ * u₁₀ * abs(u₁₀) # m² s⁻²
u_bcs = UVelocityBoundaryConditions(grid, top = FluxBoundaryCondition(Qᵘ))




# Create model
using Oceananigans.TurbulenceClosures: SmagorinskyLilly
model = IncompressibleModel(architecture = arch,
                            advection = WENO5(),
                            timestepper = :RungeKutta3,
                            grid = grid,
                            tracers = (:b,), # Important!
                            coriolis = FPlane(f=1e-4),
                            buoyancy = buoyancy,
                            closure = SmagorinskyLilly(C=0.13),
                            boundary_conditions = (u=u_bcs, b=b_bcs,))
println(model, "\n\n")


# Impose ICs
Ξ(z) = randn() * z / model.grid.Lz * (1 + z / model.grid.Lz) # noise
bᵢ(x, y, z) = dbdz * (z+grid.Lz) + 1e-8 * Ξ(z)
uᵢ(x, y, z) = sqrt(abs(Qᵘ)) * 1e-4 * Ξ(z)
set!(model, u=uᵢ, w=uᵢ, b=bᵢ)


# Create Simulation
wizard = TimeStepWizard(cfl=1.0, Δt=10.0, max_change=1.1, max_Δt=1minute, min_Δt=0.1second)
start_time = 1e-9*time_ns() # so we can print the total elapsed wall time
using Oceanostics: SingleLineProgressMessenger
simulation = Simulation(model, Δt=wizard, stop_time=40minutes, iteration_interval=10,
                        progress=SingleLineProgressMessenger(LES=true, initial_wall_time_seconds=start_time))
println(simulation, "\n\n")


# Create diagnostics
u, v, w = model.velocities
b = model.tracers.b
νₑ = model.diffusivities.νₑ
ζ = ∂x(v) - ∂y(u) # Vertical vorticity
δ = ∂x(u) + ∂y(v) # Horizontal convergence
ε = ComputedField(νₑ*(∂x(u)^2 + ∂x(v)^2 + ∂x(w)^2 + # Pseudo energy dissipation rate
                      ∂y(u)^2 + ∂y(v)^2 + ∂y(w)^2 + 
                      ∂z(u)^2 + ∂z(v)^2 + ∂z(w)^2))

XYAverage(field) = AveragedField(field, dims=(1,2))
ū = XYAverage(u)
v̄ = XYAverage(v)

using Oceanostics: TurbulentKineticEnergy
tke = TurbulentKineticEnergy(model, u, v, w; U=ū, V=v̄, W=0)


operations = (; ζ, δ, ε, tke,)
computed_fields = map(ComputedField, operations)
averaged_fields = map(XYAverage, operations)

fields = (; u, v, w, b, νₑ)

# Set-up outputs
output_slices = merge(fields, computed_fields)
simulation.output_writers[:slices] =
    NetCDFOutputWriter(model, output_slices,
                           filepath = "ocean_wind_mixing_and_convection2.nc",
                     field_slicer = FieldSlicer(j=Int(grid.Ny/2)),
                         schedule = TimeInterval(1minute),
                            mode = "c")

output_average = merge(fields, averaged_fields)
simulation.output_writers[:averages] =
    NetCDFOutputWriter(model, output_average,
                           filepath = "ocean_wind_mixing_and_convection2_avg.nc",
                         schedule = TimeInterval(20seconds),
                            mode = "c")

# Run!
@info "Running simulation"
run!(simulation)

