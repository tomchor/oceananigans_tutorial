using Pkg
Pkg.instantiate()
using ArgParse
using Printf
using Oceananigans
using Oceananigans.Units
using Oceananigans.Advection: WENO5
using Oceananigans.OutputWriters, Oceananigans.Fields
using SpecialFunctions: erf
using CUDA: has_cuda


# Read and parse initial arguments
#++++
"Returns a dictionary of command line arguments."
function parse_command_line_arguments()
    settings = ArgParseSettings()
    @add_arg_table! settings begin

        "--factor"
            help = "Factor to divide Nh and Nz for"
            default = 32
            arg_type = Int

        "--arch"
            help = "CPU or GPU"
            default = "CPU"
            arg_type = String

        "--fullname"
            help = "Setup and name of jet in jetinfo.jl"
            default = "S2d_CIjet1"
            arg_type = String
    end
    return parse_args(settings)
end
args = parse_command_line_arguments()
factor = args["factor"]
arch = eval(Meta.parse(args["arch"]*"()"))
fullname = args["fullname"]

setup, jet = split(fullname, "_")
ndims = parse(Int, strip(setup, ['S', 'd']))
jet = Symbol(jet)

if ndims ∉ [2,3] # Validade input
    throw(AssertionError("Dimensions must be 2 or 3"))
end

@printf("Starting %sd jet %s with a dividing factor of %d and a %s architecture\n", 
        ndims,
        jet,
        factor,
        arch,)
#-----


# Get simulation parameters
#++++
as_background=false
include("jetinfo.jl")

if ndims==3
    simulation_nml = getproperty(SurfaceJetSimulations(Ny=400*2^4, Nz=2^7), jet)
    prefix = "PNN"
else
    simulation_nml = getproperty(SurfaceJetSimulations(), jet)
    prefix = "FNN"
end
LES = false
@unpack name, f0, u₀, N2_inf, N2_pyc, Ny, Nz, Ly, Lz, σy, σz, y₀, z₀, νz, sponge_frac = simulation_nml

simname = "$(prefix)_$(name)"
pickup = any(startswith("chk.$simname"), readdir("data"))
#-----

# Set GRID
#++++  GRID
if ndims==3
    Nx = Ny÷64
    Lx = (Ly / Ny) * Nx
else
    Nx = factor
    Lx = 6 * (Ly / Ny) * Nx
end
topology = (Periodic, Bounded, Bounded)

grid = RegularRectilinearGrid(size=(Nx÷factor, Ny÷factor, Nz÷factor),
                              x=(0, Lx),
                              y=(0, Ly),
                              z=(-Lz, 0), 
                              topology=topology)
println("\n", grid, "\n")
#-----


# Calculate secondary parameters
#++++
b₀ = u₀ * f0
ρ₀ = 1027
T_inertial = 2*π/f0
y_r = y₀ + √2/4 * σy
z_r = 0
Ro_r = - √2 * u₀ * (z₀/σz-1) * exp(-1/8) / (2*f0*σy)
Ri_r = N2_inf * σz^2 * exp(1/4) / u₀^2
νh = νz * (grid.Δy / grid.Δz)^(4/3)

secondary_params = merge((LES=Int(LES), u_0=u₀, y_0=y₀, z_0=z₀, b0=b₀), 
                         (;y_r, z_r, Ro_r, Ri_r, T_inertial, νh))

global_attributes = merge(simulation_nml, secondary_params)
println("\n", global_attributes, "\n")
#-----



# Set up Geostrophic flow
#++++++
const n2_inf = N2_inf
const n2_pyc = N2_pyc
const Hz = grid.Lz
const Hy = grid.Ly
const sig_z = σz
const sig_y = σy
const u_0 = u₀
const y_0 = y₀
const z_0 = z₀
const z_c = -40
const z_m = z_c - n2_pyc/n2_inf*(z_c+Hz)
const f_0 = f0
#-----


# Define model!
#++++
if LES
    import Oceananigans.TurbulenceClosures: SmagorinskyLilly, AnisotropicMinimumDissipation
    closure = SmagorinskyLilly(C=0.16)
    #closure = AnisotropicMinimumDissipation()
else
    import Oceananigans.TurbulenceClosures: AnisotropicDiffusivity, IsotropicDiffusivity
    closure = AnisotropicDiffusivity(νh=νh, κh=νh, νz=νz, κz=νz)
end
model_kwargs = (architecture = arch,
                grid = grid,
                advection = WENO5(),
                timestepper = :RungeKutta3,
                coriolis = FPlane(f=f0),
                tracers = (:b,),
                buoyancy = BuoyancyTracer(),
                )
model = IncompressibleModel(; model_kwargs..., closure=nothing)
println("\n", model, "\n")
#-----


# Finally define Simulation!
#++++
start_time = 1e-9*time_ns()
using Oceanostics: SingleLineProgressMessenger
simulation = Simulation(model, Δt=10seconds, 
                        stop_time=10*T_inertial,
                        wall_time_limit=23.5hours,
                        iteration_interval=5,
                        progress=SingleLineProgressMessenger(LES=LES, initial_wall_time_seconds=start_time),
                        stop_iteration=Inf,)
println("\n", simulation, "\n")
#-----



# Run the simulation!
#+++++
println("\n", simulation,
        "\n",)

@printf("---> Starting run!\n")
run!(simulation, pickup=false)
#-----
