# CLSolvers.jl

This is the code generating the simulations done in: https://arxiv.org/abs/2105.02735

To run the code you need a version of Julia installed, then you can make separate scripts or
follow the Main.jl file which can be run line by line using the Julia vscode extension.

## Instantiate

To initialize the project run these comments inside the Julia REPL (From inside the project directory)
```julia
    import Pkg
    Pkg.activate(".")
    Pkg.instantiate()
```
For more information see: https://docs.julialang.org/en/v1/stdlib/Pkg/

Now all dependencies should be downloaded and the code is ready to be run.

## Example

Under is an example simulating the Anharmonic oscillator on the canonical Schwinger-Keldysh contour
and plotting the two-point function compared to the true solution.

```julia
using CLSolvers
using DifferentialEquations, StochasticDiffEq
using Random
rng = MersenneTwister(12345);

# In this file the true solutions used in the paper are stored as
# dictionares
include("Scro_Data.jl")

# Setup the contour
t_steps = 64
contour = SchwingerKeldyshContour(1.0,0.5,0.0,0.0)
discContour = discretizeContour(contour,t_steps/4,t_steps/4,2*t_steps/4)

# Get out the needed information from the contour
a = getContourDistances(discContour, contour) # Points separation
CC = getContour(discContour, contour)         # Contour points (Complex numbers)
tp = CLSolvers.spread_timepoints(discContour,contour) # Parameterization points

# initialize the field
y0 = vcat(randn(rng, Float64, t_steps),zeros(Float64,t_steps))

# Setup parameters for the simulation
args = (N_Tr = 50, # Number of trajectories 
          θ = 0.6, # Impliciteness of solver
          dt = 1e-4, tol = 5e-2, dtmax = 1e-3, #Langevin stepsize paramater
          tspan=(0.0,10.0), # Langevin time span
          adaptive=true) # Select adaptive step-size or fixed step-size

# Get parameters and functions on a form that the StochasticDiffEq can read
params = getParamArray(p,a)
sde_a, sde_b, _ = CLSolvers.get_sde_funcs(p)

ff = SDEFunction(sde_a,sde_b)
prob_sde2 = SDEProblem(ff,sde_b,y0,args.tspan,params)

# Funciton that is run before each trajectory (Initialize the trajectory)
function prob_func(prob,i,repeat)
    println("Starting ensamble: ",i)
    remake(prob,u0=vcat(randn(rng, Float64, t_steps),zeros(Float64,t_steps)))
end

# Setup a ensamble simulation (See https://diffeq.sciml.ai/stable/tutorials/sde_example/#Ensemble-Simulations)
ensemble_prob = EnsembleProblem(prob_sde2,prob_func=prob_func)

# Run the simulation
@time sol = solve(ensemble_prob,ImplicitEM(theta=args.θ, symplectic=false),
            EnsembleThreads(), trajectories=args.N_Tr,
            progress = true, saveat = 0.01, save_start = false, 
            dtmax=args.dtmax, dt=args.dt,
            adaptive=args.adaptive,
            abstol=args.tol,reltol=args.tol,maxiter=10^6)

# Compute the observables
avg_Re, avg_Im, avg2_Re, avg2_Im, corr0t_Re, corr0t_Im, corr0t_2_Re, corr0t_2_Im, corrt0_Re, corrt0_Im = get_statistics(sol)

# Plot the two-point function
plot_realtime_corr0t(corr0t_Re, corr0t_Im, CC; true_solution=dfScroSolMinkBeta1)
```