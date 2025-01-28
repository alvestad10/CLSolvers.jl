using CLSolvers
using DifferentialEquations, StochasticDiffEq
using Random
using Statistics
using LaTeXStrings

# In this file the true solutions used in the paper are stored as
# dictionares
include("Scro_Data.jl")


###############################################
###############################################
# SETTING UP THE CONTOUR                      #
###############################################
###############################################

## SELECT THE MODEL
p = AHO_Param(1, 1, 24)
#p = NonEquil_AHO_Param(1,1,1,1,1,1,0,1)


######################
## SELECT CONTOUR  ###
######################
#contour = EucledianContour(1.0)
#contour = NonEquilTiltedSchwingerKeldyshContour(1.0,0.0001)
#contour = NonEquilSchwingerKeldyshContour(0.5,0.0)
contour = SchwingerKeldyshContour(1.0, 0.5, 0.0, 0.0)
#contour = SchwingerKeldyshContourRounded1(1.0,0.6,0.5,0.005,0.0025)
#contour = SchwingerKeldyshContourRounded2(1.0,0.6,0.25,0.1,0.025,0.05)
#contour = SchwingerKeldyshContourRounded3(1.0,0.6,0.25,0.4,0.1,0.15,0.02)


#########################
### DISCRETIZE CONTOUR ##
#########################
# Comment out the lines not used below if this file is used as
# a runnable script

# Set total number of points in contour
t_steps = 64


# Let the library choose how the contour should be distributed
# (Not optimal for non standard contours)
discContour = discretizeContour(contour, t_steps)


# script to distribute equally in a tilted contour
begin
    L1 = CLSolvers.L1(contour) / CLSolvers.T(contour)
    L2 = CLSolvers.L2(contour) / CLSolvers.T(contour)

    t_step1 = ceil(Int64, t_steps * L1)
    t_step2 = floor(Int64, t_steps * L2)
    discContour = discretizeContour(contour, t_step1, t_step2, 0)
end

# Possible to specify number of points in each of the contour parts
# Some examples:
#       - Equal number of points in FW, BW and Eucl part.
discContour = discretizeContour(contour, t_steps / 3, t_steps / 3, t_steps / 3)
#       - Twice as many points along the Eucl. part as the FW and BW parts.
discContour = discretizeContour(contour, t_steps / 4, t_steps / 4, 2 * t_steps / 4)
#       - Twice as many points along the FW and BW parts as the Eucl. part.
discContour = discretizeContour(contour, 2t_steps / 5, 2t_steps / 5, t_steps / 5)


# Get usefull parameters for simulation and plotting
begin
    a = getContourDistances(discContour, contour) # Points separation
    CC = getContour(discContour, contour)         # Contour points (Complex numbers)
    tp = CLSolvers.spread_timepoints(discContour, contour) # Parameterization points
end

# Plot contour to be safe
fig = CLSolvers.scatter(CC, legend=false)
#CLSolvers.scatter(a)


################################
## INITIAL CONDITION ###########
################################
begin
    RANDOM_INIT = true

    if RANDOM_INIT
        rng = MersenneTwister(12345)
        y0 = vcat(randn(rng, Float64, t_steps), zeros(Float64, t_steps))
    else
        y0 = zeros(Float64, 2 * t_steps)
    end
end


######################
# RUN THE SIMULATION #
######################
args = (N_Tr=50, θ=0.6, dt=1e-4, tol=5e-2, dtmax=1e-3, tspan=(0.0, 10.0), adaptive=true)
begin
    params = getParamArray(p, a)
    sde_a, sde_b, _ = CLSolvers.get_sde_funcs(p)

    ff = SDEFunction(sde_a, sde_b)
    prob_sde2 = SDEProblem(ff, sde_b, y0, args.tspan, params)

    function prob_func(prob, i, repeat)
        println("Starting ensamble: ", i)
        remake(prob, u0=vcat(randn(rng, Float64, t_steps), zeros(Float64, t_steps)))
    end

    if RANDOM_INIT
        ensemble_prob = EnsembleProblem(prob_sde2, prob_func=prob_func)
    else
        ensemble_prob = EnsembleProblem(prob_sde2)
    end


    @time sol = solve(
        # prob_sde2,
        ensemble_prob,
        ImplicitEM(theta=args.θ, symplectic=false),
        EnsembleThreads(), trajectories=args.N_Tr,
        progress=true, saveat=0.01, save_start=false,
        dtmax=args.dtmax, dt=args.dt,
        adaptive=args.adaptive,
        abstol=args.tol, reltol=args.tol, maxiters=10^6)
end


#######  Check that simulation completed without errors
#######  removes warning trajectories if remove_warning_tr=true
check_errors(sol, remove_warning_tr=true)

##########################
## CALCULATE OBSERVABLES #
##########################
avg_Re, avg_Im, avg2_Re, avg2_Im, corr0t_Re, corr0t_Im, Gpm_Re, Gpm_Im, Gmp_Re, Gmp_Im = get_statistics(sol, thermal=32)
# (thermal argument should correspond to the last point in the real-time axes to use the Gpm and Gmp correlators)

#######################
## Plotting functions #
#######################
# The properties of the plot can be changed in the src/plot_observables.jl file
# ex. you can set the legend position here

####
# Thermal

# <ϕ(0)ϕ(t)> correlators
plot_realtime_corr0t(corr0t_Re, corr0t_Im, CC; true_solution=dfScroSolMinkBeta1)
plot_full_corr0t(corr0t_Re, corr0t_Im, CC; true_solution=dfScroSolMinkBeta1)

# <ϕ_+(0)ϕ_-(t)> and <ϕ_-(0)ϕ_+(t)> correlators
plot_Gpm_Gmp(Gpm_Re, Gpm_Im, Gmp_Re, Gmp_Im, CC; true_solution=dfScroSolMinkBeta1)

# Plot corr0t along the Eucledian part of the contour
plot_eucledian_corr0t(corr0t_Re, corr0t_Im, CC[1:end-1]; true_solution=dfScroSolEuclBeta1)

# plot <x> and <x^2>
plot_avg_avg2(avg_Re, avg_Im, avg2_Re, avg2_Im, tp[1:end-1]; true_solution=dfScroSolMinkBeta1)

#####
# Non-Equilibrium

# <ϕ(0)ϕ(t)>
plot_full_corr0t(corr0t_Re, corr0t_Im, CC; true_solution=dfScroSolMinkNonEquil)

# plot <x> and <x^2>
plot_avg_avg2_nonEquil(avg_Re, avg_Im, avg2_Re, avg2_Im, CC; true_solution=dfScroSolMinkNonEquil)




######################
## Corrections #######
######################
SS_Re, SS2_Re, SS_Im, SS2_Im = get_corrections(sol, params, args)
plot_corrections(SS_Re, SS2_Re, SS_Im, SS2_Im, avg2_Re, avg2_Im, tp; true_solution=dfScroSolMinkBeta1["avg2Re"])




#################################
### SAVE/LOAD CONFIGS AND SETUP #
#################################

filename = "filename1" # Set the filename to the file be saved/loaded

# Save solution + all information needed to continue simulating or plotting observables
begin
    CLSolvers.save_sol(string("results/", filename, ".jld"), sol, p, contour, a, CC, tp)
end

# Load solution etc.
sol, p, contour, a, CC, tp = load_sol(string("results/", filename, ".jld"))
