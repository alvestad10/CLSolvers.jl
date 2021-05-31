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
p = AHO_Param(1,1,24)
#p = NonEquil_AHO_Param(1,1,1,1,1,1,0,1)


######################
## SELECT CONTOUR  ###
######################
#contour = EucledianContour(1.0)
#contour = NonEquilTiltedSchwingerKeldyshContour(1.0,0.0001)
#contour = NonEquilSchwingerKeldyshContour(0.5,0.0)
contour = SchwingerKeldyshContour(1.0,0.5,0.0,0.0)
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
discContour = discretizeContour(contour,t_steps)


# script to distribute equally in a tilted contour
begin
    L1 = CLSolvers.L1(contour)/CLSolvers.T(contour)
    L2 = CLSolvers.L2(contour)/CLSolvers.T(contour)

    t_step1 = ceil(Int64,t_steps*L1)
    t_step2 = floor(Int64,t_steps*L2)
    discContour = discretizeContour(contour,t_step1,t_step2,0)
end

# Possible to specify number of points in each of the contour parts
# Some examples:
#       - Equal number of points in FW, BW and Eucl part.
discContour = discretizeContour(contour,t_steps/3,t_steps/3,t_steps/3) 
#       - Twice as many points along the Eucl. part as the FW and BW parts.
discContour = discretizeContour(contour,t_steps/4,t_steps/4,2*t_steps/4)
#       - Twice as many points along the FW and BW parts as the Eucl. part.
discContour = discretizeContour(contour,2t_steps/5,2t_steps/5,t_steps/5)


# Get usefull parameters for simulation and plotting
begin
    a = getContourDistances(discContour, contour)
    CC = getContour(discContour, contour)
    tp = CLSolvers.spread_timepoints(discContour,contour)
end

# Plot contour to be safe
fig = CLSolvers.scatter(CC,legend=false,xlim=[0,1])
#CLSolvers.scatter(a)


################################
## INITIAL CONDITION ###########
################################
begin
RANDOM_INIT = true

if RANDOM_INIT
    rng = MersenneTwister(12345);
    y0 = vcat(randn(rng, Float64, t_steps),zeros(Float64,t_steps))
else
    y0 = zeros(Float64,2*t_steps)
end
end


######################
# RUN THE SIMULATION #
######################
args = (N_Tr = 50, θ = 0.6, dt = 1e-4, tol = 5e-2, dtmax = 1e-3, tspan=(0.0,10.0), adaptive=true)
begin
params = getParamArray(p,a)
sde_a, sde_b, _ = CLSolvers.get_sde_funcs(p)

ff = SDEFunction(sde_a,sde_b)
prob_sde2 = SDEProblem(ff,sde_b,y0,args.tspan,params)

function prob_func(prob,i,repeat)
    println("Starting ensamble: ",i)
    remake(prob,u0=vcat(randn(rng, Float64, t_steps),zeros(Float64,t_steps)))
end

if RANDOM_INIT
    ensemble_prob = EnsembleProblem(prob_sde2,prob_func=prob_func)
else
    ensemble_prob = EnsembleProblem(prob_sde2)
end


@time sol = solve(ensemble_prob,ImplicitEM(theta=args.θ, symplectic=false),
            EnsembleThreads(), trajectories=args.N_Tr,
            progress = true, saveat = 0.01, save_start = false, 
            dtmax=args.dtmax, dt=args.dt,
            adaptive=args.adaptive,
            abstol=args.tol,reltol=args.tol,maxiter=BigInt(10)^18)

end


#######  Check that simulation completed without errors
#######  removes warning trajectories if remove_warning_tr=true
check_errors(sol,remove_warning_tr=true)

##########################
## CALCULATE OBSERVABLES #
##########################
avg_Re, avg_Im, avg2_Re, avg2_Im, corr0t_Re, corr0t_Im, corr0t_2_Re, corr0t_2_Im, corrt0_Re, corrt0_Im = get_statistics(sol)


#######################
## Plotting functions #
#######################
plot_realtime_corr0t(corr0t_Re, corr0t_Im, CC; true_solution=dfScroSolMinkBeta1)

plot_full_corr0t(corr0t_Re, corr0t_Im, CC; true_solution=dfScroSolMinkBeta1, save=false, fig_name="thermal_corrx0xt")
plot_full_corr0t_corrt0(corr0t_Re, corr0t_Im, corr0t_2_Re, corr0t_2_Im, corrt0_Re, corrt0_Im, CC; true_solution=dfScroSolMinkNonEquil, save=false,fig_name="thermal_G+-_G-+")

plot_full_avg_realtime(avg_Re, avg_Im, CC; true_solution=dfScroSolMinkNonEquil)
plot_eucledian_corr0t(corr0t_Re, corr0t_Im, CC[1:end-1]; true_solution=dfScroSolEuclBeta1, save=true, fig_name="eucledian_corrx0xt_noTilt")
plot_eucledian_avg(avg_Re,avg_Im,avg2_Re,avg2_Im,tp[1:end-1]; true_solution=dfScroSolMinkBeta1,save=true,fig_name="tparam_x_x2_thermal")
plot_eucledian_avg_nonEquil(avg_Re,avg_Im,avg2_Re,avg2_Im,CC[1:end]; true_solution=dfScroSolMinkNonEquil,save=false,fig_name="tparam_x_x2_nonEquil")




######################
## Corrections #######
######################
SS_Re, SS2_Re, SS_Im, SS2_Im = get_corrections(sol,params,args)
plot_corrections(SS_Re, SS2_Re, SS_Im, SS2_Im, avg2_Re, avg2_Im, tp; true_solution=dfScroSolMinkBeta1["avg2Re"])




#################################
### SAVE/LOAD CONFIGS AND SETUP #
#################################

filename = "filename1" # Set the filename to the file be saved/loaded

# Save solution + all information needed to continue simulating or plotting observables
begin
    CLSolvers.save_sol(string("results/",filename,".jld"),sol, p, contour, a, CC, tp)
end

# Load solution etc.
sol, p, contour, a, CC, tp = load_sol(string("results/",filename,".jld"))
