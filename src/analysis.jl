export autocor_sol, autocor_sol_corr0t, thermalization, thermalization_corr0t, get_statistics,
        distribution_avg,distribution_avg2,distribution_corr0t, get_corrections

##################################################

using StatsBase,Statistics,Measurements,Bootstrap

"""
plot AUTO-CORRELATION from solution of DifferentialEquation.jl solver
"""
function autocor_sol(sol,TT;therm=1, trajectories=[1],inxs::Array{Int64}=[1])
    fig = plot(legend=false)
    for tr in trajectories
        D = transpose(hcat(sol[tr].u[therm:TT:end]...))
        plot!(fig,autocor(D),label=string("Tr: ",tr))
    end
    display(fig)
end

function autocor_sol_corr0t(sol,TT;therm=1, trajectories=[1])
    t_steps = div(size(sol[1])[1],2)
    fig = plot(legend=false)
    for tr in trajectories
        SMu = hcat(sol[tr].u[therm:TT:end]...)
        corr0t_Re_1 = SMu[1,:] .* SMu[1,:] .- SMu[t_steps+1,:] .* SMu[t_steps+1,:]

        plot!(fig,mean(autocor(corr0t_Re_1),dims=2))
    end
    display(fig)
end

   
function thermalization(sol,TT::Integer,therm::Integer; trajectory=1)
    thermalization(sol,TT,therm,collect(1:size(sol[1])[1]); trajectory=trajectory)
end

function thermalization(sol,TT::Integer,therm::Integer, inxs::Array{Int64}; trajectory=1)
    D = hcat(sol[trajectory].u[therm:TT:end]...)[inxs,:]
    println("Size of solution after thermalization: ", size(D))
    plot(transpose(D),legend=true)#,xlim=[0,100])
end

function thermalization_corr0t(sol,TT::Integer,therm::Integer,inxs::Array{Int64})
    t_steps = div(size(sol[1])[1],2)
    SMu = hcat(sol[trajectory].u[therm:TT:end]...)
    corr0t_Re_1 = SMu[1,:] .* SMu[1,:] .- SMu[t_steps+1,:] .* SMu[t_steps+1,:]
    corr0t_Re_2 = SMu[1,:] .* SMu[4,:] .- SMu[t_steps + 1 ,:] .* SMu[t_steps+4,:]
    corr0t_Re_3 = SMu[1,:] .* SMu[floor(Int64,t_steps/3),:] .- SMu[t_steps + 1 ,:] .* SMu[floor(Int64,4*t_steps/3),:]
    corr0t_Re_4 = SMu[1,:] .* SMu[floor(Int64,t_steps/2),:] .- SMu[t_steps + 1 ,:] .* SMu[floor(Int64,3*t_steps/2),:]
    plot(corr0t_Re_1)#,legend=false)#,xlim=[0,10])
    plot!(corr0t_Re_2)
    plot!(corr0t_Re_3)
    plot!(corr0t_Re_4)#,xlim=[0,100])
end

function distribution_avg(sol,TT::Integer,therm::Integer, inxs::Array{Int64})
    N_Tr = size(sol)[end]

    D = hcat(sol[1].u[therm:TT:end]...)[inxs,:]
    for tr in 2:N_Tr
        D = hcat(D,hcat(sol[tr].u[therm:TT:end]...)[inxs,:])
    end
    histogram(transpose(D),bins=100,normalize=true)
end

function distribution_avg2(sol,TT::Integer,therm::Integer, inxs::Array{Int64}; Re=true, Im=false)
    t_steps = div(size(sol[1])[1],2)
    N_Tr = size(sol)[end]

    D = hcat(sol[1].u[therm:TT:end]...)[inxs,:]
    D2 =  hcat(sol[1].u[therm:TT:end]...)[t_steps .+ inxs,:]
    for tr in 2:N_Tr
        D  = hcat(D, hcat(sol[tr].u[therm:TT:end]...)[inxs,:])
        D2 = hcat(D2,hcat(sol[tr].u[therm:TT:end]...)[t_steps .+ inxs,:])
    end

    if Re
        V = transpose(D.^2 .+ D2.^2)
    elseif Im
        V = transpose(2*D .* D2)
    end

    bin_size = 0.1
    bins = floor(Int64,(maximum(V)-minimum(V))/bin_size)

    #histogram(V,bins=bins,normalize=true, xlim=[-3,3], seriestype=:stephist)
    StatsPlots.density(V,xlim=[-2,2])
end

function distribution_corr0t(sol,TT::Integer,therm::Integer, inxs::Array{Int64})#, inxs::Array{Int64})
    t_steps = div(size(sol[1])[1],2)
    N_Tr = size(sol)[end]

    D = hcat(sol[1].u[therm:TT:end]...)[inxs,:]
    D2 =  hcat(sol[1].u[therm:TT:end]...)[t_steps .+ inxs,:]
    for tr in 2:N_Tr
        D  = hcat(D, hcat(sol[tr].u[therm:TT:end]...)[inxs,:])
        D2 = hcat(D2,hcat(sol[tr].u[therm:TT:end]...)[t_steps .+ inxs,:])
    end
    histogram(D,bins=70,normalize=true,legend=false)
end
    




## Statistics
function get_statistics(sol;TT=1,therm=1,thermal = 0)

    N_Tr = length(sol)
    t_steps = div(size(sol[1])[1],2)
    n_boot  = 50

    avg_Re = zeros(Float64,N_Tr,t_steps)
    avg_Im = zeros(Float64,N_Tr,t_steps)
    avg2_Re = zeros(Float64,N_Tr,t_steps)
    avg2_Im = zeros(Float64,N_Tr,t_steps)
    corr0t_Re = zeros(Float64,N_Tr,t_steps)
    corr0t_Im = zeros(Float64,N_Tr,t_steps)
    Gpm_Re = zeros(Float64,N_Tr,t_steps)
    Gpm_Im = zeros(Float64,N_Tr,t_steps)
    Gmp_Re = zeros(Float64,N_Tr,t_steps)
    Gmp_Im = zeros(Float64,N_Tr,t_steps)

    

    Threads.@threads for tr in 1:N_Tr
    #for tr in 1:N_Tr
        SMu = hcat(sol[tr].u[therm:TT:end]...)
        
        for i in 1:t_steps-thermal
            avg_Re[tr,i] = mean(SMu[i,:])        
            avg_Im[tr,i] = mean(SMu[t_steps+i,:])
            
            avg2_Re[tr,i] = mean(SMu[i,:].^2 - SMu[t_steps+i,:].^2)
            avg2_Im[tr,i] = mean(2*SMu[i,:] .* SMu[t_steps+i,:])
            
            corr0t_Re_tmp = SMu[1,:] .* SMu[i,:] .- SMu[t_steps+1,:] .* SMu[t_steps+i,:] 
            corr0t_Re[tr,i] = mean(corr0t_Re_tmp)
        
            corr0t_Im_tmp =  SMu[1,:] .* SMu[t_steps+i,:] .+ SMu[i,:] .* SMu[t_steps+1,:] 
            corr0t_Im[tr,i] = mean(corr0t_Im_tmp)

            Gpm_Re_tmp = SMu[1,:] .* SMu[t_steps-i-thermal+1,:] .- SMu[t_steps+1,:] .* SMu[2*t_steps-i-thermal+1,:] 
            Gpm_Re[tr,i] = mean(Gpm_Re_tmp)
        
            Gpm_Im_tmp =  SMu[1,:] .* SMu[2t_steps-i-thermal+1,:] .+ SMu[t_steps-i-thermal+1,:] .* SMu[t_steps+1,:] 
            Gpm_Im[tr,i] = mean(Gpm_Im_tmp)

            Gmp_Re_tmp = SMu[t_steps-thermal,:] .* SMu[i,:] .- SMu[2*t_steps-thermal,:] .* SMu[t_steps+i,:] 
            Gmp_Re[tr,i] = mean(Gmp_Re_tmp)

            Gmp_Im_tmp =  SMu[t_steps-thermal,:] .* SMu[t_steps+i,:] .+ SMu[i,:] .* SMu[2*t_steps-thermal,:] 
            Gmp_Im[tr,i] = mean(Gmp_Im_tmp)

        end
    end

    d = 1
    
    mean_d(x) = mean(x,dims=d) 
    avg_Re_bts = bootstrap(mean_d, avg_Re, BalancedSampling(n_boot))
    avg_Im_bts = bootstrap(mean_d, avg_Im, BalancedSampling(n_boot))
    avg2_Re_bts = bootstrap(mean_d, avg2_Re, BalancedSampling(n_boot))
    avg2_Im_bts = bootstrap(mean_d, avg2_Im, BalancedSampling(n_boot))
    corr0t_Re_bts = bootstrap(mean_d, corr0t_Re, BalancedSampling(n_boot))
    corr0t_Im_bts = bootstrap(mean_d, corr0t_Im, BalancedSampling(n_boot))
    Gpm_Re_bts = bootstrap(mean_d, Gpm_Re, BalancedSampling(n_boot))
    Gpm_Im_bts = bootstrap(mean_d, Gpm_Im, BalancedSampling(n_boot))
    Gmp_Re_bts = bootstrap(mean_d, Gmp_Re, BalancedSampling(n_boot))
    Gmp_Im_bts = bootstrap(mean_d, Gmp_Im, BalancedSampling(n_boot))

    std_err(x) = std(x,dims=d,corrected=true) ./ sqrt(size(x)[2])
    avg_Re_bts_std = bootstrap(std_err, avg_Re, BalancedSampling(n_boot))
    avg_Im_bts_std = bootstrap(std_err, avg_Im, BalancedSampling(n_boot))
    avg2_Re_bts_std = bootstrap(std_err, avg2_Re, BasicSampling(n_boot))
    avg2_Im_bts_std = bootstrap(std_err, avg2_Im, BalancedSampling(n_boot))
    corr0t_Re_bts_std = bootstrap(std_err, corr0t_Re, BalancedSampling(n_boot))
    corr0t_Im_bts_std = bootstrap(std_err, corr0t_Im, BalancedSampling(n_boot))
    Gpm_Re_bts_std = bootstrap(std_err, Gpm_Re, BalancedSampling(n_boot))
    Gpm_Im_bts_std = bootstrap(std_err, Gpm_Im, BalancedSampling(n_boot))
    Gmp_Re_bts_std = bootstrap(std_err, Gmp_Re, BalancedSampling(n_boot))
    Gmp_Im_bts_std = bootstrap(std_err, Gmp_Im, BalancedSampling(n_boot))

    avg_Re = collect(original(avg_Re_bts)) .± original(avg_Re_bts_std)
    avg_Im = collect(original(avg_Im_bts)) .± original(avg_Im_bts_std)
    avg2_Re = collect(original(avg2_Re_bts)) .± (original(avg2_Re_bts_std))
    avg2_Im = collect(original(avg2_Im_bts)) .± original(avg2_Im_bts_std)
    corr0t_Re = (collect(original(corr0t_Re_bts)) .± original(corr0t_Re_bts_std)) .- avg_Re.*avg_Re[1]
    corr0t_Im = (collect(original(corr0t_Im_bts)) .± original(corr0t_Im_bts_std)) .- avg_Im.*avg_Im[1]
    Gpm_Re = (collect(original(Gpm_Re_bts)) .± original(Gpm_Re_bts_std)) .- avg_Re.*avg_Re[1]
    Gpm_Im = (collect(original(Gpm_Im_bts)) .± original(Gpm_Im_bts_std)) .- avg_Im.*avg_Im[1]
    Gmp_Re = (collect(original(Gmp_Re_bts)) .± original(Gmp_Re_bts_std)) .- avg_Re.*avg_Re[t_steps]
    Gmp_Im = (collect(original(Gmp_Im_bts)) .± original(Gmp_Im_bts_std)) .- avg_Im.*avg_Im[t_steps]
    
    return avg_Re, avg_Im, avg2_Re, avg2_Im, corr0t_Re, corr0t_Im, Gpm_Re, Gpm_Im, Gmp_Re, Gmp_Im  #twoPoint_Re, twoPoint_Im
    
end

function get_corrections(sol,params,args)

    N_Tr = length(sol)
    t_steps = div(size(sol[1])[1],2)

    du = zeros(Float64,t_steps*2)
    SS_Re_tmp = zeros(Float64,length(sol),t_steps)
    SS2_Re_tmp = zeros(Float64,length(sol),t_steps)
    SS_Im_tmp = zeros(Float64,length(sol),t_steps)
    SS2_Im_tmp = zeros(Float64,length(sol),t_steps)
    for tr in 1:N_Tr
        for i in 1:length(sol[tr][1,:])
            a_AHO_real(du,sol[tr][:,i],params,1)
            SS_Re_tmp[tr,:]  += du[1:t_steps].^2 .- du[t_steps+1:2*t_steps].^2 
            SS2_Re_tmp[tr,:] +=(du[1:t_steps] .* sol[tr][1:t_steps,i] 
                               .- du[t_steps+1:2*t_steps] .* sol[tr][t_steps+1:2*t_steps,i])
            SS_Im_tmp[tr,:] += 2*du[1:t_steps] .* du[t_steps+1:2*t_steps]
            SS2_Im_tmp[tr,:] += (du[1:t_steps] .* sol[tr][t_steps+1:2*t_steps,i] 
                               .+ du[t_steps+1:2*t_steps] .* sol[tr][1:t_steps,i])
        end
    end
    

    SS_Re = transpose(mean(((args.dt/2)*(0.5-args.θ))^2 .* SS_Re_tmp ./ length(sol[1][1,:]),dims=1))
    SS2_Re = -transpose(mean((args.dt*(0.5-args.θ)) .* SS2_Re_tmp ./ length(sol[1][1,:]),dims=1))
    SS_Im = transpose(mean(((args.dt/2)*(0.5-args.θ))^2 .* SS_Im_tmp ./ length(sol[1][1,:]),dims=1))
    SS2_Im = -transpose(mean((args.dt*(0.5-args.θ)) .* SS2_Im_tmp ./ length(sol[1][1,:]),dims=1))
    
    return SS_Re, SS2_Re, SS_Im, SS2_Im
end