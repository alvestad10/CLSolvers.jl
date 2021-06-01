export plot_realtime_corr0t,plot_realtime_corr0t!, plot_full_corr0t,
       plot_Gpm_Gmp, plot_eucledian_corr0t, plot_xi_param_corr0t,
       plot_avg_avg2, plot_avg_avg2_nonEquil, 
       plot_corrections, plot_avg2_with_corrections

####################################################

using Plots,LaTeXStrings
using StatsPlots
using Measures


M_SIZE = 1.8
SOL_LW = 0.8

plot_setup = Dict(:legend     => :topright, #(0.1, 0.8), #(1.08, 0.6),
                  :grid       => false,
                  :foreground_color_legend=> nothing,
                  :background_color_legend=> nothing,
                  :framestyle => :box,
                  :thickness_scaling => 1.5)

markers_dict(color,markershape) = Dict(
                        :color => color,
                        :markershape => markershape,
                        :markersize  => M_SIZE,
                        :markerstrokecolor => color,
                        :lw => false)

markers_dict(color) = markers_dict(color,:square) 

solution_line_dict = Dict(:color => "black",:lw => SOL_LW)

color_1 = 1
color_2 = 2
color_3 = "black"
color_4 = 3
 
 
## PLOT Real-Time CORR0T
function plot_realtime_corr0t(corr0t_Re, corr0t_Im, CC;
        true_solution=0,true_solution2=0,label="",theta=0.1,
        save=false,fig_name="realtime_corrx0xt")

    x = real(CC)
    tpmax = argmax(x)

    fig = plot(xlabel=L"$x_0$",ylabel= true ? "" : L"$\langle \phi(0) \phi(t) \rangle$";plot_setup...)#,ylim=[-0.3,0.32])#,xlim=[0.3,0.501])#,ylim=[0.15,0.25],xlim=[0.3,0.501])

    if true_solution != 0
        tmax_inx = argmin(abs.(true_solution["x"] .- maximum(x[1:tpmax])))
        plot!(fig,true_solution["x"][1:tmax_inx],true_solution["corr0tRe"][1:tmax_inx],label=L"$\textrm{Solution}$", color=color_3, lw=SOL_LW)
        plot!(fig,true_solution["x"][1:tmax_inx],true_solution["corr0tIm"][1:tmax_inx],label=false, color=color_3, lw=SOL_LW)
    end

    if true_solution2 != 0
        tmax_inx = argmin(abs.(true_solution2["x"] .- maximum(x[1:tpmax])))
        plot!(fig,true_solution2["x"][1:tmax_inx],true_solution2["corr0tRe"][1:tmax_inx],label=L"$\textrm{Solution Tilted}$",color=color_4)
        plot!(fig,true_solution2["x"][1:tmax_inx],true_solution2["corr0tIm"][1:tmax_inx],label=false,color=color_4)
    end

    fig = plot_realtime_corr0t!(fig, corr0t_Re, corr0t_Im, CC; label=label,theta=theta)
    #lens!([0.39,0.505], [-0.205,-0.165], inset = (1,bbox(0.65, 0.4, 0.3, 0.3) ),xticks= nothing,yticks=nothing,framestyle = :box,subplot = 2)
    #lens!([0.39,0.505], [-0.205,-0.165], inset = (1,bbox(0.68, 0.4, 0.3, 0.3) ),xticks= nothing,yticks=nothing,framestyle = :box,subplot = 2)
    #lens!([0.41,0.505], [0.155,0.21], inset = (1, bbox(0.20, 0.13, 0.3, 0.3)), xticks= nothing,yticks=nothing, framestyle = :box, subplot=3)
    #lens!([0.41,0.505], [0.155,0.21], inset = (1, bbox(0.37, 0.22, 0.3, 0.3)), xticks= nothing,yticks=nothing, framestyle = :box, subplot=3)
    
    if save
        savefig(fig, string("Imgs/",fig_name,".pdf"))
    end
    return fig
end

function plot_realtime_corr0t!(fig,corr0t_Re, corr0t_Im, CC; label="",theta=0.1)
    x = real(CC)
    tpmax = argmax(x)

    scatter!(fig,x[1:tpmax],corr0t_Re[1:tpmax],label= (label=="" ? L"$\textrm{Re} \; G_{++}(x_0)$" : label);markers_dict(color_1,:square)...)#color=color_1, markershape=:square,markersize=M_SIZE,markerstrokecolor = color_1,lw=0.2)
    scatter!(fig,x[1:tpmax],corr0t_Im[1:tpmax],label= (label=="" ? L"$\textrm{Im} \; G_{++}(x_0)$" : label);markers_dict(color_2,:utriangle)...)#color=color_2, markershape=:square,markersize=M_SIZE,markerstrokecolor = color_2,lw=0.2)
    return fig
end


## PLOT Eucledian CORR0T
function plot_full_corr0t(corr0t_Re, corr0t_Im, CC;
    true_solution=0,label="",
    save=false,fig_name="realtime_corrx0xt0t_full")
    x = real(CC)
    tpmax = argmax(x)


    fig = plot(xlabel = L"$\xi$",
    ylabel = true ? "" : L"$ \langle \phi(0) \phi(\xi) \rangle - \langle \phi(\xi) \rangle\langle \phi(0) \rangle $";plot_setup...)#, ylim=[-0.25,0.32])
    tmax_inx = argmin(abs.(true_solution["x"] .- maximum(x[1:tpmax])))
    plot!(fig,true_solution["x"][1:tmax_inx],true_solution["corr0tRe"][1:tmax_inx],label=L"$\textrm{Solution}$",color=color_3, lw=SOL_LW)
    plot!(fig,true_solution["x"][1:tmax_inx],true_solution["corr0tIm"][1:tmax_inx],label=false,color=color_3, lw=SOL_LW)

    plot!(fig,true_solution["x"][tmax_inx] .+ true_solution["x"][1:tmax_inx],
    reverse(true_solution["corr0tRe"][1:tmax_inx]),label=false,color=color_3,lw=SOL_LW)
    plot!(fig,true_solution["x"][tmax_inx] .+ true_solution["x"][1:tmax_inx],
    reverse(true_solution["corr0tIm"][1:tmax_inx]),label=false,color=color_3,lw=SOL_LW)

    # Forward
    fw_interval = 1:tpmax
    scatter!(fig,x[fw_interval],corr0t_Re[fw_interval],
    label= (label=="" ? L"$\textrm{Re}\; G_{++}(\xi)$" : label); markers_dict(color_1)...)
    scatter!(fig,x[fw_interval],corr0t_Im[fw_interval],
    label= (label=="" ? L"$\textrm{Im}\; G_{++}(\xi)$" : label); markers_dict(color_2,:utriangle)...)


    # Backward
    bw_interval = tpmax+1:2*tpmax-1
    scatter!(fig,x[tpmax] .+ (x[tpmax] .- x[bw_interval]),corr0t_Re[bw_interval],
    label= L"$\textrm{Re}\; G_{+-}(\xi)$";markers_dict(3,:dtriangle)...)
    scatter!(fig,x[tpmax] .+ (x[tpmax] .- x[bw_interval]),corr0t_Im[bw_interval],
    label= L"$\textrm{Im}\; G_{+-}(\xi)$";markers_dict(4,:star)...)

    lens!([0.4,0.505], [-0.21,-0.17], inset = (1, bbox(0.55, 0.48, 0.25, 0.25)),xticks= nothing,yticks=nothing,framestyle = :box,subplot = 2)
    lens!([0.4,0.505], [0.15,0.21], inset = (1, bbox(0.7, 0.15, 0.25, 0.25)), xticks= nothing,yticks=nothing, framestyle = :box, subplot=3)


    if save
        savefig(fig, string("Imgs/",fig_name,".pdf"))
    end
    display(fig)
end

## PLOT Eucledian CORR0T
function plot_Gpm_Gmp(Gpm_Re, Gpm_Im, Gmp_Re, Gmp_Im, CC;
    true_solution=0,label="",
    save=false,fig_name="realtime_corrx0xt0t_full")
    
    x = real(CC)
    tpmax = argmax(x)


    fig = plot(xlabel = L"$x_0$", 
    ylabel = true ? "" : L"$ \langle \phi(0) \phi(t) \rangle - \langle \phi(x_0) \rangle\langle \phi(0) \rangle $";plot_setup...)#, ylim=[-0.25,0.32])
    tmax_inx = argmin(abs.(true_solution["x"] .- maximum(x[1:tpmax])))
    plot!(fig,true_solution["x"][1:tmax_inx],true_solution["corr0tRe"][1:tmax_inx],label=L"$\textrm{Solution } G_{+-}$",color=color_3, lw=SOL_LW)
    plot!(fig,true_solution["x"][1:tmax_inx],true_solution["corr0tIm"][1:tmax_inx],label=false,color="black", lw=SOL_LW)
    plot!(fig,true_solution["x"][1:tmax_inx],(-1) .* true_solution["corr0tIm"][1:tmax_inx],label=L"$\textrm{Solution } G_{-+}$",color="gray", lw=SOL_LW)

    # Forward
    fw_interval = 1:tpmax
    scatter!(fig,x[fw_interval],Gpm_Re[fw_interval],
    label= (label=="" ? L"$\textrm{Re} \;G_{+-}(x_0)$" : label);markers_dict(1,:square)...)
    scatter!(fig,x[fw_interval],Gpm_Im[fw_interval],
    label= (label=="" ? L"$\textrm{Im} \;G_{+-}(x_0)$" : label);markers_dict(3,:square)...)

    # Forward
    fw_interval = 1:tpmax
    scatter!(fig,x[fw_interval],Gmp_Re[fw_interval],
    label= (label=="" ? L"$\textrm{Re} \; G_{-+}(x_0)$" : label);markers_dict(2,:dtriangle)...)
    scatter!(fig,x[fw_interval],Gmp_Im[fw_interval],
    label= (label=="" ? L"$\textrm{Im} \; G_{-+}(x_0)$" : label);markers_dict(4,:dtriangle)...)

    if save
        savefig(fig, string("Imgs/",fig_name,".pdf"))
    end
    display(fig)
end


## PLOT Eucledian CORR0T
function plot_eucledian_corr0t(corr0t_Re, corr0t_Im, CC; true_solution=0,
            save=false,fig_name="eucledian_corrx0xt")

    x = -1 .* imag(CC)
    fig = plot(xlabel=L"$\xi$", ylabel=true ? "" : L"$\langle \phi(0)\phi(\xi) \rangle - \langle \phi(0) \rangle \langle\phi(\xi) \rangle$";
    ylim=[-Inf,Inf],plot_setup...)
    plot!(1 .+ true_solution["x"],true_solution["corr0tRe"],label=L"$\textrm{Solution}$", color=color_3,lw=SOL_LW)
    plot!(fig,1 .+ true_solution["x"],true_solution["corr0tIm"],label=false, color=color_3,lw=SOL_LW)
    scatter!(fig,1 .+ vcat(x[ceil(Int64,length(x)/2)+1:end],[1]),vcat(corr0t_Re[ceil(Int64,length(x)/2)+1:end],corr0t_Re[1]),label=L"$\textrm{Re}\; G_E(\xi)$",color=color_1, markershape=:square,markersize=M_SIZE,markerstrokecolor = color_1,lw=0)
    scatter!(fig,1 .+ vcat(x[ceil(Int64,length(x)/2)+1:end],[1]),vcat(corr0t_Im[ceil(Int64,length(x)/2)+1:end],corr0t_Im[1]),label=L"$\textrm{Im}\; G_E(\xi)$",color=color_2, markershape=:utriangle,markersize=M_SIZE,markerstrokecolor = color_2,lw=0)

    if save
    savefig(fig, string("Imgs/",fig_name,".pdf"))
    end
    display(fig)
end

## PLOT corr0t along the xi contour parmeter
function plot_xi_param_corr0t(corr0t_Re, corr0t_Im, tp; true_solution=0)
    fig = plot(legend=:bottomright)#,ylim=[-0.3,0.4])
    #plot!(true_solution["x"],true_solution["corr0tRe"],label="ScroSol Re")
    #plot!(fig,true_solution["x"],true_solution["corr0tIm"],label="ScroSol Im")
    plot!(fig,tp,corr0t_Re,label=L"$Re(\phi(0)\phi(t))$",color=1, markershape=:square,markerstrokecolor = 1,lw=0.1)
    plot!(fig,tp,corr0t_Im,label=L"$Im(\phi(0)\phi(t))$",color=2, markershape=:square,markerstrokecolor = 2,lw=0.1)
end


## PLOT AVG and AVG2
function plot_avg_avg2(avg_Re, avg_Im, avg2_Re, avg2_Im, tp; true_solution=0,save=false,fig_name="param_corrx")

    fig = plot(xlabel=L"$\xi$"; plot_setup...)#,ylim=[-0.02,0.325])

    tpmax = argmax(tp)
    tmax_inx = argmin(abs.(true_solution["x"] .- maximum(tp[1:tpmax])))
    plot!(fig,2*true_solution["x"][1:tmax_inx],true_solution["avgRe"][1:tmax_inx],label=L"$\textrm{Solution}$"; solution_line_dict...)
    plot!(fig,2*true_solution["x"][1:tmax_inx],true_solution["avg2Re"][1:tmax_inx],label=false; solution_line_dict...)

    scatter!(fig,2*tp,avg_Re,label=L"$\textrm{Re}\langle \phi \rangle $"; markers_dict(1, :square)...)
    scatter!(fig,2*tp,avg_Im,label=L"$\textrm{Im}\langle \phi \rangle $"; markers_dict(2, :utriangle)...)
    scatter!(fig,2*tp,avg2_Re,label=L"$\textrm{Re}\langle \phi^2\rangle $"; markers_dict(3, :star)...)
    scatter!(fig,2*tp,avg2_Im,label=L"$\textrm{Im}\langle \phi^2 \rangle $"; markers_dict(4, :xcross)...)

    if save
        savefig(fig, string("Imgs/",fig_name,".pdf"))
    end
    display(fig)
end


## PLOT AVG and AVG
function plot_avg_avg2_nonEquil(avg_Re, avg_Im, avg2_Re, avg2_Im, CC; true_solution=0,save=false,fig_name="param_corrx")

    x = real(CC)
    tpmax = argmax(x)

    fig = plot(xlabel=L"$\xi$"; plot_setup...)

    x_true = true_solution["x"]
    if true_solution != 0
        tmax_inx = argmin(abs.(true_solution["x"] .- maximum(x[1:tpmax])))

        plot!(fig,x_true[1:tmax_inx],true_solution["avgRe"][1:tmax_inx],
        label=L"\textrm{Solution}"; solution_line_dict...)
        plot!(fig,x_true[1:tmax_inx],true_solution["avgIm"][1:tmax_inx],
        label=false;solution_line_dict...)

        plot!(fig,x_true[tmax_inx] .+ x_true[1:tmax_inx],
        reverse(true_solution["avgRe"][1:tmax_inx]),
        label=false;solution_line_dict...)
        plot!(fig,x_true[tmax_inx] .+ x_true[1:tmax_inx],
        reverse(true_solution["avgIm"][1:tmax_inx]),
        label=false;solution_line_dict...)

        plot!(fig,x_true[1:tmax_inx],true_solution["avg2Re"][1:tmax_inx],
        label=false; solution_line_dict...)
        plot!(fig,x_true[1:tmax_inx],true_solution["avg2Im"][1:tmax_inx],
        label=false;solution_line_dict...)

        plot!(fig,x_true[tmax_inx] .+ x_true[1:tmax_inx],
        reverse(true_solution["avg2Re"][1:tmax_inx]),
        label=false;solution_line_dict...)
        plot!(fig,x_true[tmax_inx] .+ x_true[1:tmax_inx],
        reverse(true_solution["avg2Im"][1:tmax_inx]),
        label=false;solution_line_dict...)
    end


    # Forward
    fw_interval = 1:tpmax
    scatter!(fig,x[fw_interval],avg_Re[fw_interval],
    label= L"$\textrm{Re}\langle \phi \rangle$"; markers_dict(1,:square)...)
    scatter!(fig,x[fw_interval],avg_Im[fw_interval],
    label= L"$\textrm{Im}\langle \phi \rangle$"; markers_dict(2,:utriangle)...)

    # Forward
    fw_interval = 1:tpmax
    scatter!(fig,x[fw_interval],avg2_Re[fw_interval],
    label= L"$\textrm{Re}\langle \phi^2 \rangle$"; markers_dict(3,:star)...)
    scatter!(fig,x[fw_interval],avg2_Im[fw_interval],
    label= L"$\textrm{Im}\langle \phi^2 \rangle$"; markers_dict(4,:cross)...)


    # Backward
    bw_interval = tpmax+1:2*tpmax-1
    plot!(fig,x[tpmax] .+ (x[tpmax] .- x[bw_interval]),avg_Re[bw_interval],
    label= false; markers_dict(1,:square)...)
    plot!(fig,x[tpmax] .+ (x[tpmax] .- x[bw_interval]),avg_Im[bw_interval],
    label= false; markers_dict(2,:utriangle)...)

    # Backward
    bw_interval = tpmax+1:2*tpmax-1
    plot!(fig,x[tpmax] .+ (x[tpmax] .- x[bw_interval]),avg2_Re[bw_interval],
    label= false; markers_dict(3,:star)...)
    plot!(fig,x[tpmax] .+ (x[tpmax] .- x[bw_interval]),avg2_Im[bw_interval],
    label= false; markers_dict(4,:cross)...)


    if save
        savefig(fig, string("Imgs/",fig_name,".pdf"))
    end
    display(fig)
end




function plot_corrections(SS_Re, SS2_Re, SS_Im, SS2_Im, avg2_Re, avg2_Im, tp;true_solution=0, save=false)
    
    fig_Re = plot(xlabel="";plot_setup...,xticks=nothing,bottom_margin=-3mm)#,ylim=[-0.15,0.015])
    
    if true_solution != 0
        scatter!(fig_Re,2*tp[1:end-1],(avg2_Re .- true_solution[1] ); markers_dict(2,:square)...,
                 label=L"$\textrm{Re}\langle \Delta \phi^2 \rangle $")
    end

    scatter!(fig_Re,2*tp[1:end-1], SS2_Re; markers_dict(1,:circle)...,
             label=L"$\textrm{Re} \; \Sigma$")
    scatter!(fig_Re,2*tp[1:end-1],SS_Re .+ SS2_Re; markers_dict(5,:utriangle)...,
             label=L"$\textrm{Re}\left\langle \left(\tilde\phi - \phi\right)^2 \right\rangle $")

    fig_Im = plot(xlabel=L"$\xi$"; plot_setup...,top_margin=0mm)#,ylim=[-0.015,0.009])

    if true_solution != 0
        scatter!(fig_Im,2*tp[1:end-1],avg2_Im; markers_dict(4,:square)...,
                 label=L"$\textrm{Im}\langle \Delta \phi^2 \rangle $")
    end
    scatter!(fig_Im,2*tp[1:end-1],SS2_Im; markers_dict(3,:circle)...,
             label=L"$\textrm{Im} \; \Sigma$")
    scatter!(fig_Im,2*tp[1:end-1], SS_Im .+ SS2_Im; markers_dict(6,:utriangle)...,
             label=L"$\textrm{Im}\left\langle \left(\tilde\phi - \phi\right)^2 \right\rangle $")
    
    fig = plot(fig_Re, fig_Im, layout = (2, 1),size=(750, 400),right_margin=25mm)
    if save
        savefig(fig, string("Imgs/Correcting_vs_difference_HO_Avg2.pdf"))
    end
    display(fig)
end

function plot_avg2_with_corrections(true_solution=0)
    fig = plot(xlabel=L"$\xi$";plot_setup...)
    
    plot!(fig,2*tp[1:end-1],SS_Re .+ SS2_Re;
                    markers_dict(1,:square)...,label=L"$\textrm{Re}\langle \left(\tilde\phi - \phi\right)^2 \rangle $")
    plot!(fig,2*tp[1:end-1], SS2_Re;
                    markers_dict(1,:square)...,alpha=0.3,label=L"$\textrm{Re}\left\langle \phi (iS) \right\rangle$")
    plot!(fig,2*tp[1:end-1],(avg2_Re .- dfScroSolMinkHO["avg2Re"][1] );
                    markers_dict(2,:utriangle)...,label=L"$\textrm{Re}\langle \Delta \phi^2 \rangle $")
    
    plot!(fig,2*tp[1:end-1], SS_Im .+ SS2_Im;
                    markers_dict(3,:star)...,label=L"$\textrm{Im}\langle \left(\tilde\phi - \phi\right)^2 \rangle $")
    plot!(fig,2*tp[1:end-1],SS2_Im;
                    markers_dict(3,:star)...,alpha=0.3,label=L"$\textrm{Im}\left\langle \phi (iS) \right\rangle$")
    plot!(fig,2*tp[1:end-1],avg2_Im;
                    markers_dict(4,:utriangle)...,label=L"$\textrm{Im}\langle \Delta \phi^2 \rangle $")

    if save
        savefig(fig, string("Imgs/Correcting_vs_difference_AHO_Avg2.pdf"))
    end
    
    display(fig)
end