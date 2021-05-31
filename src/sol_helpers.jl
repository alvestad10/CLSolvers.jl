export check_errors, save, load_sol

######################################

using JLD2
using FileIO
using UUIDs
using Dates
using DataFrames

######################
# Helper methods for the solution of the solve method in stochasticDiffEe.jl library
######################

"""
   CHECK FOR ERRORS/WARNINGS DURING RUN

   can also remove warning trajectories by setting
                remove_warning_tr=true
"""
function check_errors(sol;remove_warning_tr=false)
    allok = true
    warnings_inxes = []
    for (i,s) in enumerate(sol) 
        if (s.retcode != :Success)
            allok = false
            println("Trajectory ",i," has retcode: ", string(s.retcode))

            push!(warnings_inxes,i) 

        end
    end
        
    if remove_warning_tr && !allok
        deleteat!(sol.u,warnings_inxes)
        println("Removing trajectories ", warnings_inxes)
    else
        println("All ok")
    end 


end



function save(path, sol, args, model::AbstractModel, contour::AbstractContour, a, CC, tp, discretizationContour::AbstractContourDiscretization)
    ID = uuid1()

    if !isdir(path)
        mkdir(path)
    end

    save_sol(ID,path,sol,model,contour,a,CC,tp)
    make_config_file(ID,path,sol,args,contour,discretizationContour,model)
end


"""
   Load the sol object from filename

   returns
        sol: sol object from stochastifDiffEq.jl solver
        p: parameters used
        contour: contour
        a: lattice spacings (complex)
        CC: lattice points (complex)
        tp: t_param points (real)
"""
function save_sol(ID,path, sol, p, contour, a, CC, tp)
    

    path_sol = joinpath(path,"sol")
    if !isdir(path_sol)
        mkdir(path_sol)
    end

    filename = joinpath(path_sol,string(ID))
    
    @save filename :sol=sol :p=p :contour=contour :a=a  :CC=CC :tp=tp
end

"""
   Save the sol object with filename

   returns
        sol: sol object from stochastifDiffEq.jl solver
        p: parameters used
        contour: contour
        a: lattice spacings (complex)
        CC: lattice points (complex)
        tp: t_param points (real)
"""
function save_sol(filename, sol, p, contour, a, CC, tp)
    
    @save filename :sol=sol :p=p :contour=contour :a=a  :CC=CC :tp=tp
end


"""
    Make simulaiton configuration file
"""
function make_config_file(ID,path, sol, args, contour::AbstractContour, discretizationContour::AbstractContourDiscretization, model::AbstractModel)

    path_config = joinpath(path,"config")

    if !isdir(path_config)
        mkdir(path_config)
    end

    filename = joinpath(path_config,string(ID))

    line_shift = "\n"

    open(filename, "w") do io
        write(io, string("ID: ",        ID,                line_shift,
                  "Time: ",      Time(Dates.now()), line_shift,
                  "Date: ",      Date(Dates.now()), line_shift,
                  "Time_used: ", sol.elapsedTime,         line_shift)
            )
        write(io,"************ CONTOUR *************",line_shift)
        write(io, "contour: ", string(typeof(contour)), line_shift)
        write(io, string(contour),              line_shift)

        write(io,"************ DISCRETIZATION *************",line_shift)
        write(io, "discretizationContour: ", string(typeof(discretizationContour)), line_shift)
        write(io, string(discretizationContour), line_shift)

        write(io, "*********** MODEL ***********",line_shift)
        write(io, "model: ", string(typeof(model)), line_shift)
        write(io, string(model), line_shift)
            
        write(io, "*********** SOLVER ***********",line_shift)
        write(io, "algorithm: ", split(string(typeof(sol[1].alg)),"{")[1], line_shift)
        if hasproperty(sol[1].alg,:theta)
            write(io, "theta: ", string(sol[1].alg.theta), line_shift)
            write(io, "tol: ", string(args.tol), line_shift)
            write(io, "dtmax: ", string(args.dtmax), line_shift)
        end
        write(io, "tspan: ", string(sol[1].prob.tspan), line_shift)
        write(io, "dt: ", string(args.dt), line_shift)
        write(io, "N_Tr: ", string(args.N_Tr), line_shift)


            
    end
end


"""
   Load the sol object from filename

   returns
        sol: sol object from stochastifDiffEq.jl solver
        p: parameters used
        contour: contour
        a: lattice spacings (complex)
        CC: lattice points (complex)
        tp: t_param points (real)
"""
function load_sol(filename)
    lo = load(filename)
    sol = lo[":sol"]
    p = lo[":p"]
    contour = lo[":contour"]
    a = lo[":a"]
    CC = lo[":CC"]
    tp = lo[":tp"]

    return sol, p, contour, a, CC, tp
end



"""
   Read config file

   returns
        dictionary of all the config options
"""
function read_config(file_path)
    #df = DataFrame(ID = String[], Time = String[], Date = String[], Time_used = Float64[],
    #               contour = String[], tmax=Float64[],β=Float64[],A_β=Float64[],F_β=Float64[],
    #               discretizationContour=String[], NR_Points=Int[],NR_RT_FW_Points=Int[],
    #               NR_RT_BW_Points=Int[],NR_ET_Points=Int[],
    #               model=String[],g=Float64[],m=Float64[],λ=Float64[],
    #               algorithm=)
    D = Dict{String,Any}("filename"=>file_path)
    open(file_path) do f
        for l in eachline(f)


            # Check for int
            if !isnothing(try D[split(l,":")[1]] = parse(Int,split(l,":")[2]) catch ArgumentError end)
                continue
            end

            # Check for int
            if !isnothing(try D[split(l,":")[1]] = parse(Float64,split(l,":")[2]) catch ArgumentError end)
                continue
            end

            try
                D[split(l,":")[1]] = strip(join(split(l,":")[2:end],":"))
            catch BoundsError
            end

        end
    end

    return D
end


"""
   Get all configurations in path

   returns
        dataframe of all config options as columns and 
"""
function get_configs(path)

    path_configs = joinpath(path,"configs")
    
    df = DataFrame()
    for file in readdir(path_configs,join=true)
        D = read_config(file)

        push!(df,D,cols=:union)
    end

    return df
end


"""
   Load from ID

   returns
        dataframe of all config options as columns and 
"""
function load_sol_from_ID(ID,path)
    
    filename = joinpath(joinpath(path,"sol"),ID)

    return load_sol(filename)
end


