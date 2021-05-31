export EucledianContourDiscretization,
       discretizeContour,
       getContour,
       getContourDistances,
       getContourDerivatives,
       spread_timepoints

#############################################################################

abstract type AbstractContourDiscretization end
abstract type AbstractSK_ContourDiscretization <: AbstractContourDiscretization end


"""
    Get distance between points in contour
"""
function getContourDistances(CD::AbstractContourDiscretization,
                             C::AbstractContour)
    contour = getContour(CD,C)
    N = length(contour)
    return [contour[i+1] - contour[i]  for i in 1:N-1]
end

"""
    Get the discretized contour
"""
function getContour(CD::AbstractContourDiscretization,
                    C::AbstractContour)
    ts = spread_timepoints(CD,C) 
    return map(t -> get_timepoint(t,C), ts)
end

"""
    Get the discretized contour derivatives
"""
function getContourDerivatives(CD::AbstractContourDiscretization,
                    C::AbstractContour)
    ts = spread_timepoints(CD,C) 
    return map(t -> get_derivative_timepoint(t,C), ts)
end


"""
    Discretizing the contour
"""
function discretizeContour(CD::EucledianContour, N_points)
    return EucledianContourDiscretization(N_points)
end

function discretizeContour(CD::NonEquilSchwingerKeldyshContour, points)
    return NonEquilSchwingerKeldyshContourDiscretization(points)
end

function discretizeContour(CD::NonEquilTiltedSchwingerKeldyshContour, points)
    return NonEquilTiltedSchwingerKeldyshContourDiscretization(points)
end

function discretizeContour(CD::SchwingerKeldyshContour, N_points)
    return SchwingerKeldyshContourDiscretization(N_points)
end


function discretizeContour(CD::SchwingerKeldyshContour, NR_RT_points, NR_ET_Points)
    return SchwingerKeldyshContourDiscretization(NR_RT_points,NR_ET_Points)
end

function discretizeContour(CD::SchwingerKeldyshContour, NR_RT_points_1, NR_RT_points_2, NR_ET_Points)
    return SchwingerKeldyshContourDiscretization(NR_RT_points_1, NR_RT_points_2,NR_ET_Points)
end

function discretizeContour(CD::SchwingerKeldyshContourRounded1, N_points_1, N_points_rounded)
    return SchwingerKeldyshRounded1ContourDiscretization(N_points_1,N_points_rounded)
end

function discretizeContour(CD::SchwingerKeldyshContourRounded1, N_points)
    N_points_1 = 3*N_points/4
    N_points_rounded = 1*N_points/4
    return SchwingerKeldyshRounded1ContourDiscretization(N_points_1,N_points_rounded)
end


function discretizeContour(CD::SchwingerKeldyshContourRounded2, N_points_1, N_points_rounded1, N_points_rounded2)
    return SchwingerKeldyshRounded2ContourDiscretization(N_points_1,N_points_rounded1,N_points_rounded2)
end

function discretizeContour(CD::SchwingerKeldyshContourRounded2, N_points)
    N_points_1 = N_points/2
    N_points_rounded1 = N_points/4
    N_points_rounded2 = N_points/4
    return SchwingerKeldyshRounded2ContourDiscretization(N_points_1,N_points_rounded1,N_points_rounded2)
end

function discretizeContour(CD::SchwingerKeldyshContourRounded3, N_points_1, N_points_rounded1, N_points_rounded2, N_points_rounded3)
    return SchwingerKeldyshRounded3ContourDiscretization(N_points_1,N_points_rounded1,N_points_rounded2,N_points_rounded3)
end

function discretizeContour(CD::SchwingerKeldyshContourRounded3, N_points)
    N_points_1 = N_points/2
    N_points_rounded1 = 2*N_points/8
    N_points_rounded2 = 1*N_points/8
    N_points_rounded3 = 1*N_points/8
    return SchwingerKeldyshRounded3ContourDiscretization(N_points_1,N_points_rounded1,N_points_rounded2,N_points_rounded3)
end

################# Eucledian contour #################
struct EucledianContourDiscretization <: AbstractContourDiscretization
    NR_Points::Integer
end

Base.string(CD::EucledianContourDiscretization) = string("NR_Points: ", CD.NR_Points)

function spread_timepoints(CD::EucledianContourDiscretization,
                           C::EucledianContour)
    return collect(range(0,stop=1,length=CD.NR_Points+1))
end


################# Non-Equilibrium Schwinger-Keldysh Contour ####################

"""
    Discretization of the Non-Equilibrium Schwinger Keldysh Contour
"""
struct NonEquilSchwingerKeldyshContourDiscretization <: AbstractSK_ContourDiscretization
    points::Integer
    function NonEquilSchwingerKeldyshContourDiscretization(points)
        @assert (points-1)%2 == 0
        new(points)
    end
end

Base.string(CD::NonEquilSchwingerKeldyshContourDiscretization) = string("NR_Points: ", CD.points)


"""
    Create the parameter array, i.e. the points corresponding to
    t in Reals in the parameterize curve

    This will include the periodic boundary value, i.e. the t=1.0
"""
function spread_timepoints(CD::NonEquilSchwingerKeldyshContourDiscretization,
                          C::NonEquilSchwingerKeldyshContour)
    timePoints = []
    N = floor(Int64,CD.points/2)
    for i in 0:CD.points-1
        if i < CD.points/2
            push!(timePoints, i*t1(C)/N)
        else
            i_eff = i - N
            push!(timePoints, t1(C) + i_eff*t2(C)/N)
        end

    end  
    return timePoints
end

################# Non-Equilibrium Schwinger-Keldysh Contour ####################

"""
    Discretization of the Non-Equilibrium Schwinger Keldysh Contour
"""
struct NonEquilTiltedSchwingerKeldyshContourDiscretization <: AbstractSK_ContourDiscretization
    points::Integer
    function NonEquilTiltedSchwingerKeldyshContourDiscretization(points)
        @assert (points-1)%2 == 0
        new(points)
    end
end

Base.string(CD::NonEquilTiltedSchwingerKeldyshContourDiscretization) = string("NR_Points: ", CD.points)


"""
    Create the parameter array, i.e. the points corresponding to
    t in Reals in the parameterize curve

    This will include the periodic boundary value, i.e. the t=1.0
"""
function spread_timepoints(CD::NonEquilTiltedSchwingerKeldyshContourDiscretization,
                          C::NonEquilTiltedSchwingerKeldyshContour)
    timePoints = []
    N = floor(Int64,CD.points/2)
    for i in 0:CD.points-1
        if i < CD.points/2
            push!(timePoints, i*t1(C)/N)
        else
            i_eff = i - N
            push!(timePoints, t1(C) + i_eff*t2(C)/N)
        end

    end  
    return timePoints
end


################## Schwinger-Keldysh Contour ####################

"""
    Discretization of the Schwinger Keldysh Contour
"""
struct SchwingerKeldyshContourDiscretization <: AbstractSK_ContourDiscretization
    NR_RT_FW_Points::Integer
    NR_RT_BW_Points::Integer
    NR_RT_Points::Integer
    NR_ET_Points::Integer
    function SchwingerKeldyshContourDiscretization(RTP1,RTP2,ETP)
        RTP = RTP1 + RTP2
        new(RTP1,RTP2,RTP,ETP)
    end
end

function Base.string(CD::SchwingerKeldyshContourDiscretization) 
        string("NR_Points: ", get_total_points(CD), "\n",
              "NR_RT_FW_Points: ", CD.NR_RT_FW_Points, "\n",
              "NR_RT_BW_Points: ", CD.NR_RT_BW_Points, "\n",
              "NR_ET_Points: ", CD.NR_ET_Points)
end

"""
Short initialization of the SchwingerKeldyshContourDiscretization struct

Separates the nr of points into 2/3 RT-points and 1/3 ET-points.
If not possible the RT-points will be rounded up and the nr of
ET points rounded down
"""
function SchwingerKeldyshContourDiscretization(NR_Points)
    @assert NR_Points>2
    
    # Special case for 4 points
    if NR_Points == 4 
        return SchwingerKeldyshContourDiscretization(2,2)
    end

    RTP = 2*ceil(Int64,NR_Points/3)
    RTP1 = ceil(Int64,RTP/2)
    ETP = floor(Int64, NR_Points/3 - (NR_Points%2==0 ? (NR_Points/3)%2 : (NR_Points%3)%2))
    

    SchwingerKeldyshContourDiscretization(RTP1,RTP1,ETP)
end

"""
Short initialization of the SchwingerKeldyshContourDiscretization struct

Separates the nr of points into 2/3 RT-points and 1/3 ET-points.
If not possible the RT-points will be rounded up and the nr of
ET points rounded down
"""
function SchwingerKeldyshContourDiscretization(RTP,ETP)
    @assert RTP>0

    RTP1 = floor(Int64,RTP/2)
    RTP2 = ceil(Int64,RTP/2)

    SchwingerKeldyshContourDiscretization(RTP1,RTP2,ETP)
end

"""
    Create the parameter array, i.e. the points corresponding to
    t in Reals in the parameterize curve

    This will include the periodic boundary value, i.e. the t=1.0
"""
function spread_timepoints(CD::SchwingerKeldyshContourDiscretization,
                          C::SchwingerKeldyshContour)
    timePoints = []
    for i in 0:get_total_points(CD)
        if i <= CD.NR_RT_FW_Points
            N = CD.NR_RT_FW_Points
            push!(timePoints, i*t1(C)/N)
        elseif i < CD.NR_RT_Points
            i_eff = i - CD.NR_RT_FW_Points
            N = CD.NR_RT_BW_Points  
            push!(timePoints, t1(C) + i_eff*t2(C)/N)
        else
            i_eff = i - CD.NR_RT_Points
            N = CD.NR_ET_Points

            if N == 0
                push!(timePoints, t1(C) + t2(C) + t3(C))
            else
                push!(timePoints, t1(C) + t2(C) + i_eff*t3(C)/N)
            end
        end

    end  
    return timePoints
end



################## Schwinger-Keldysh Contour with one rounded corner ####################

"""
    Schwinger-Keldysh contour with one rounded corners
"""
struct SchwingerKeldyshRounded1ContourDiscretization <: AbstractSK_ContourDiscretization
    NR_RT_Points::Integer
    NR_ET_Points::Integer
    NR_Rounded1_Points::Integer
    function SchwingerKeldyshRounded1ContourDiscretization(RTP,ETP,Rounded1)
        @assert RTP%2 == 0
        @assert RTP > 0
        @assert ETP > 0
        @assert Rounded1 > 2
        new(RTP,ETP,Rounded1)
    end
end

function Base.string(CD::SchwingerKeldyshRounded1ContourDiscretization) 
    string("NR_Points: ", get_total_points(CD), "\n",
          "NR_RT_Points: ", CD.NR_RT_FW_Points, "\n",
          "NR_ET_Points: ", CD.NR_ET_Points, "\n",
          "NR_Rounded1_Points: ", CD.NR_Rounded1_Points)
end

"""
Short initialization of the SchwingerKeldyshRounded1ContourDiscretization struct

Separates the nr of points into 2/3 RT-points and 1/3 ET-points.
If not possible the RT-points will be rounded up and the nr of
ET points rounded down
"""
function SchwingerKeldyshRounded1ContourDiscretization(NR_Points, NR_Rounded1)
    LinearPoints = SchwingerKeldyshContourDiscretization(NR_Points)
    RTP = LinearPoints.NR_RT_Points
    ETP = LinearPoints.NR_ET_Points
    SchwingerKeldyshRounded1ContourDiscretization(RTP,ETP,NR_Rounded1)
end

"""
    Create the parameter array, i.e. the points corresponding to
    t in Reals in the parameterize curve

    This will include the periodic boundary value, i.e. the t=1.0
"""
function spread_timepoints(CD::SchwingerKeldyshRounded1ContourDiscretization,
                          C::SchwingerKeldyshContourRounded1)
    timePoints = []
    for i in 0:get_total_points(CD)
        if i <= CD.NR_RT_Points / 2
            N = CD.NR_RT_Points / 2
            push!(timePoints, i*t1(C)/N)

        elseif i <= CD.NR_RT_Points / 2 + CD.NR_Rounded1_Points / 2
            i_eff = i - (CD.NR_RT_Points / 2)
            N = CD.NR_Rounded1_Points / 2
            push!(timePoints, t1(C) + i_eff*t2(C)/N)

        elseif i <= CD.NR_RT_Points / 2 + CD.NR_Rounded1_Points
            i_eff = i - (CD.NR_RT_Points / 2) - (CD.NR_Rounded1_Points / 2)
            N = CD.NR_Rounded1_Points / 2
            push!(timePoints, t1(C) + t2(C) + i_eff*t3(C)/N)

        elseif i <= CD.NR_RT_Points + CD.NR_Rounded1_Points
            i_eff = i - (CD.NR_RT_Points / 2) - CD.NR_Rounded1_Points
            N = CD.NR_RT_Points / 2
            push!(timePoints, t1(C) + t2(C) + t3(C) + i_eff*t4(C)/N)

        else
            i_eff = i - CD.NR_RT_Points - CD.NR_Rounded1_Points
            N = CD.NR_ET_Points + 1
            push!(timePoints, t1(C) + t2(C) + t3(C) + t4(C) + i_eff*t5(C)/N)
        end

    end  
    return timePoints
end

################## Schwinger-Keldysh Contour with two rounded corner ####################

"""
    Schwinger-Keldysh contour with one rounded corners
"""
struct SchwingerKeldyshRounded2ContourDiscretization <: AbstractSK_ContourDiscretization
    NR_RT_Points
    NR_ET_Points
    NR_Rounded1_Points
    NR_Rounded2_Points
    function SchwingerKeldyshRounded2ContourDiscretization(RTP,ETP,Rounded1,Rounded2)
        @assert RTP%2 == 0
        @assert RTP > 0
        @assert ETP > 0
        @assert Rounded1 > 2
        new(RTP,ETP,Rounded1,Rounded2)
    end
end

function Base.show(io::IO, CD::SchwingerKeldyshRounded2ContourDiscretization) 
    string("NR_Points: ", get_total_points(CD), "\n",
          "NR_RT_Points: ", CD.NR_RT_FW_Points, "\n",
          "NR_ET_Points: ", CD.NR_ET_Points, "\n",
          "NR_Rounded1_Points: ", CD.NR_Rounded1_Points, "\n",
          "NR_Rounded2_Points: ", CD.NR_Rounded2_Points)
end

"""
Short initialization of the SchwingerKeldyshRounded1ContourDiscretization struct

Separates the nr of points into 2/3 RT-points and 1/3 ET-points.
If not possible the RT-points will be rounded up and the nr of
ET points rounded down
"""
function SchwingerKeldyshRounded2ContourDiscretization(NR_Points, NR_Rounded1, NR_Rounded2)
    LinearPoints = SchwingerKeldyshContourDiscretization(NR_Points)
    RTP = LinearPoints.NR_RT_Points
    ETP = LinearPoints.NR_ET_Points
    SchwingerKeldyshRounded2ContourDiscretization(RTP,ETP,NR_Rounded1, NR_Rounded2)
end

"""
    Create the parameter array, i.e. the points corresponding to
    t in Reals in the parameterize curve

    This will include the periodic boundary value, i.e. the t=1.0
"""
function spread_timepoints(CD::SchwingerKeldyshRounded2ContourDiscretization,
                          C::SchwingerKeldyshContourRounded2)
    timePoints = []
    for i in 0:get_total_points(CD)
        if i <= CD.NR_RT_Points / 2
            N = CD.NR_RT_Points / 2
            push!(timePoints, i*t1(C)/N)

        elseif i <= CD.NR_RT_Points / 2 + CD.NR_Rounded1_Points / 2
            i_eff = i - (CD.NR_RT_Points / 2)
            N = CD.NR_Rounded1_Points / 2
            push!(timePoints, t1(C) + i_eff*t2(C)/N)

        elseif i <= CD.NR_RT_Points / 2 + CD.NR_Rounded1_Points
            i_eff = i - (CD.NR_RT_Points / 2) - (CD.NR_Rounded1_Points / 2)
            N = CD.NR_Rounded1_Points / 2
            push!(timePoints, t1(C) + t2(C) + i_eff*t3(C)/N)

        elseif i <= CD.NR_RT_Points + CD.NR_Rounded1_Points
            i_eff = i - (CD.NR_RT_Points / 2) - CD.NR_Rounded1_Points
            N = CD.NR_RT_Points / 2
            push!(timePoints, t1(C) + t2(C) + t3(C) + i_eff*t4(C)/N)

        elseif i <= CD.NR_RT_Points + CD.NR_Rounded1_Points + CD.NR_Rounded2_Points
            i_eff = i - CD.NR_RT_Points - CD.NR_Rounded1_Points
            N = CD.NR_Rounded2_Points
            push!(timePoints, t1(C) + t2(C) + t3(C) + t4(C) + i_eff*t5(C)/N)

        else
            i_eff = i - CD.NR_RT_Points - CD.NR_Rounded1_Points - CD.NR_Rounded2_Points
            N = CD.NR_ET_Points + 1
            push!(timePoints, t1(C) + t2(C) + t3(C) + t4(C) + t5(C) + t5(C) + i_eff*t6(C)/N)
        end

    end  
    return timePoints
end


################## Schwinger-Keldysh Contour with three rounded corner ####################

"""
    Schwinger-Keldysh contour with one rounded corners
"""
struct SchwingerKeldyshRounded3ContourDiscretization <: AbstractSK_ContourDiscretization
    NR_RT_Points
    NR_ET_Points
    NR_Rounded1_Points
    NR_Rounded2_Points
    NR_Rounded3_Points
    function SchwingerKeldyshRounded3ContourDiscretization(RTP,ETP,Rounded1,Rounded2,Rounded3)
        @assert RTP%2 == 0
        @assert RTP > 0
        @assert ETP > 0
        @assert Rounded1 > 2
        new(RTP,ETP,Rounded1,Rounded2,Rounded3)
    end
end

function Base.string(CD::SchwingerKeldyshRounded3ContourDiscretization) 
    string("NR_Points: ", get_total_points(CD), "\n",
          "NR_RT_Points: ", CD.NR_RT_FW_Points, "\n",
          "NR_ET_Points: ", CD.NR_ET_Points, "\n",
          "NR_Rounded1_Points: ", CD.NR_Rounded1_Points, "\n",
          "NR_Rounded2_Points: ", CD.NR_Rounded2_Points, "\n",
          "NR_Rounded3_Points: ", CD.NR_Rounded3_Points)
end

"""
Short initialization of the SchwingerKeldyshRounded1ContourDiscretization struct

Separates the nr of points into 2/3 RT-points and 1/3 ET-points.
If not possible the RT-points will be rounded up and the nr of
ET points rounded down
"""
function SchwingerKeldyshRounded3ContourDiscretization(NR_Points, NR_Rounded1, NR_Rounded2, NR_Rounded3)
    LinearPoints = SchwingerKeldyshContourDiscretization(NR_Points)
    RTP = LinearPoints.NR_RT_Points
    ETP = LinearPoints.NR_ET_Points
    SchwingerKeldyshRounded3ContourDiscretization(RTP,ETP,NR_Rounded1, NR_Rounded2, NR_Rounded3)
end

"""
    Create the parameter array, i.e. the points corresponding to
    t in Reals in the parameterize curve

    This will include the periodic boundary value, i.e. the t=1.0
"""
function spread_timepoints(CD::SchwingerKeldyshRounded3ContourDiscretization,
                          C::SchwingerKeldyshContourRounded3)
    timePoints = []
    for i in 0:get_total_points(CD)
        if i <= CD.NR_RT_Points / 2
            N = CD.NR_RT_Points / 2
            push!(timePoints, i*t1(C)/N)

        elseif i <= CD.NR_RT_Points / 2 + CD.NR_Rounded1_Points / 2
            i_eff = i - (CD.NR_RT_Points / 2)
            N = CD.NR_Rounded1_Points / 2
            push!(timePoints, t1(C) + i_eff*t2(C)/N)

        elseif i <= CD.NR_RT_Points / 2 + CD.NR_Rounded1_Points
            i_eff = i - (CD.NR_RT_Points / 2) - (CD.NR_Rounded1_Points / 2)
            N = CD.NR_Rounded1_Points / 2
            push!(timePoints, t1(C) + t2(C) + i_eff*t3(C)/N)

        elseif i <= CD.NR_RT_Points + CD.NR_Rounded1_Points
            i_eff = i - (CD.NR_RT_Points / 2) - CD.NR_Rounded1_Points
            N = CD.NR_RT_Points / 2
            push!(timePoints, t1(C) + t2(C) + t3(C) + i_eff*t4(C)/N)

        elseif i <= CD.NR_RT_Points + CD.NR_Rounded1_Points + CD.NR_Rounded2_Points
            i_eff = i - CD.NR_RT_Points - CD.NR_Rounded1_Points
            N = CD.NR_Rounded2_Points
            push!(timePoints, t1(C) + t2(C) + t3(C) + t4(C) + i_eff*t5(C)/N)

        elseif i < CD.NR_RT_Points + CD.NR_Rounded1_Points + CD.NR_Rounded2_Points + CD.NR_ET_Points
            i_eff = i - CD.NR_RT_Points - CD.NR_Rounded1_Points - CD.NR_Rounded2_Points
            N = CD.NR_ET_Points
            push!(timePoints, t1(C) + t2(C) + t3(C) + t4(C) + t5(C) + i_eff*t6(C)/N)
        
        else
            i_eff = i - CD.NR_RT_Points - CD.NR_Rounded1_Points - CD.NR_Rounded2_Points - CD.NR_ET_Points
            N = CD.NR_Rounded3_Points
            push!(timePoints, t1(C) + t2(C) + t3(C) + t4(C) + t5(C) + t6(C) + i_eff*t7(C)/N)
        end

    end  
    return timePoints
end


################## Help functions to get nr of points ######################

function get_NR_FW_BW_points(CD::AbstractSK_ContourDiscretization)
    return CD.NR_RT_Points / 2
end

function get_NR_rounded1_points(CD::SchwingerKeldyshRounded1ContourDiscretization)
    return CD.NR_Rounded1_Points / 2
end

function get_total_points(CD::SchwingerKeldyshContourDiscretization)
    return CD.NR_RT_Points + CD.NR_ET_Points
end

function get_total_points(CD::SchwingerKeldyshRounded1ContourDiscretization)
    return CD.NR_RT_Points + CD.NR_Rounded1_Points + CD.NR_ET_Points
end

function get_total_points(CD::SchwingerKeldyshRounded2ContourDiscretization)
    return CD.NR_RT_Points + CD.NR_Rounded1_Points + CD.NR_Rounded2_Points + CD.NR_ET_Points
end

function get_total_points(CD::SchwingerKeldyshRounded3ContourDiscretization)
    return CD.NR_RT_Points + CD.NR_Rounded1_Points + CD.NR_Rounded2_Points + CD.NR_ET_Points + CD.NR_Rounded3_Points
end








