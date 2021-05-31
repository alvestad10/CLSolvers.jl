export NonEquilSchwingerKeldyshContour,
    NonEquilTiltedSchwingerKeldyshContour,
    SchwingerKeldyshContour,
    SchwingerKeldyshContourRounded1,
    SchwingerKeldyshContourRounded2,
    SchwingerKeldyshContourRounded3,
    get_timepoint,
    get_derivative_timepoint,
    EucledianContour
       
################################################################

abstract type AbstractContour end
abstract type AbstractSchwingerKeldyshContour <: AbstractContour end

abstract type AbstractNonEquilSchwingerKeldyshContour <: AbstractContour end


####################################################
###         DEFINING the contour objects
####################################################

################ Eucledian ###################

struct EucledianContour <: AbstractContour
    β
end

"""
    Calculate the time point on the contour from the parametetric variable u in [0,1]
    for the Schwinger-Keldysh contour
"""
function get_timepoint(u,C::EucledianContour)
    return -u * C.β * im
end



############### Non-Equilibrium Schwinger-Keldysh contour ###########
struct NonEquilSchwingerKeldyshContour <: AbstractNonEquilSchwingerKeldyshContour
    tmax::Float32
    ϵ::Float64
end

Base.string(c::NonEquilSchwingerKeldyshContour) = string("tmax: ", c.tmax, "\n", "ϵ: ", c.ϵ)

L1(contour::NonEquilSchwingerKeldyshContour) = contour.tmax
L2(contour::NonEquilSchwingerKeldyshContour) = contour.tmax
T(contour::AbstractNonEquilSchwingerKeldyshContour) = L1(contour) + L2(contour)
t1(contour::AbstractNonEquilSchwingerKeldyshContour) = L1(contour) / T(contour)
t2(contour::AbstractNonEquilSchwingerKeldyshContour) = L2(contour) / T(contour)

############### Non-Equilibrium tilted Schwinger-Keldysh contour ###########
struct NonEquilTiltedSchwingerKeldyshContour <: AbstractNonEquilSchwingerKeldyshContour
    tmax::Float32
    θ::Float64      # Tilt angel
end

Base.string(c::NonEquilTiltedSchwingerKeldyshContour) = print("tmax: ", c.tmax, "\n", "θ: ", c.θ)

L1(contour::NonEquilTiltedSchwingerKeldyshContour) = contour.tmax/cos(contour.θ)
L2(contour::NonEquilTiltedSchwingerKeldyshContour) = contour.tmax/cos(contour.θ)




############### Schwinger-Keldysh contour ###############

struct SchwingerKeldyshContour <: AbstractSchwingerKeldyshContour
    β::Float32
    tmax::Float32
    A_β::Float64
    F_β::Float64
end

Base.string(c::SchwingerKeldyshContour) = string("tmax: ", c.tmax, "\n", "β: ", c.β, "\n", "A_β: ", c.A_β, "\n", "F_β: ", c.F_β)


# Defining functions related to the calculation of the SchwingerKeldyshContour
β1(contour::AbstractSchwingerKeldyshContour) = contour.A_β * contour.F_β * contour.β
β2(contour::AbstractSchwingerKeldyshContour) = (1-contour.A_β) * contour.F_β * contour.β
L1(contour::SchwingerKeldyshContour) = sqrt(β1(contour)^2 + contour.tmax^2)
L2(contour::SchwingerKeldyshContour) = sqrt(β2(contour)^2 + contour.tmax^2)
L3(contour::SchwingerKeldyshContour) = (1-contour.F_β)*contour.β
T(contour::SchwingerKeldyshContour)  = L1(contour) + L2(contour) + L3(contour)
t1(contour::SchwingerKeldyshContour) = L1(contour) / T(contour)
t2(contour::SchwingerKeldyshContour) = L2(contour) / T(contour)
t3(contour::SchwingerKeldyshContour) = L3(contour) / T(contour)

# Default parameters 
SchwingerKeldyshContour() = 
    SchwingerKeldyshContour(1.0,0.8,0.01,0.5) 

############### Schwinger-Keldysh contour Rounded 1 ###############
#
#   The Schwinger-Keldysh contour for one rounded corner, i.e.,
#   rounded corner at the $t_{max}$ point.
#   
#   
#   NR_RT_Points, NR_ET_Points, β, t_max, A_β and F_β is the same parameters
#   as in normal Schwinger-Keldysh contour
#   
#   NR_Rounded1_Points: nr of points in the arc
#   r1: radius of the rounded corner
#   
#
struct SchwingerKeldyshContourRounded1 <: AbstractSchwingerKeldyshContour
    β::Float32
    tmax::Float32
    A_β::Float64
    F_β::Float64
    r1::Float64
end

Base.string(c::SchwingerKeldyshContourRounded1) = string("tmax: ", c.tmax, "\n", "β: ", c.β, "\n", "A_β: ", x.A_β, "\n", "F_β:", c.F_β, "\n", "r1:", c.r1)


# Defining functions related to the calculation of the SchwingerKeldyshContour
#β1(contour::SchwingerKeldyshContour) = contour.A_β * contour.F_β * contour.β
#β2(contour::SchwingerKeldyshContourRounded1) = (1-contour.A_β) * contour.F_β * contour.β
β3(contour::AbstractSchwingerKeldyshContour) = (1-contour.F_β) * contour.β

H1(contour::SchwingerKeldyshContourRounded1) = sqrt((contour.tmax-contour.r1)^2 + β1(contour)^2)
L1(contour::SchwingerKeldyshContourRounded1) = sqrt(H1(contour)^2 - contour.r1^2)
γ1(contour::SchwingerKeldyshContourRounded1) = pi/2 - asin(contour.r1/H1(contour)) - acos(β1(contour)/H1(contour))

H2(contour::SchwingerKeldyshContourRounded1) = sqrt((contour.tmax-contour.r1)^2 + β2(contour)^2)
L2(contour::SchwingerKeldyshContourRounded1) = sqrt(H2(contour)^2 - contour.r1^2)
γ2(contour::SchwingerKeldyshContourRounded1) = pi/2 - asin(contour.r1/H2(contour)) - acos(β2(contour)/H2(contour))

Gy1(contour::SchwingerKeldyshContourRounded1) = L1(contour)*sin(γ1(contour))
Gx1(contour::SchwingerKeldyshContourRounded1) = L1(contour)*cos(γ1(contour))
Gy2(contour::SchwingerKeldyshContourRounded1) = L2(contour)*sin(γ2(contour))
Gx2(contour::SchwingerKeldyshContourRounded1) = L2(contour)*cos(γ2(contour))

η1(contour::SchwingerKeldyshContourRounded1) = pi/2 - γ1(contour)
η2(contour::SchwingerKeldyshContourRounded1) = pi/2 - γ2(contour)

CC1(contour::SchwingerKeldyshContourRounded1) = η1(contour)*contour.r1
CC2(contour::SchwingerKeldyshContourRounded1) = η2(contour)*contour.r1
L3(contour::SchwingerKeldyshContourRounded1) = (1-contour.F_β)*contour.β

T(contour::SchwingerKeldyshContourRounded1)  = L1(contour) + CC1(contour) + CC2(contour) + L2(contour) + L3(contour)
t1(contour::SchwingerKeldyshContourRounded1) = L1(contour) / T(contour)
t2(contour::SchwingerKeldyshContourRounded1) = CC1(contour) / T(contour)
t3(contour::SchwingerKeldyshContourRounded1) = CC2(contour) / T(contour)
t4(contour::SchwingerKeldyshContourRounded1) = L2(contour) / T(contour)
t5(contour::SchwingerKeldyshContourRounded1) = L3(contour) / T(contour)


############### Schwinger-Keldysh contour Rounded 2 ###############
#
#   The Schwinger-Keldysh contour for one rounded corner, i.e.,
#   rounded corner at the $t_{max}$ point.
#   
#   
#   NR_RT_Points, NR_ET_Points, β, t_max, A_β and F_β is the same parameters
#   as in normal Schwinger-Keldysh contour
#   
#   NR_Rounded1_Points: nr of points in the arc
#   r1: radius of the rounded corner
#   
#
struct SchwingerKeldyshContourRounded2 <: AbstractSchwingerKeldyshContour
    β::Float32
    tmax::Float32
    A_β::Float64
    F_β::Float64
    r1::Float64
    r2::Float64
end

Base.string(c::SchwingerKeldyshContourRounded2) = string("tmax: ", c.tmax, "\n", "β: ", c.β, "\n", "A_β: ", x.A_β, "\n", "F_β:", c.F_β, "\n",
                                                              "r1:", c.r1,"\n", "r2:", c.r2)


# Defining functions related to the calculation of the SchwingerKeldyshContour
#β1(contour::SchwingerKeldyshContourRounded2) = contour.A_β * contour.F_β * contour.β
#β2(contour::SchwingerKeldyshContourRounded2) = (1-contour.A_β) * contour.F_β * contour.β
#β3(contour::SchwingerKeldyshContourRounded2) = (1-contour.F_β) * contour.β

H1(contour::SchwingerKeldyshContourRounded2) = sqrt((contour.tmax-contour.r1)^2 + β1(contour)^2)
L1(contour::SchwingerKeldyshContourRounded2) = sqrt(H1(contour)^2 - contour.r1^2)
γ1(contour::SchwingerKeldyshContourRounded2) = pi/2 - asin(contour.r1/H1(contour)) - acos(β1(contour)/H1(contour))

H2(contour::SchwingerKeldyshContourRounded2) = sqrt((contour.tmax-contour.r1-contour.r2)^2 + β2(contour)^2)
L2(contour::SchwingerKeldyshContourRounded2) = sqrt(H2(contour)^2 - (contour.r1+contour.r2)^2)
γ2(contour::SchwingerKeldyshContourRounded2) = pi/2 - asin((contour.r1+contour.r2)/H2(contour)) - acos(β2(contour)/H2(contour))

Gy1(contour::SchwingerKeldyshContourRounded2) = L1(contour)*sin(γ1(contour))
Gx1(contour::SchwingerKeldyshContourRounded2) = L1(contour)*cos(γ1(contour))
Gy2(contour::SchwingerKeldyshContourRounded2) = L2(contour)*sin(γ2(contour))
Gx2(contour::SchwingerKeldyshContourRounded2) = L2(contour)*cos(γ2(contour))

η1(contour::SchwingerKeldyshContourRounded2) = pi/2 - γ1(contour)
η2(contour::SchwingerKeldyshContourRounded2) = pi/2 - γ2(contour)

CC1_1(contour::SchwingerKeldyshContourRounded2) = η1(contour)*contour.r1
CC1_2(contour::SchwingerKeldyshContourRounded2) = η2(contour)*contour.r1
CC2(contour::SchwingerKeldyshContourRounded2) = η2(contour)*contour.r2
L3(contour::SchwingerKeldyshContourRounded2) = β3(contour)

T(contour::SchwingerKeldyshContourRounded2)  = L1(contour) + CC1_1(contour) + CC1_2(contour) + L2(contour) + CC2(contour) + L3(contour)
t1(contour::SchwingerKeldyshContourRounded2) = L1(contour) / T(contour)
t2(contour::SchwingerKeldyshContourRounded2) = CC1_1(contour) / T(contour)
t3(contour::SchwingerKeldyshContourRounded2) = CC1_2(contour) / T(contour)
t4(contour::SchwingerKeldyshContourRounded2) = L2(contour) / T(contour)
t5(contour::SchwingerKeldyshContourRounded2) = CC2(contour) / T(contour)
t6(contour::SchwingerKeldyshContourRounded2) = L3(contour) / T(contour)

############### Schwinger-Keldysh contour Rounded 3 ###############
#
#   The Schwinger-Keldysh contour for one rounded corner, i.e.,
#   rounded corner at the $t_{max}$ point.
#   
#   
#   NR_RT_Points, NR_ET_Points, β, t_max, A_β and F_β is the same parameters
#   as in normal Schwinger-Keldysh contour
#   
#   NR_Rounded1_Points: nr of points in the arc
#   r1: radius of the rounded corner
#   
#
struct SchwingerKeldyshContourRounded3 <: AbstractSchwingerKeldyshContour
    β::Float32
    tmax::Float32
    A_β::Float64
    F_β::Float64
    r1::Float64
    r2::Float64
    r3::Float64
end

Base.string(io::IO, c::SchwingerKeldyshContourRounded3) = string("tmax: ", c.tmax, "\n", "β: ", c.β, "\n", "A_β: ", x.A_β, "\n", "F_β:", c.F_β, "\n",
                                                              "r1:", c.r1, "\n", "r2:", c.r2, "\n", "r3:", c.r3)

# Defining functions related to the calculation of the SchwingerKeldyshContour
#β1(contour::SchwingerKeldyshContourRounded2) = contour.A_β * contour.F_β * contour.β
#β2(contour::SchwingerKeldyshContourRounded2) = (1-contour.A_β) * contour.F_β * contour.β
#β3(contour::SchwingerKeldyshContourRounded2) = (1-contour.F_β) * contour.β

H1(contour::SchwingerKeldyshContourRounded3) = sqrt((contour.tmax-contour.r1-contour.r3)^2 + (β1(contour)+contour.r3)^2)
L1(contour::SchwingerKeldyshContourRounded3) = sqrt(H1(contour)^2 - (contour.r1+contour.r3)^2)
γ1(contour::SchwingerKeldyshContourRounded3) = pi/2 - asin((contour.r1+contour.r3)/H1(contour)) - acos((β1(contour)+contour.r3)/H1(contour))

H2(contour::SchwingerKeldyshContourRounded3) = sqrt((contour.tmax-contour.r1-contour.r2)^2 + β2(contour)^2)
L2(contour::SchwingerKeldyshContourRounded3) = sqrt(H2(contour)^2 - (contour.r1+contour.r2)^2)
γ2(contour::SchwingerKeldyshContourRounded3) = pi/2 - asin((contour.r1+contour.r2)/H2(contour)) - acos(β2(contour)/H2(contour))

Gy1(contour::SchwingerKeldyshContourRounded3) = L1(contour)*sin(γ1(contour))
Gx1(contour::SchwingerKeldyshContourRounded3) = L1(contour)*cos(γ1(contour))
Gy2(contour::SchwingerKeldyshContourRounded3) = L2(contour)*sin(γ2(contour))
Gx2(contour::SchwingerKeldyshContourRounded3) = L2(contour)*cos(γ2(contour))

η1(contour::SchwingerKeldyshContourRounded3) = pi/2 - γ1(contour)
η2(contour::SchwingerKeldyshContourRounded3) = pi/2 - γ2(contour)

CC1_1(contour::SchwingerKeldyshContourRounded3) = η1(contour)*contour.r1
CC1_2(contour::SchwingerKeldyshContourRounded3) = η2(contour)*contour.r1
CC2(contour::SchwingerKeldyshContourRounded3) = η2(contour)*contour.r2

L3(contour::SchwingerKeldyshContourRounded3) = β3(contour) - contour.r3*sin(η1(contour))
CC3(contour::SchwingerKeldyshContourRounded3) = η1(contour)*contour.r3

T(contour::SchwingerKeldyshContourRounded3)  = L1(contour) + CC1_1(contour) + CC1_2(contour) + L2(contour) + CC2(contour) + L3(contour) + CC3(contour)
t1(contour::SchwingerKeldyshContourRounded3) = L1(contour) / T(contour)
t2(contour::SchwingerKeldyshContourRounded3) = CC1_1(contour) / T(contour)
t3(contour::SchwingerKeldyshContourRounded3) = CC1_2(contour) / T(contour)
t4(contour::SchwingerKeldyshContourRounded3) = L2(contour) / T(contour)
t5(contour::SchwingerKeldyshContourRounded3) = CC2(contour) / T(contour)
t6(contour::SchwingerKeldyshContourRounded3) = L3(contour) / T(contour)
t7(contour::SchwingerKeldyshContourRounded3) = CC3(contour) / T(contour)


################################################
### Defining the related functions
################################################

"""
    Calculate the time point on the contour from the parametetric variable u in [0,1]
    for the Non-Equilibrium Schwinger-Keldysh contour
"""
function get_timepoint(u,C::NonEquilSchwingerKeldyshContour)
    if u <= t1(C)
        return complex(C.tmax*(u/t1(C)), 0)
    elseif t1(C) < u <= t1(C) + t2(C)
        return complex(C.tmax - C.tmax*((u-t1(C))/t2(C)))
    end
    return nothing
end


"""
    Calculate the time point on the contour from the parametetric variable u in [0,1]
    for the Non-Equilibrium tilted Schwinger-Keldysh contour
"""
function get_timepoint(u,C::NonEquilTiltedSchwingerKeldyshContour)
    if u <= t1(C)
        return complex(C.tmax*(u/t1(C)), -(C.tmax * C.θ)*(u/t1(C)))
    elseif t1(C) < u <= t1(C) + t2(C)
        return complex(C.tmax - C.tmax*((u-t1(C))/t2(C)),-(C.tmax * C.θ) -  (C.tmax * C.θ)* ((u-t1(C))/t2(C)))
    end
    return nothing
end


"""
    Calculate the time point on the contour from the parametetric variable u in [0,1]
    for the Schwinger-Keldysh contour
"""
function get_timepoint(u,C::SchwingerKeldyshContour)
    if u <= t1(C)
        return complex(C.tmax*(u/t1(C)), -β1(C)*(u/t1(C)))
    elseif t1(C) < u <= t1(C) + t2(C)
        return complex(C.tmax - C.tmax*((u-t1(C))/t2(C)), -β1(C) - β2(C)*((u-t1(C))/t2(C)))
    elseif u > t1(C) + t2(C)
        return complex(0, -β1(C) - β2(C) - L3(C)*((u-t1(C)-t2(C))/t3(C)))
    end
    return nothing
end

"""
    Calculate the time point on the contour from the parametetric variable u in [0,1]
    for the Schwinger-Keldysh contour
"""
function get_timepoint(u,C::SchwingerKeldyshContourRounded1)
    if (u <= t1(C))
        ui = u / t1(C)
        return complex(Gx1(C)*ui, - Gy1(C)*ui);
    elseif t1(C) < u <= t1(C) + t2(C)
        ui = (u - t1(C))/t2(C)
        return complex(C.tmax - C.r1 + C.r1*cos(η1(C) - η1(C)*ui), -β1(C) + C.r1*sin( η1(C) - η1(C)*ui));
    elseif t1(C) + t2(C) < u <= t1(C) + t2(C) + t3(C)
        ui = (u - t1(C) - t2(C))/t3(C)
        return complex(C.tmax - C.r1 + C.r1*cos(-η2(C)*ui), -β1(C)+C.r1*sin(-η2(C)*ui));
    elseif t1(C) + t2(C) + t3(C) < u <= t1(C) + t2(C) + t3(C) + t4(C)
        ui = (u - t1(C) - t2(C) - t3(C))/t4(C)
        return complex(Gx2(C) - Gx2(C)*ui, - (β1(C)+β2(C) - Gy2(C)) - Gy2(C)*ui);
    else
        ui = (u - t1(C) - t2(C) - t3(C) - t4(C))/t5(C)
        return complex(0, - β1(C)-β2(C) - ui*β3(C))
    end

    return nothing
end

"""
    Calculate the time point on the contour from the parametetric variable u in [0,1]
    for the Schwinger-Keldysh contour
"""
function get_timepoint(u,C::SchwingerKeldyshContourRounded2)
    if (u <= t1(C))
        ui = u / t1(C)
        return complex(Gx1(C)*ui, - Gy1(C)*ui);
    elseif t1(C) < u <= t1(C) + t2(C)
        ui = (u - t1(C))/t2(C)
        return complex(C.tmax - C.r1 + C.r1*cos(η1(C) - η1(C)*ui), -β1(C) + C.r1*sin( η1(C) - η1(C)*ui));
    elseif t1(C) + t2(C) < u <= t1(C) + t2(C) + t3(C)
        ui = (u - t1(C) - t2(C))/t3(C)
        return complex(C.tmax - C.r1 + C.r1*cos(-η2(C)*ui), -β1(C)+C.r1*sin(-η2(C)*ui));
    elseif t1(C) + t2(C) + t3(C) < u <= t1(C) + t2(C) + t3(C) + t4(C)
        ui = (u - t1(C) - t2(C) - t3(C))/t4(C)
        return complex(Gx2(C) + C.r2*(1-cos(pi/2-γ2(C))) - Gx2(C)*ui, - (β1(C)+β2(C) - Gy2(C) - C.r2*sin(pi/2 -γ2(C))) - Gy2(C)*ui);
    elseif t1(C) + t2(C) + t3(C) + t4(C) < u <= t1(C) + t2(C) + t3(C) + t4(C) + t5(C)
        ui = (u - t1(C) - t2(C) - t3(C) - t4(C))/t5(C)
        return complex(C.r2*(1 + cos(pi-η2(C)*(1-ui))), -β1(C)-β2(C) + C.r2*sin(pi-η2(C)*(1-ui)) ); 
    else
        ui = (u - t1(C) - t2(C) - t3(C) - t4(C) - t5(C))/t6(C)
        return complex(0, - β1(C)-β2(C) - ui*β3(C))
    end

    return nothing
end

"""
    Calculate the time point on the contour from the parametetric variable u in [0,1]
    for the Schwinger-Keldysh contour with three rounded corners
"""
function get_timepoint(u,C::SchwingerKeldyshContourRounded3)
    if (u <= t1(C))
        ui = u / t1(C)
        return complex(C.r3*(1-cos(η1(C))) + Gx1(C)*ui, C.r3*(1-sin(η1(C))) - Gy1(C)*ui);
    elseif t1(C) < u <= t1(C) + t2(C)
        ui = (u - t1(C))/t2(C)
        return complex(C.tmax - C.r1 + C.r1*cos(η1(C) - η1(C)*ui), -β1(C) + C.r1*sin( η1(C) - η1(C)*ui));
    elseif t1(C) + t2(C) < u <= t1(C) + t2(C) + t3(C)
        ui = (u - t1(C) - t2(C))/t3(C)
        return complex(C.tmax - C.r1 + C.r1*cos(-η2(C)*ui), -β1(C)+C.r1*sin(-η2(C)*ui));
    elseif t1(C) + t2(C) + t3(C) < u <= t1(C) + t2(C) + t3(C) + t4(C)
        ui = (u - t1(C) - t2(C) - t3(C))/t4(C)
        return complex(Gx2(C) + C.r2*(1-cos(pi/2-γ2(C))) - Gx2(C)*ui, - (β1(C)+β2(C) - Gy2(C) - C.r2*sin(pi/2 -γ2(C))) - Gy2(C)*ui);
    elseif t1(C) + t2(C) + t3(C) + t4(C) < u <= t1(C) + t2(C) + t3(C) + t4(C) + t5(C)
        ui = (u - t1(C) - t2(C) - t3(C) - t4(C))/t5(C)
        return complex(C.r2*(1 + cos(pi-η2(C)*(1-ui))), -β1(C)-β2(C) + C.r2*sin(pi-η2(C)*(1-ui)) ); 
    elseif t1(C) + t2(C) + t3(C) + t4(C) + t5(C) < u <= t1(C) + t2(C) + t3(C) + t4(C) + t5(C) + t6(C)
        ui = (u - t1(C) - t2(C) - t3(C) - t4(C) - t5(C))/t6(C)
        return complex(0, - β1(C)-β2(C) - ui*L3(C))
    else
        ui = (u - t1(C) - t2(C) - t3(C) - t4(C) - t5(C) - t6(C))/t7(C)
        return complex(C.r3*(1+cos(pi + η1(C)*ui)), -C.β + C.r3*( 1 + sin(pi + η1(C)*ui) ))
    end

    return nothing
end



##################### Calculating the derivative ############################

"""
    Caluclate the derivative of the timepoint in the contour
    given the parameter variable u in [0,1]

    The derivative is taken using the ForwardDiff.derivative

    This will work for all contours
"""
function get_derivative_timepoint(u,C::AbstractContour)
    return complex(
            ForwardDiff.derivative(x->real(get_timepoint(x,C)),u),
            ForwardDiff.derivative(x->imag(get_timepoint(x,C)),u)
        )
end

