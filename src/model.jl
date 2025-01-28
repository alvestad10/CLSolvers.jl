export LM_HO_Param, LM_AHO_Param, AHO_Param, NonEquil_AHO_Param, Imag_Sig_Param, Imag_Sig_Lmd_Param,
    get_action, get_D_action, get_sde_funcs, getParamArray

#########################################################

####### MODEL ###########
abstract type AbstractModel end
abstract type AbstractOscillator <: AbstractModel end
abstract type Abstract_LM_Oscillator <: AbstractOscillator end

"""
    Imagenary sigma model parameters
"""
struct Imag_Sig_Param <: Abstract_LM_Oscillator
    σ::Complex
end

Base.string(p::Imag_Sig_Param) = string("σ: ", p.σ)

"""
    Imagenary sigma and lambda model parameters
"""
struct Imag_Sig_Lmd_Param <: Abstract_LM_Oscillator
    σ::Complex
    λ::Complex
end

Base.string(p::Imag_Sig_Lmd_Param) = string("σ: ", p.σ, "\n", "λ: ", p.λ)

"""
    Larg-Mass Harmonic oscillator parameters
"""
struct LM_HO_Param <: Abstract_LM_Oscillator
    m::Float64
end

Base.string(p::LM_HO_Param) = string("m: ", p.m)


"""
    Larg-Mass Anharmonic oscillator parameters
"""
struct LM_AHO_Param <: Abstract_LM_Oscillator
    m::Float64
    λ::Float64
end

Base.string(p::LM_AHO_Param) = string("m: ", p.m, "\n", "λ: ", p.λ)


"""
    Anharmonic oscillator parameters
"""
struct AHO_Param <: AbstractOscillator
    g::Float64
    m::Float64
    λ::Float64
end

Base.string(p::AHO_Param) = string("g: ", p.g, "\n", "m: ", p.m, "\n", "λ: ", p.λ)


"""
    Non-Equilibrium Anharmonic oscillator parameters
"""
struct NonEquil_AHO_Param <: AbstractOscillator
    g::Float64
    m::Float64
    λ::Float64

    # Source term parameters
    ϕ::Float64
    ξ::Float64
    σ::Float64
    ϕ_dot::Float64
    η::Float64
end

NonEquil_AHO_Param(m::Float64, λ::Float64, ϕ::Float64, ξ::Float64, σ::Float64, η::Float64) = NonEquil_AHO_Param(1, m, λ, ϕ, ξ, σ, 0, η)
NonEquil_AHO_Param(m::Float64, λ::Float64, ϕ::Float64, ξ::Float64, σ::Float64) = NonEquil_AHO_Param(1, m, λ, ϕ, ξ, σ, 0, 0)
NonEquil_AHO_Param(m::Float64, λ::Float64) = NonEquil_AHO_Param(1, m, λ, 1, 1, 1, 0, 0)

function Base.string(p::NonEquil_AHO_Param)
    string("g: ", p.g, "\n", "m: ", p.m, "\n", "λ: ", p.λ, "\n",
        "Source_term: Gaussian\n",
        "ϕ: ", p.ϕ, "\n", "ξ: ", p.ξ, "\n", "ϕ_dot: ", p.ϕ_dot, "\n",
        "η: ", p.η, "\n", "σ: ", p.σ)
end




###################################
# CONSTRUCTING THE PARAMETER ARRAY
###################################
"""
    Constructing the parameter array

    Arguments:
        p: Parameters of model
        a: Contour points spacing
    
    Return:
        LVector of arguments
"""
function getParamArray(p::T, a::AbstractArray) where {T<:AbstractModel}
    LVector(p=p, a=a,
        vals=fill(zeros(length(a)), 6),
        as=(circshift(a, 1), circshift(a, -1), circshift(a, 2))
    )
end
#########################
####### IMAG SIGMA ######
#########################
"""
(σ1+iσ2)(x1+ix2) = σ1*x1-σ2*x2 + i(σ1*x2+σ2*x2)
"""
function a_Imag_Sig_real(du, u, param, t)

    p = param.p

    du[1] = -(real(p.σ) * u[1] - imag(p.σ) * u[2])
    du[2] = -(real(p.σ) * u[2] + imag(p.σ) * u[1])
end

function b_Imag_Sig_real(du, u, param, t)
    du[1] = sqrt(2)
end

function get_sde_funcs(::Imag_Sig_Param)
    return a_Imag_Sig_real, b_Imag_Sig_real, nothing
end

################################
####### IMAG SIGMA LAMBDA ######
################################
"""
(σ1+iσ2)(x1+ix2) = σ1*x1-σ2*x2 + i(σ1*x2+σ2*x2)
(x1+ix2)(x1+ix2)(x1+ix2) = (x1*x1-x2*x2 + i(x1*x2 + x2*x2))(x1+ix2) 
                         = x1^3 - 2*x1*x2^2 - x2^3 + i(x1^2 * x2 - x2^3 + x1^2 x2 + x1 x2^2)
"""
function a_Imag_Sig_Lmd_real(du, u, param, t)

    p = param.p

    du[1] = -(real(p.σ) * u[1] - imag(p.σ) * u[2]) - (1 / 6) * p.λ(u[1] * u[1] * u[1] - 2 * u[1] * u[2] * u[2] - u[2] * u[2] * u[2])
    du[2] = -(real(p.σ) * u[2] + imag(p.σ) * u[1]) - (1 / 6) * p.λ(2 * u[1] * u[1] * u[2] + u[1] * u[1] * u[2] - u[2] * u[2] * u[2])
end

function b_Imag_Sig_Lmd_real(du, u, param, t)
    du[1, 1] = sqrt(2)
end

function get_sde_funcs(::Imag_Sig_Lmd_Param)
    return a_Imag_Sig_Lmd_real, b_Imag_Sig_Lmd_real, nothing
end


####################################
#######  AHO Non-Equilibrium  ######
####################################
"""
    Returns the a term in the SDE for the AHO model
"""
function a_AHO_real_nonEquil(du, u, param, t)


    p = param.p
    a = param.a

    a_R = real(a)
    a_I = imag(a)

    xR, xI, xR_diff_i, xR_diff_pi, xI_diff_i, xI_diff_pi = param.vals
    a_m1, a_p1, _ = param.as

    a_m1_R = real(a_m1)
    a_p1_I = imag(a_p1)

    xR = @view u[1:div(end, 2)]
    xI = @view u[div(end, 2)+1:end]

    xR_diff_i = xR .- xR[vcat(end, 1:end-1)]#circshift(xR,1)
    xR_diff_pi = xR[vcat(2:end, 1)] .- xR
    xI_diff_i = xI .- xI[vcat(end, 1:end-1)]#circshift(xI,1)
    xI_diff_pi = xI[vcat(2:end, 1)] .- xI

    ####
    # Non-Equilibrium factors
    ####
    #ϕ = 1
    #σ = 1
    #ξ = 1
    nonEqFac = ((p.σ^2 + 1) / (8 * p.ξ^2))
    nonEqFac2 = -((p.σ^2 - 1) / (4 * p.ξ^2))
    nonEqFac3 = p.η / (2p.ξ)

    #########
    # i=0
    ########

    preFac_1 = 1 / (abs(a[1]))

    a_R_1 = a_R[1]
    a_I_1 = a_I[1]

    x_R_1 = xR[1]
    x_I_1 = xI[1]

    x_R_N = xR[end]
    x_I_N = xI[end]

    a2_1 = (a_R_1^2 + a_I_1^2)

    du[1] = preFac_1 * (
        (1 / a2_1) * (a_R_1 * xI_diff_pi[1] - a_I_1 * xR_diff_pi[1])
        + (1 / 12) * (a_I_1 * x_R_1 * (6p.m + p.λ * (-3 * x_I_1^2 + x_R_1^2)) + a_R_1 * x_I_1 * (6p.m^2 + p.λ * (-x_I_1^2 + 3 * x_R_1^2)))
        + (2 * nonEqFac * (p.ϕ - x_R_1) + nonEqFac2 * (p.ϕ - x_R_N) + 2nonEqFac3 * x_I_1) # Gaussian distribution initial condition
    )

    du[div(end, 2)+1] = preFac_1 * (
        -(1 / a2_1) * (a_I_1 * xI_diff_pi[1] + a_R_1 * xR_diff_pi[1])
        -
        (1 / 12) * (a_R_1 * x_R_1 * (6p.m^2 + p.λ * (-3 * x_I_1^2 + x_R_1^2)) - a_I_1 * x_I_1 * (6p.m^2 + p.λ * (-x_I_1^2 + 3 * x_R_1^2)))
        -
        (2 * nonEqFac * x_I_1 + nonEqFac2 * x_I_N - 2nonEqFac3 * (p.ϕ - x_R_1)) # Gaussian distribution initial condition
    )



    #########
    # i=N
    ########
    preFac_N = 1 / (abs(a[end]))

    a_R_Nm1 = a_R[end]
    a_I_Nm1 = a_I[end]


    a2_Nm1 = (a_R_Nm1^2 + a_I_Nm1^2)

    du[div(end, 2)] = preFac_N * (
        (1 / a2_Nm1) * (a_I_Nm1 * xR_diff_i[end] - a_R_Nm1 * xI_diff_i[end]) #* - a_R_Nm1*xI_diff_i[end] )
        + (1 / 12) * (a_I_Nm1 * x_R_N * (6p.m + p.λ * (-3 * x_I_N^2 + x_R_N^2)) + a_R_Nm1 * x_I_N * (6p.m^2 + p.λ * (-x_I_N^2 + 3 * x_R_N^2)))
        + (2 * nonEqFac * (p.ϕ - x_R_N) + nonEqFac2 * (p.ϕ - x_R_1) - 2nonEqFac3 * x_I_N) # Gaussian distribution initial condition
    )

    du[end] = preFac_N * (
        (1 / a2_Nm1) * (a_I_Nm1 * xI_diff_i[end] + a_R_Nm1 * xR_diff_i[end])
        -
        (1 / 12) * (a_R_Nm1 * x_R_N * (6p.m + p.λ * (-3 * x_I_N^2 + x_R_N^2)) - a_I_Nm1 * x_I_N * (6p.m^2 + p.λ * (-x_I_N^2 + 3 * x_R_N^2)))
        -
        (2 * nonEqFac * x_I_N + nonEqFac2 * x_I_1 + 2nonEqFac3 * (p.ϕ - x_R_N))# Gaussian distribution initial condition
    )

    ###############
    # i=1 to i=N-1
    ###############
    preFac = 1 ./ (0.5 .* (abs.(a[2:end]) .+ abs.(a_m1[2:end])))

    a_R_j = a_R[2:end]
    a_I_j = a_I[2:end]
    a_R_jm1 = a_R[1:end-1]
    a_I_jm1 = a_I[1:end-1]

    x_R_j = xR[2:end-1]
    x_I_j = xI[2:end-1]
    xR_diff_j = xR_diff_i[2:end-1]
    xR_diff_pj = xR_diff_pi[2:end-1]
    xI_diff_j = xI_diff_i[2:end-1]
    xI_diff_pj = xI_diff_pi[2:end-1]



    du[2:div(end, 2)-1] .= preFac .* (
        1 .* (a_I_jm1 .* xR_diff_j .- a_R_jm1 .* xI_diff_j) ./ (a_I_jm1 .^ 2 .+ a_R_jm1 .^ 2)
        .+
        1 .* (a_R_j .* xI_diff_pj .- a_I_j .* xR_diff_pj) ./ (a_I_j .^ 2 .+ a_R_j .^ 2)
        .+
        (1.0 / 12.0) .* (
            (a_I_jm1 .+ a_I_j) .* x_R_j .* (6 * p.m^2 .+ p.λ .* (.-3 .* x_I_j .^ 2 .+ x_R_j .^ 2))
            .+
            (a_R_jm1 .+ a_R_j) .* x_I_j .* (6 * p.m^2 .+ p.λ .* (.-x_I_j .^ 2 .+ 3 .* x_R_j .^ 2))
        )
    )

    du[div(end, 2)+2:end-1] .= preFac .* (
        1 .* (a_I_jm1 .* xI_diff_j .+ a_R_jm1 .* xR_diff_j) ./ (a_I_jm1 .^ 2 .+ a_R_jm1 .^ 2)
        .-
        1 .* (a_I_j .* xI_diff_pj .+ a_R_j .* xR_diff_pj) ./ (a_I_j .^ 2 .+ a_R_j .^ 2)
        .+
        (1.0 / 12.0) .* (
            (a_I_jm1 .+ a_I_j) .* x_I_j .* (6 * p.m^2 .+ p.λ .* (.-x_I_j .^ 2 .+ 3 .* x_R_j .^ 2))
            .-
            (a_R_jm1 .+ a_R_j) .* x_R_j .* (6 * p.m^2 .+ p.λ .* (.-3 .* x_I_j .^ 2 .+ x_R_j .^ 2))
        )
    )


end

function b_AHO_real_nonEquil(du, u, param, t)
    a = param.a

    du[1] = sqrt(2) / sqrt(abs(a[1]))
    du[2:div(end, 2)-1] .= sqrt.(2 ./ abs.(a[2:end]))
    du[div(end, 2)] = sqrt(2) / sqrt(abs(a[end]))

end

function get_sde_funcs(::NonEquil_AHO_Param)
    return a_AHO_real_nonEquil, b_AHO_real_nonEquil, nothing
end


##################
####### AHO ######
##################
"""
    Returns the a term in the SDE for the AHO model
"""
function a_AHO_real(du, u, param, t)


    p = param.p
    a = param.a

    xR, xI, xR_diff_i, xR_diff_pi, xI_diff_i, xI_diff_pi = param.vals
    a_m1, a_p1, _ = param.as

    xR = @view u[1:div(end, 2)]
    xI = @view u[div(end, 2)+1:end]

    xR_diff_i = xR .- xR[vcat(end, 1:end-1)]#circshift(xR,1)
    xR_diff_pi = xR[vcat(2:end, 1)] .- xR
    xI_diff_i = xI .- xI[vcat(end, 1:end-1)]#circshift(xI,1)
    xI_diff_pi = xI[vcat(2:end, 1)] .- xI

    preFac = (1.0 ./ (0.5 .* (abs.(a) + abs.(a_m1)))) .* 0.5

    du[1:div(end, 2)] .= -preFac .* (
        2.0 .* p.g .* (real.(a_m1) .* xI_diff_i .- imag.(a_m1) .* xR_diff_i) ./ abs.(a_m1) .^ 2
        .-
        2.0 .* p.g .* (real.(a) .* xI_diff_pi .- imag.(a) .* xR_diff_pi) ./ abs.(a) .^ 2
        .-
        real.(a_m1 .+ a) .* (p.m .* xI .+ (1 / 6) * p.λ .* (-(xI .^ 3) .+ 3xI .* (xR .^ 2)))
        .-
        imag.(a_m1 .+ a) .* (p.m .* xR .+ (1 / 6) * p.λ .* (xR .^ 3 .- 3xR .* (xI .^ 2)))
        #.+ 1e-8
    )

    du[div(end, 2)+1:end] .= preFac .* (
        2.0 .* p.g .* (real.(a_m1) .* xR_diff_i .+ imag.(a_m1) .* xI_diff_i) ./ abs.(a_m1) .^ 2
        .-
        2.0 .* p.g .* (real.(a) .* xR_diff_pi .+ imag.(a) .* xI_diff_pi) ./ abs.(a) .^ 2
        .-
        real.(a_m1 .+ a) .* (p.m .* xR .+ (1 / 6) * p.λ .* ((xR .^ 3) .- 3xR .* (xI .^ 2)))
        .+
        imag.(a_m1 .+ a) .* (p.m .* xI .+ (1 / 6) * p.λ .* (-xI .^ 3 .+ 3xI .* (xR .^ 2)))
    )
end


"""
    Returns the jacobian of the a term in the SDE for the AHO model
"""
function jac_AHO_real(J, u, param, t)

    p = param.p
    a = param.a


    t_steps = div(length(u), 2)

    xR = u[1:t_steps]
    xI = u[t_steps+1:end]
    a_m1 = circshift(a, 1)
    a_p1 = circshift(a, -1)



    for i in 1:t_steps
        preFac = (1.0 / (0.5 * (abs(a[i]) + abs(a_m1[i])))) * 0.5
        J[i, i] = -preFac * (
            -2.0 * p.g * imag(a_m1[i]) / abs(a_m1[i])^2
            -
            2.0 * p.g * imag(a[i]) / abs(a[i])^2
            -
            real(a_p1[i] + a[i]) * ((1 / 6) * p.λ * (6 * xI[i] * xR[i]))
            -
            imag(a_p1[i] + a[i]) * (p.m + (1 / 6) * p.λ * (3 * (xR[i])^2 - 3 * (xI[i]) .^ 2))
        )

        J[t_steps+i, t_steps+i] = preFac * (
            2.0 * p.g * imag(a_m1[i]) / abs(a_m1[i])^2
            + 2.0 * p.g * imag(a[i]) / abs(a[i])^2
            + real(a_p1[i] + a[i]) * (1 / 6) * p.λ * (6 * xR[i] * xI[i])
            + imag(a_p1[i] + a[i]) * (p.m + (1 / 6) * p.λ * (-3 * xI[i]^2 + 3(xR[i]^2)))
        )

        J[i, t_steps+i] = -preFac * (
            2.0 * p.g * real(a_m1[i]) / abs(a_m1[i])^2
            +
            2.0 * p.g * real(a[i]) / abs(a[i])^2
            -
            real(a_p1[i] + a[i]) * (p.m + (1 / 6) * p.λ * (-3(xI[i]^2) + 3(xR[i]^2)))
            .-
            imag(a_p1[i] + a[i]) * ((1 / 6) * p.λ * (6 * xR[i] * xI[i]))
        )

        J[t_steps+i, i] = preFac * (
            2.0 * p.g * real(a_m1[i]) / abs(a_m1[i])^2
            +
            2.0 * p.g * real(a[i]) / abs(a[i])^2
            -
            real(a_p1[i] + a[i]) * (p.m + (1 / 6) * p.λ * 3 * ((xR[i]^2) - 3 * (xI[i]^2)))
            +
            imag(a_p1[i] + a[i]) * (1 / 6) * p.λ * (6 * xI[i] * xR[i])
        )

        J[i, mod1(i + 1, t_steps)] = -preFac * (
            +2.0 * p.g * imag(a[i]) / abs(a[i])^2
        )

        J[i, t_steps+mod1(i + 1, t_steps)] = -preFac * (
            -2.0 * p.g * real(a[i]) / abs(a[i])^2
        )

        J[i, mod1(i - 1, t_steps)] = -preFac * (
            2.0 * p.g * imag(a_m1[i]) / abs(a_m1[i])^2
        )

        J[i, t_steps+mod1(i - 1, t_steps)] = -preFac * (
            -2.0 * p.g * real(a_m1[i]) / abs(a_m1[i])^2
        )

        J[t_steps+i, t_steps+mod1(i + 1, t_steps)] = preFac * (
            -2.0 * p.g * imag(a[i]) / abs(a[i])^2
        )

        J[t_steps+i, mod1(i + 1, t_steps)] = preFac * (
            -2.0 * p.g * real(a[i]) / abs(a[i])^2
        )

        J[t_steps+i, t_steps+mod1(i - 1, t_steps)] = preFac * (
            -2.0 * p.g * imag(a_m1[i]) / abs(a_m1[i])^2
        )

        J[t_steps+i, mod1(i - 1, t_steps)] = preFac * (
            -2.0 * p.g * real(a_m1[i]) / abs(a_m1[i])^2
        )

    end
end

"""
    Returns the b term in the SDE for the AHO model
"""
function b_AHO_real(du, u, param, t)
    a = param.a
    a_m1, a_p1, _ = param.as
    du[1:div(end, 2)] .= sqrt.(2 ./ (0.5 .* (abs.(a) + abs.(a_m1))))
end

function get_sde_funcs(::AHO_Param)
    return a_AHO_real, b_AHO_real, jac_AHO_real
end