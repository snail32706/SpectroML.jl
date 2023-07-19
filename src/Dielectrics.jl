module Dielectics

# export LestSquare

# = = = = = = = = = = = = = = = = = = = = = #
# Loren2z Model                             #
# = = = = = = = = = = = = = = = = = = = = = #

#=
ϵr12 function:
--------------
	c1 + (c4 - c2 * ω^2) / ((c3 - ω^2)^2 + c5 * ω^2) + i (c6 * ω) / ((c3 - ω^2)^2 + c5 * ω^2)

Parameters:
--------------
c0 ~ c5 are converted from the parameters of the Lorentz Model.

Returns:
--------------
Tuple{Float64, Float64}. Representing real part and image part of dielectric constant

Examples:
--------------
params = (ϵ_∞, ωₚ, fⱼ, Ωⱼ, γj)
ϵr12(3.77e15, ϵrParams(params...)...)
> (4.32577, 1.21313e-16)
=#
function ϵr12(ω::Float64, ϵinf::Float64, c1::Float64, c2::Float64, c3::Float64, c4::Float64, c5::Float64)
    return ϵinf + (c3 - c1 * abs2(ω)) / ( abs2(c2 -  abs2(ω)) + c4 * abs2(ω)), 
                             (c5 * ω) / ( abs2(c2 -  abs2(ω)) + c4 * abs2(ω))
end

function ϵrParams(ϵinf::Real, ωp::Real, fj::Real, Ωj::Real, γj::Real)
	arg1 = ϵinf
	arg2 = fj * abs2(ωp)
	arg3 = abs2(Ωj)
	arg4 = arg2 * arg3
	arg5 = abs2(γj)
	arg6 = arg2 * γj
	return (arg1, arg2, arg3, arg4, arg5, arg6)
end

# = = = = = = = = = = = = = = = = = = = = = #
# Refractive index                          #
# = = = = = = = = = = = = = = = = = = = = = #

function nκ(ω, ϵr1, ϵr2)
	return sqrt(0.5 * (sqrt(ϵr1*ϵr1 + ϵr2*ϵr2) + ϵr1)) + im * sqrt(0.5 * (sqrt(ϵr1*ϵr1 + ϵr2*ϵr2) - ϵr1)) 
end

# = = = = = = = = = = = = = = = = = = = = = #
# Transmission and Reflection coefficien2s  #
# = = = = = = = = = = = = = = = = = = = = = #

function t_12(n1::ComplexF64, n2::ComplexF64)
	return 2 * n1 / (n1 + n2)
end

function r_12(n1::ComplexF64, n2::ComplexF64)
	return (n1 - n2) / (n1 + n2)
end

function phase(ω::Float64, ns::ComplexF64, d::Float64) 
	return ω * ns / 299_792_458 * d
end

function t_123(ω::Float64, d::Float64, n1::ComplexF64, n2::ComplexF64, n3::ComplexF64)
	return t_12(n1, n2) * t_12(n2, n3) * exp(im * phase(ω, n2, d)) / ( 1 - r_12(n2, n3) * r_12(n2, n1) * exp(im * 2 * phase(ω, n2, d)))
end

function r_123(ω::Float64, d::Float64, n1::ComplexF64, n2::ComplexF64, n3::ComplexF64)
	return r_12(n1, n2) + t_12(n1, n2) * r_12(n2, n3) * t_12(n2, n1) * exp(im * 2 * phase(ω, n2, d)) / ( 1 - r_12(n2, n3) * r_12(n2, n1) * exp(im * 2 * phase(ω, n2, d)))
end

function t_1234(ω::Float64, d2::Float64, d3::Float64, n1::ComplexF64, n2::ComplexF64, n3::ComplexF64, n4::ComplexF64)
	return t_123(ω, d2, n1, n2, n3) * t_12(n3, n4) * exp(im * phase(ω, n3, d3)) / ( 1 - r_12(n3, n4) * r_123(ω, d2, n3, n2, n1) * exp(im * 2 * phase(ω, n3, d3)))
end

function r_1234(ω::Float64, d2::Float64, d3::Float64, n1::ComplexF64, n2::ComplexF64, n3::ComplexF64)
	return r_123(ω, d1, n1, n2, n3) + t_123(n3, n2) * r_12(n2, n1) * t_12(n2, n3) * exp(im * 2 * phase(ω, n2, d)) / ( 1 - r_12(n3, n4) * r_123(ω, d2, n3, n2, n1) * exp(im * 2 * phase(ω, n3, d3)))
end

function T(ω::Float64, d2::Float64, n1::ComplexF64, n2::ComplexF64, n3::ComplexF64)
	return abs2(t_123(ω, d2, n1, n2, n3))
end
 
function T(ω::Float64, d2::Float64, d3::Float64, n1::ComplexF64, n2::ComplexF64, n3::ComplexF64, n4::ComplexF64)
	return abs2(t_1234(ω, d2, d3, n1, n2, n3, n4))
end

function R(ω::Float64, d2::Float64, n1::ComplexF64, n2::ComplexF64, n3::ComplexF64)
	return abs2(r_123(ω, d2, n1, n2, n3))
end

function R(ω::Float64, d2::Float64, d3::Float64, n1::ComplexF64, n2::ComplexF64, n3::ComplexF64, n4::ComplexF64)
	return abs2(r_1234(ω, d2, d3, n1, n2, n3, n4))
end

function Tmodel(ω::Float64, d1::Float64, params1::NTuple{N}) where N
	# air_refractive = 1.000_293 + 0 * im
	return T(
		ω, d1, 
		1.000_293 + 0 * im,
		nκ(ω, ϵr12(ω, ϵrParams(params1...)...)...),
		1.000_293 + 0 * im
	)
end

function Tmodel(ω::Float64, d1::Float64, d2::Float64, params1::NTuple{N},  params2::NTuple{N}) where N
	# air_refractive = 1.000_293 + 0 * im
	return T(
		ω, d1, d2,
		1.000_293 + 0 * im,
		nκ(ω, ϵr12(ω, ϵrParams(params1...)...)...),
		nκ(ω, ϵr12(ω, ϵrParams(params2...)...)...),
		1.000_293 + 0 * im
	)
end

function ω(λ::Vector{Float64})
	# λ un1t be nm
	# 299_792_458 * 1.000_293 / 1e-9 = 2.99880297190194e17
    return [2 * pi * 2.99880297190194e17 / i for i in λ]
end

#=
LestSquare function:
--------------
	Calculate one or two unknown reflection indices from a transmittance spectrum.

Parameters:
--------------
	1. ωs      := Spectrum frequency.
	2. Texp    := Transmittance of examperiment
	3. d1      := First unknow material thickness (unit: m)
	4. d2      := Second unknow material thickness (unit: m)
	5. params1 := First unknow material parameters (ϵ_∞, ωₚ, fⱼ, Ωⱼ, γj)
	6. params2 := Second unknow material parameters (ϵ_∞, ωₚ, fⱼ, Ωⱼ, γj)

Examples:
--------------
	const DAT
	501×2 Matrix{Float64}:
	850.0  110.464
	849.0   75.2251
		⋮ 
	350.0   67.2263

	params_ = (0.9, 4.1e15, 4.5, 8.8e15, 1.e11)
	LestSquare(ω(DAT[:,1]), DAT[:,2]*0.01, 0.55e-3, params_)
	>> 16.28426066492815
=#
function LestSquare(ωs::Vector, Texp::Vector, d1::Float64, params1::NTuple{N}) where N
	sum = 0.0
	@inbounds for i in eachindex(ωs)
		sum += abs2(Tmodel(ωs[i], d1, params1) - Texp[i])
	end
	return sum
end

function LestSquare(ωs::Vector, Texp::Vector, d1::Float64, d2::Float64, params1::NTuple{N}, params2::NTuple{N}) where N
	sum = 0.0
	@inbounds for i in eachindex(ωs)
		sum += abs2(Tmodel(ωs[i], d1, d2, params1, params2) - Texp[i])
	end
	return sum
end

end # module Dielectics