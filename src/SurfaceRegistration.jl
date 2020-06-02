module SurfaceRegistration

using LinearAlgebra: norm, tr, det, ⋅, ×
using HCubature: hcubature
using Cubature: pcubature
using PyPlot
using ForwardDiff: derivative, gradient, jacobian

include("fourier-basis.jl")
export construct_basis

include("operators.jl")
export l2_squared_norm, l2_squared_distance, l2_inner_product, tangent_vector, unscaled_tangent_vector

include("vectorfields.jl")
export construct_vectorfield, construct_jacobian_matrix
export jacobian_determinant, gradient_jacobian_determinant

include("visual.jl")
export plot_parametric_surface, plot_parametric_wireframe, plot_diffeomorphism

include("qmap.jl")
export Qmap, af, Daf, normalize_surface

include("gradient_descent.jl")
export SurfaceRegistrationAlgorithm

end # module
