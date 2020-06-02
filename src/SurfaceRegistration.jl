module SurfaceRegistration

include("sine-basis.jl")
include("operators.jl")
include("visual.jl")
include("diffeomorphism.jl")
include("registration.jl")
include("fourier-basis.jl")

export build_basis # From Basis
export ⊗, l2_squared_norm,  l2_squared_distance, l2_inner_product
export plot_parametric_surface # From visualize
export build_diffeomorphism, build_D_diffeomorphism, chain_diffeomorphism # From diffeomorphism
export update_q_map, update_D_q_map, Jγ, ∇Jγ
#export ξ, ξx, ξxx, ξxy, η, ηx, ηxx, ηxy, φ, φx, φxx, φxy
end # module
