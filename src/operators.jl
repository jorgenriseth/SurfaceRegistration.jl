using Cubature
using LinearAlgebra

# Syntactic Sugar
⊗(v::AbstractArray, u::AbstractArray) = kron(v, u')

# Integrates ||F1(x,y) - F2(x, y)|| over unit square [0,1]² 
function l2_distance(F1::Function, F2::Function; tol=1e-6)
    diff = x::Vector -> norm(F1(x[1], x[2]) - F2(x[1], x[2]))^2
    return √(pcubature(diff, [0, 0], [1, 1]; reltol=tol, abstol=tol)[1])
end

# Calculate the L2-innerproduct between two function F, G: R^2 -> R^3
function inner_product(F::Function, G::Function; tol=1e-6)
    H = x::Vector -> F(x[1], x[2]) ⋅ G(x[1], x[2])
    return pcubature(H, [0, 0], [1, 1]; reltol=tol, abstol=tol)[1]
end

# Integrate a function f(x, y) -> R^1 on the domain [xmin, xmax] x [ymin, ymax].
# Assumes rectangular domain.
function tangent_vector(q::Function, r::Function)
    factor = 1. / l2_distance(q, r)
    return (x::Real, y::Real) -> factor * (q(x, y) - r(x, y))
end

function l2_squared_norm(f, tol=1e-5)
    return pcubature(x -> norm(f(x))^2, [0, 0], [1, 1], rtol=tol, atol=tol)[1]
end

function l2_squared_distance(f, g, tol=1e-5)
    return pcubature(x -> norm(f(x) - g(x))^2, [0, 0], [1, 1], rtol=tol, atol=tol)[1]
end

# Calculate the L2-innerproduct between two function f, g: R^2 -> R^3
function l2_inner_product(f, g, tol=1e-5)
    return pcubature(x -> f(x) ⋅ g(x), [0, 0], [1, 1], rtol=tol, atol=tol)[1]
end