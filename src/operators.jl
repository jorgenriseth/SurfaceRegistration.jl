using Cubature
using LinearAlgebra

function l2_squared_norm(f; tol=1e-12)
    return pcubature(x -> norm(f(x))^2, [0., 0.], [1., 1.], abstol=tol)[1]
end

function l2_squared_distance(f, g; tol=1e-12)
    return pcubature(x -> norm(f(x) - g(x))^2, [0., 0.], [1., 1.], abstol=tol)[1] #, rtol=tol, atol=tol)[1]
end

# Calculate the L2-innerproduct between two function f, g: R^2 -> R^3
function l2_inner_product(f, g; tol=1e-12)
    function euclidean_inner_product(x)
        return f(x) ⋅ g(x)
    end
    return pcubature(euclidean_inner_product, [0., 0.], [1., 1.], abstol=tol)[1] #, rtol=tol, atol=tol)[1]
end

function tangent_vector(q, r)
    scaling = √l2_squared_distance(q, r)
    return x -> (q(x) - r(x) ) / scaling
end

function unscaled_tangent_vector(q, r)
    return x -> q(x) - r(x)
end