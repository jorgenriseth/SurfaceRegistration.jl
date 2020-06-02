"""
Area scaling factor for a function f with partial derivatives fx, fy
    * fx::Function
    * fy::Function

Returns:
    * [Function] : Area scaling factor at x.
"""
function af(fx, fy)
    return x -> norm(fx(x)×fy(x))
end


"""
Derivative of area scaling factor for a function f with derivatives 
    * fx::Function
    * fy::Function
    * fxx::Function
    * fxy::Function
    * fyy::Function

Returns:
    * [Function] : Derivative of area scaling factor at x.
"""
function Daf(af, fx, fy, fxx, fxy, fyy)
    return x -> [
        (fxx(x) × fy(x) - fxy(x) × fx(x)) ⋅ (fx(x) × fy(x)) / af(x)
        (fxy(x) × fy(x) - fyy(x) × fx(x)) ⋅ (fx(x) × fy(x)) / af(x)
    ]
end


"""
Get q-map and derivative of a given function
    * f::Function

Returns:
    * q::Function, Dq::Function
"""
function Qmap(f)
    # Find partial derivatives
    fx(x) = derivative(z -> f([z, x[2]]), x[1])
    fy(x) = derivative(z -> f([x[1], z]), x[2])
    
    # Compute area scaling factor
    a = af(fx, fy)
    
    # Return q-map and derivative 
    q(x) = √a(x) * f(x)
    Dq(x) = jacobian(q, x)
    return q, Dq
end


"""
Get q-map and derivative of a given function with derivatives
    * f::Function
    * fx::Function
    * fy::Function

Returns:
    * q::Function, Dq::Function
"""
function Qmap(f, fx, fy)
    a = af(fx, fy)
    
    q(x) = √a(x) * f(x)
    Dq(x) = jacobian(q, x)
    return q, Dq
end


"""
Get q-map and derivative of a given function with derivatives
    * f::Function
    * fx::Function
    * fy::Function
    * fxx::Function
    * fxy::Function
    * fyy::Function

Returns:
    * q::Function, Dq::Function
"""
function Qmap(f, fx, fy, fxx, fxy, fyy)
    a = af(fx, fy)
    Da = Daf(a, fx, fy, fxx, fxy, fyy)
    
    Df(x) = [fx(x) fy(x)]
    q(x) = √a(x) * f(x)
    Dq(x) = 1. / (2. * √a(x)) * f(x) * Da(x)' + √a(x) * Df(x)
    return q, Dq
end

"""
Compute the centroid of the parametric surface.
    * f::Function
    * af:: Function

Returns: 
    * [Array{Float}] : Vector value of surface centre
"""
function surface_centroid(f, af)
    return hcubature(x -> f(x)*af(x), [0, 0], [1, 1])[1] / hcubature(af, [0, 0], [1, 1])[1]
end


"""
Translate function so that centroid coincides with the origin.
    * f::Function
    * af:: Function

Returns: 
    * [Function] : Translated Surface
"""
function centered(f, af)
    c = surface_centroid(f, af)
    return x -> f(x) - c
end


"""
Compute the parametric surface area.
    * f::Function
    * af:: Function

Returns: 
    * [Float] : Surface Area
"""
function surface_area(f, af)
    return hcubature(af, [0, 0], [1, 1])[1]
end


"""
Scale surface to unit area.
    * f::Function
    * af:: Function

Returns: 
    * [Function] : Scaled Surface
"""
function scale_surface(f, af)
    s = surface_area(f, af)
    return x -> f(x) / s
end


"""
Translate and scale surface 
    * f::Function

Returns: 
    * [Function] : Scaled Surface
"""
function normalize_surface(f)
    # Find partial derivatives
    fx(x) = derivative(z -> f([z, x[2]]), x[1])
    fy(x) = derivative(z -> f([x[1], z]), x[2])
    
    # Compute area scaling factor
    a = af(fx, fy)
    
    c = surface_centroid(f, a)
    s = surface_area(f, a)
    return x -> (f(x) - c) / s
end


"""
Translate and scale surface 
    * f::Function
    * af:: Function

Returns: 
    * [Function] : Scaled Surface
"""
function normalize_surface(f, af)
    c = surface_centroid(f, af)
    s = surface_area(f, af)
    return x -> (f(x) - c) / s
end