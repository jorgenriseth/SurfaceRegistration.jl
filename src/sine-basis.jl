"""
Implement a fourier basis for boundary preserving vector fields on
[0, 1]². Based off the fourier series for the unit square [-1, 1]², but
restricted to only include terms that are 0 for x,y= 0, 1
"""

# ====
# TYPE 1 x-direction 
# ====

"""
Return pure sine basis functions preserving boundaries 
in x-direction.

Args:
    k::Int - Frequency of function
Returns:
    Function
"""
function sin_k_x(k::Int)::Function
    return (x::Number, y::Number) -> [√2 * sin(π*k*x), 0.]
end

"""
Partial derivative in x-direction of above function.
"""
function ∂x_sin_k_x(k::Int)::Function
    return (x::Number, y::Number) -> √2 * π*k * cos(π*k*x)
end


"""
Partial derivative in y-direction of above function.
"""
function ∂y_sin_k_x(k::Int)::Function
    return (x::Number, y::Number) -> 0.
end

"""
Jacobian of type 1 basis functions.
"""
function D_sin_k_x(k::Int)::Function
    return (x::Number, y::Number) -> [[√2*π*k*cos(π*k*x) 0.]; [0. 0.]] 
end

# =====
# TYPE 1 y-direction
# =====

"""
Return pure sine basis functions preserving boundaries 
in y-direction.

Args:
    k::Int - Frequency of function
Returns:
    Function
"""
function sin_k_y(k::Int)::Function
    return (x::Number, y::Number) -> [0., √2*sin(π*k*y)]
end


"""
Partial derivative in y-direction of above function.
"""
function ∂x_sin_k_y(k::Int)::Function
    return (x::Number, y::Number) -> 0.
end


"""
Partial derivative in y-direction of above function.
"""
function ∂y_sin_k_y(k::Int)::Function
    return (x::Number, y::Number) -> √2*π*k * cos(π*k*y)
end

"""
Jacobian of type 1 basis functions in x-direction.
"""
function D_sin_k_y(k::Int)::Function
    return (x::Number, y::Number) -> [[0. 0.]; [0 √2*π*k*sin(π*k*y)]]
end

# ==== 
# TYPE 2 x-direction
# ====


"""
Return sin * sin basis functions preserving boundaries 
in both x- and y-direction.

Args:
    k::Int - Frequency of sine-x function.
    l::Int - Frequency of sine-y function.
Returns:
    Function
"""
function sin_k_sin_l_x(k::Int, l::Int)::Function
    return (x::Number, y::Number) -> [2 * sin(π*k*x) * sin(2π*l*y), 0.]
end


"""
Partial derivative in x-direction of above function.
"""
function ∂x_sin_k_sin_l_x(k::Int, l::Int)::Function
    return (x::Number, y::Number) -> 2π*k * cos(π*k*x) * sin(2π*l*y)
end


"""
Partial derivative in y-direction of above function.
"""
function ∂y_sin_k_sin_l_x(k::Int, l::Int)::Function
    return (x::Number, y::Number) -> 4π*l * sin(π*k*x) * cos(2π*l*y)
end


function D_sin_k_sin_l_x(k::Int, l::Int)::Function
    return (x::Number, y::Number) -> [[2π*k*cos(π*k*x)*sin(2π*l*y) 4π*l*sin(π*k*x)*cos(2π*l*y)]; [0. 0.]] 
end

# ==== 
# TYPE 2 y-directionn
# ====


"""
Return sin * sin basis functions preserving boundaries 
in both x- and y-direction.

Args:
    k::Int - Frequency of sine-x function.
    l::Int - Frequency of sine-y function.
Returns:
    Function
"""
function sin_k_sin_l_y(k::Int, l::Int)::Function
    return (x::Number, y::Number) -> [0., 2sin(π*k*y) * sin(2π*l*x)] # y-direction
end


"""
Partial derivative in x-direction of above function.
"""
function ∂x_sin_k_sin_l_y(k::Int, l::Int)::Function
    return (x::Number, y::Number) -> 4π*l * sin(π*k*y) * cos(2π*l*x)
end


"""
Partial derivative in y-direction of above function.
"""
function ∂y_sin_k_sin_l_y(k::Int, l::Int)::Function
    return (x::Number, y::Number) -> 2π*k * cos(π*k*y) * sin(2π*l*x)
end

"""
Derivative/Jacobian of type 2 functions in y-direction
"""
function D_sin_k_sin_l_y(k::Int, l::Int)::Function
    return (x::Number, y::Number) -> [[0. 0.]; [4π*l*sin(π*k*y)*cos(2π*l*x) 2π*k*cos(π*k*y)*sin(π*l*x)]]
end

# ==== 
# TYPE 3 x-direction
# ==== 

"""
Return sin * cos basis function preserving boundaries
in x-diretion.

Args:
    k::Int - Frequency of sine function.
    l::Int - Frequency of cos function.
Returns:
    Function
"""
function sin_k_cos_l_x(k::Int, l::Int)::Function
    return (x::Number, y::Number) -> [2 * sin(π*k*x) * cos(2π*l*y), 0.]
end


"""
Partial derivative in x-direction of above function.
"""
function ∂x_sin_k_cos_l_x(k::Int, l::Int)::Function
    return (x::Number, y::Number) -> 2π*k * cos(π*k*x) * cos(2π*l*y)
end


"""
Partial derivative in y-direction of above function.
"""
function ∂y_sin_k_cos_l_x(k::Int, l::Int)::Function
    return (x::Number, y::Number) -> -4π*l * sin(π*k*x) * sin(2π*l*y)
end


"""
Derivative of Type 3 functions in x-direction
"""
function D_sin_k_cos_l_x(k::Int, l::Int)::Function
    return (x::Number, y::Number) -> [[2π*k*cos(π*k*x)*cos(2π*l*y) -4π*l*sin(π*k*x)*sin(2π*l*y)]; [0. 0.]]
end

# ==== 
# TYPE 3 y-direction
# ====

"""
Return sin * cos basis function preserving boundaries
in y-direction.

Args:
    k::Int - Frequency of sin function.
    l::Int - Frequency of cos function.
Returns:
    Function
"""
function sin_k_cos_l_y(k::Int, l::Int)::Function
    return (x::Number, y::Number) -> [0., 2 * sin(π*k*y) * cos(2π*l*x)]
end


"""
Partial derivative in x-direction of above function.
"""
function ∂x_sin_k_cos_l_y(k::Int, l::Int)::Function
    return (x::Number, y::Number) -> -4π*l* sin(π*k*y) * sin(2π*l*x)
end


"""
Partial derivative in y-direction of above function.
"""
function ∂y_sin_k_cos_l_y(k::Int, l::Int)::Function
    return (x::Number, y::Number) -> 2π*k * cos(π*k*y) * cos(2π*l*x)
end


"""
Derivative/Jacobian of type 3 functions in y-direction.
"""
function D_sin_k_cos_l_y(k::Int, l::Int)::Function
    return (x::Number, y::Number) -> [[0. 0.]; [-4π*l*sin(π*k*y)*sin(2π*l*x) 2π*k*cos(π*k*y)*cos(2π*l*y)]]
end

function build_basis(K::Int)
    N = K + 2K^2 # Number of bvasis functions per direction
    
    # Declare output function arrays.
    B = Array{Function}(undef, 2N) # Basis functions
    DB = Array{Function}(undef, 2N) # Jacobians of basis functions
    divB = Array{Function}(undef, 2N) # Divergence of basis functions.
    
    # Loop through type 1 basis functions (only sines)
    for i in 1:K
        k = i
        
        # Basis functions
        B[i] = sin_k_x(k)
        B[N + i] = sin_k_y(k)
        
        # Derivatives/Jacobian
        DB[i] = D_sin_k_x(k)
        DB[N + i] = D_sin_k_y(k)
        
        # Divergence operator
        divB[i] = ∂x_sin_k_x(k)
        divB[N + i] = ∂y_sin_k_y(k)
    end
    
    # Loop through type 2 basis functions (sin * sin)
    for i in (K+1):(K+K^2)
        # Remap indices i -> j, k, l
        j = i - K
        k = ((j-1) ÷ K) + 1
        l = j - K * ((j-1) ÷ K)
        
        # Basis fucntions
        B[i] = sin_k_sin_l_x(k, l)
        B[N + i] = sin_k_sin_l_y(k, l)
        
        # Derivatives 
        DB[i] = D_sin_k_sin_l_x(k, l)
        DB[N + i] = D_sin_k_sin_l_y(k, l)
        
        # Divergence operators
        divB[i] = ∂x_sin_k_sin_l_x(k, l)
        divB[N + i] = ∂y_sin_k_sin_l_y(k, l)
    end
    
    # Loop through type 3 basis functions (sin * cos)
    for i in (K + K^2 + 1):N
        # Remap indices i -> j, k, l
        j = i - (K + K^2)
        k = ((j-1) ÷ K) + 1
        l = j - K * ((j-1) ÷ K)
        
        # Basis functions
        B[i] = sin_k_cos_l_x(k, l)
        B[N + i] = sin_k_cos_l_y(k, l)
        
        # Derivatives
        DB[i] = D_sin_k_cos_l_x(k, l)
        DB[N + i] = D_sin_k_cos_l_y(k, l)
        
        # Divergence operators
        divB[i] = ∂x_sin_k_cos_l_x(k, l)
        divB[N + i] = ∂y_sin_k_cos_l_y(k, l)
    end 

    return B, DB, divB
end    