"""
Construct a basis for identity tangent space of differomorphisms, using the fourier
basis. The notation for the is in accordance with the project thesis, with "_"
 denoting the "flipping"  [bi(x, y), 0] -> [0, bi(y, x)]
"""

"""
TYPE 1 Basis Function: Sines
"""
function ξ(k)
    return x ->  [√2*sin(π*k*x[1]), 0.]
end

function ξx(k)
    return x -> √2π*k*cos(π*k*x[1]) 
end

function Dξ(k)
    return x -> [[√2π*k*cos(π*k*x[1]) 0.]; [0. 0.]]
end



"""
TYPE 2 Basis Function: Sines*Cosines
"""
function η(k, l)
    return x -> [2*sin(π*k*x[1])*cos(2π*l*x[2]), 0.]
end

function ηx(k, l)
    return x -> 2π*k*cos(π*k*x[1])*cos(2π*l*x[2])
end

function Dη(k, l)
    return x -> [[2π*k*cos(π*k*x[1])*cos(2π*l*x[2]) -4π*l*sin(π*k*x[1])*sin(2π*l*x[2])]; [0. 0.]]
end

"""
TYPE 3 Basis Function: Sines*Sines
"""
function φ(k, l)
    return x -> [2*sin(π*k*x[1])*sin(2π*l*x[2]), 0.]
end

function φx(k, l)
    return x -> 2π*k*cos(π*k*x[1])*sin(2π*l*x[2])
end

function Dφ(k, l)
    return x -> [[2π*k*cos(π*k*x[1])*sin(2π*l*x[2]) 4π*l*sin(π*k*x[1])*cos(2π*l*x[2])];[0. 0.]]
end


"""
***********************************
Now for the y-driection vector fields
*************************************
"""

"""
TYPE 1 Basis Function: Sines
"""
function _ξ(l)
    return x ->  [0., √2*sin(π*l*x[2])]
end

function _ξy(l)
    return x -> √2π*l*cos(π*l*x[2]) 
end

function _Dξ(l)
    return x -> [[0. 0.]; [0. √2π*l*cos(π*l*x[2])]]
end


"""
TYPE 2 Basis Function: Sines*Cosines
"""
function _η(l, k)
    return x -> [0., 2*sin(π*l*x[2])*cos(2π*k*x[1])]
end

function _ηy(l, k)
    return x -> 2π*l*cos(π*l*x[2])*cos(2π*k*x[1])
end

function _Dη(l, k)
    return x -> [[0. 0.]; [-4π*k*sin(π*l*x[2])*sin(2π*k*x[1]) 2π*l*cos(π*l*x[2])*cos(2π*k*x[1])]]
end


"""
TYPE 2 Basis Function: Sines*Sines
"""
function _φ(l, k)
    return x -> [0., 2*sin(π*l*x[2])*sin(2π*k*x[1])]
end

function _φy(l, k)
    return x -> 2π*l*cos(π*l*x[2])*sin(2π*k*x[1])
end

function _Dφ(l, k)
    return x -> [[0. 0.]; [4π*k*sin(π*l*x[2])*cos(2π*k*x[1]) 2π*l*cos(π*l*x[2])*sin(2π*k*x[1])]]
end


"""
Function for constructing bass functions with necessary derivatives
    * K::Int : Maximal Frequency of trigonometric functions

Returns:
    * B, DB, divB, H1, H2::Array{Function}
"""
function construct_basis(K::Int)
    N = K + 2K^2
    
    # Declare output function arrays
    B = Array{Function}(undef, 2N) # Basis functions
    DB = Array{Function}(undef, 2N) # Jacobians of basis functions
    divB = Array{Function}(undef, 2N) # Divergence of basis functions.
    
    # Loop through type 1 basis functions (only sines)
    for i in 1:K
        k = i
        
        # Basis Functions
        B[i] = ξ(k)
        B[N + i] = _ξ(k)
        
        # Jacobian Matrices
        DB[i] = Dξ(k)
        DB[N + i] = _Dξ(k)

        # Divergence 
        divB[i] = ξx(k)
        divB[N + i] = _ξy(k)
    end
    

    for i in (K+1):(K+K^2)
        # Remap indices i -> j, k, l
        j = i - K
        k = ((j-1) ÷ K) + 1
        l = j - K * ((j-1) ÷ K)
        
        # Basis Functions
        B[i] = η(k, l)
        B[N + i] = _η(k, l)
        
        # Jacobian Matrices
        DB[i] = Dη(k, l)
        DB[N + i] = _Dη(k, l)

        # Divergence 
        divB[i] = ηx(k, l)
        divB[N + i] = _ηy(k, l)
    end
    
    for i in (K + K^2 + 1):N
        # Remap indices i -> j, k, l
        j = i - (K + K^2)
        k = ((j-1) ÷ K) + 1
        l = j - K * ((j-1) ÷ K)
        

        # Basis Functions
        B[i] = φ(k, l)
        B[N + i] = _φ(k, l)
        
        # Jacobian Matrices
        DB[i] = Dφ(k, l)
        DB[N + i] = _Dφ(k, l)

        # Divergence 
        divB[i] = φx(k, l)
        divB[N + i] = _φy(k, l)
    end
    
    return B, DB, divB
end