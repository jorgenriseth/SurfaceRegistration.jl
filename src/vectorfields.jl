"""
Construct vectorfields and derivatives from bases and coefficients.
"""

function construct_vectorfield(coeff, B)
    @assert length(coeff) == length(B)
    a = copy(coeff)
    function vectorfield(x)
        out = [0., 0.]
        for i in 1:length(B)
            out += a[i]*B[i](x)
        end
        return out
    end
    return vectorfield
end

function construct_jacobian_matrix(coeff, DB)
    @assert length(coeff) == length(DB)
    a = copy(coeff)
    function jacobian_matrix(x)
        out = [ [0. 0.]; [0. 0.]]
        for i in 1:length(DB)
            out += a[i]*DB[i](x)
        end
        return out
    end
    return jacobian_matrix
end

function jacobian_determinant(Df)
    function J(x)
        return det(Df(x))
    end
end

function gradient_jacobian_determinant(Df)
    J = jacobian_determinant(Df)
    function DJ(x)
        return ForwardDiff.gradient(J, x)
    end
    return DJ
end