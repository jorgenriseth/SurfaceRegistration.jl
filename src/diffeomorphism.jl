# TODO: Chack that Mutability of arrays doesn't break this. 

function build_diffeomorphism(coeff::AbstractArray, B::AbstractArray, ε::Number)
    # Return a function which may be evaluated 
    function diffeomorphism(x::Number, y::Number)
        return [x, y] + ε * sum(map(bi -> bi(x, y), B) .* coeff) # γid + εΣa*bi
    end
    return diffeomorphism
end

function build_D_diffeomorphism(coeff::AbstractVector, DB::AbstractArray, ε::Number)
    # Return a function which may be evaluated 
    function D_diffeomorphism(x::Number, y::Number)
        return [[1. 0.]; [0. 1.]] + ε * sum(map(Dbi -> Dbi(x, y), DB) .* coeff) # γid + εΣa*bi
    end
    return D_diffeomorphism
end

function chain_diffeomorphism(chain::Array{Function})
    function diffeomorphism(x::Number, y::Number)
        out = [x, y]
        for γ in chain
            out .= γ(out...)
        end
        return out
    end
    return diffeomorphism
end