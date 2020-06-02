using LinearAlgebra

include("operators.jl")
include("vectorfields.jl")

function find_point_max_ε(x, M)
    a = det(M(x))
    
    # No zeros.
    if a >= 0.
        return Inf
    end
    
    # Find slope at 0
    b = tr(M(x))
    
    # Approximately linear
    if abs(a) < 1e-5
        # Positive zero
        if b < 0.
            return 1. / abs(b)
        end 
        
        # Only negative zero
        return Inf
    end
    
    # Positive Zero of quadratic
    return (-b - √(b^2 - 4a)) / (2*a)
end


function find_max_ε(M)
    xs = 0:0.01:1
    A = [find_point_max_ε([xi, yi], M) for xi in xs, yi in xs]
    return minimum([find_point_max_ε([xi, yi], M) for xi in xs, yi in xs]) * 0.9
end


function armijo(q, r, v, Dv; c1=0.2, ρ=0.5, max_iter=50, verbose=false)
    # Set starting value for varepsilon
    ε = find_max_ε(Dv)
    
    # Define functions
    γ(x) =  x + ε * v(x)
    Dγ(x) = [[1. 0.];[0. 1.]] + ε*Dv(x)
    Jγ = jacobian_determinant(Dγ)
    ϕr(x) = √Jγ(x) * r(γ(x))
    
    # Initial Values
    Eγ0 = l2_squared_distance(q, r)
    dEγ0 = l2_squared_norm(v)
    
    # Ensure max step size allowed
    run = true
    while run
        try 
            l2_squared_distance(q, ϕr)
            run = false
        catch
            ε *= 0.9
        end
    end
    
    #Set initial values for linesearch
    lhs0 = l2_squared_distance(q, ϕr)
    rhs = Eγ0 - ε * c1 *dEγ0
                
    
    # Init 
    it = 0
    
    # Print
    if verbose
        println("-"^60)
        println("[Armijo] Iter $it:\t $lhs0 \t $rhs")
    end
    
    # Loop until Armijo Condition is met
    while lhs0 > rhs && it < max_iter
        # Update step size
        ε *= ρ
        
        # Update Search values
        lhs1 = l2_squared_distance(q, ϕr)
        rhs = Eγ0 - ε * c1 *dEγ0
        
        # Increment and print results
        it += 1
        if verbose
            println("[Armijo] Iter $it: \t $lhs1 \t $rhs")
        end
        
        # Stop if value starts increasing.
        if lhs1 > lhs0
            ε /= ρ
            return ε, γ, Dγ, Jγ, ϕr
        end
        lhs0 = lhs1
            
    end
    if verbose 
        println("[Armijo] Step Size: \t $ε")
        println("-"^60)
    end

    return ε, γ, Dγ, Jγ, ϕr
end


function update_D_q_map(J::Function, ∇J::Function, γ::Function, Dγ::Function, q::Function, Dq::Function)
    return x -> 1. / (2.0 * √J(x)) * q(γ(x)) * ∇J(γ(x))' + √J(x) * Dq(γ(x)) * Dγ(x)
end



function SurfaceRegistrationAlgorithm(
    q1, q2, Dq2, B, DB, divB; max_iter=10, error_change_tol=1e-2,
    inner_product_tol=1e-12, norm_tol=1e-12, armijo_constant=0.2,
    armijo_scaling=0.5, armijo_max_iter=50, verbose=false
)

r = Array{Function}(undef, max_iter+1)
Dr = Array{Function}(undef, max_iter+1)
v = Array{Function}(undef, max_iter)
γ = Array{Function}(undef, max_iter)
Dγ = Array{Function}(undef, max_iter)
J = Array{Function}(undef, max_iter)
∇J = Array{Function}(undef, max_iter)
r[1] = q2
Dr[1] = Dq2

# Array to hold 
errors = zeros(max_iter+1)
errors[1] = l2_squared_distance(q1, r[1], tol=norm_tol)

N = length(B)
dϕ = Array{Function}(undef, N)
coeff = zeros(N)

# Start iterating
n = 1
while n < max_iter + 1
    
    println("="^60)
    print("[SurfaceRegistration] Iter ", n-1, "\tError=")
    println(errors[n])
    println("="^60)
    
    # Tangent Vector
    v = tangent_vector(q1, r[n])
#         v = x -> q1(x) - r[n](x)
    
    # Get Vector field coefficients
    for i in 1:N
        dϕ[i] = x -> 0.5 * divB[i](x) * r[n](x) + Dr[n](x) * B[i](x)
        coeff[i] = l2_inner_product(v, dϕ[i])
    end
    
    # Construct tangent space vector fields
    dγ = construct_vectorfield(coeff, B)
    Ddγ = construct_jacobian_matrix(coeff, DB)
    
    # Armijo search to find new parameters
    ε, γ[n], Dγ[n], J[n], r[n+1] = armijo(q1, r[n], dγ, Ddγ; c1=armijo_constant, ρ=armijo_scaling, max_iter=armijo_max_iter, verbose=verbose)
    
    ∇J[n] = gradient_jacobian_determinant(Dγ[n])
    Dr[n+1] = update_D_q_map(J[n], ∇J[n], γ[n], Dγ[n], r[n], Dr[n])
    
    # Compute new error, and check if change is surrficiently small.
    errors[n + 1] = l2_squared_distance(q1, r[n+1], tol=norm_tol)
    
    # Terminate if error change is sufficiently small
    if (errors[n] - errors[n+1])/errors[n] < error_change_tol
        println("[DEBUG]", (errors[n] - errors[n+1])/errors[n])
        println("="^60)
        print("[SurfaceRegistration] Termination: Error Change. \tError=")
        println(errors[n])
        println("="^60)
        return r[1:n+1], γ[1:n], errors[1:n+1]
    end
    
    # Terminate if step_size is sufficiently small
    if ε < 1e-8
        println("="^60)
        print("[SurfaceRegistration] Termination: Step Size. \tError=")
        println(errors[n])
        println("="^60)
        return r[1:n+1], γ[1:n], errors[1:n+1]
    end        
    
    # Update Iteration
    n += 1
end

println("="^60)
println("[SurfaceRegistration] Done. \tError=", errors[n])
println("="^60)

return r, γ, errors, Dr
end