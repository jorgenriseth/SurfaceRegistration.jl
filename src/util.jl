using ForwardDiff


function check_derivatives(g, gx, gy)
    p = rand(2)
    display(ForwardDiff.jacobian(g, p))
    display([gx(p) gy(p)])
    return ForwardDiff.jacobian(g, p) ≈ [gx(p) gy(p)]
end
    

function check_second_derivatives(f, fx, fy, fxx, fxy, fyy)
    p = rand(2)
    display([fxx(p) fxy(p) fyy(p)])
    display(ForwardDiff.jacobian(fx, p))
    display(ForwardDiff.jacobian(fy, p))
end