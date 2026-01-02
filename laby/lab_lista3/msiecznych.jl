# Tomasz Niedziałek 279754

function msiecznych(f, x0::Float64, x1::Float64, delta::Float64, epsilon::Float64,maxit::Int)::Tuple{Float64,Float64,Int,Int}
    a = x0
    b = x1
    fa = f(a)
    fb = f(b)
    for k in 1:maxit
        if abs(fa) > abs(fb)
            a,b = b,a
            fa,fb = fb,fa
        end
        s = (b-a)/(fb-fa)
        b = a
        fb = fa
        a = a - fa*s 
        fa = f(a)
        if abs(b-a) < delta || abs(fa) < epsilon
            return (a, fa, k, 0)
        end
    end
    return (a, fa, maxit, 1)
end
