# Tomasz Niedziałek 279754

function mstycznych(f,pf,x0::Float64, delta::Float64, epsilon::Float64, maxit::Int)::Tuple{Float64,Float64,Int,Int}
    v = f(x0)
    if abs(v) < epsilon
        return (x0, v, 0, 0)
    end
    
    if pf(x0) == zero(Float64)
            return (0, 0, 0, 2)
        end
    x1 = x0 - v/pf(x0)

    for k in 1:maxit
        if pf(x0) == zero(Float64)
            return (0, 0, 0, 2)
        end
        x1 = x0 - v/pf(x0)
        v = f(x1)
        if abs(x1 - x0) < delta || abs(v) < epsilon
            return (x1, v, k, 0)
        end
        x0 = x1
    end
    return (x1, v, maxit, 1)
end

