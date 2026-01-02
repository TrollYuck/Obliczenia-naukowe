#Tomasz Niedziałek 279754

using Printf

function f(x::Float64) # funkcja z zadania
    return sin(x) + cos(3*x)
end

function df(x::Float64) # pochodna rzeczywista funkcji
    return cos(x) - 3*sin(3*x)
end

function approx_df(x::Float64, h::Float64) # pochodna przybliżona wzorem z zadania
    return (f(x+h) - f(x))/h
end

function test(x::Float64; maxn::Integer = 54, outpath::AbstractString = "zad7_results.csv")
    exact = df(x)
    open(outpath, "w") do io
        println(io, "n,h,approx,exact,error")
        for n in 0:maxn
            h = 2.0^-n
            approx = approx_df(x, h)
            err = abs(exact - approx)
            @printf(io, "%d,%.17g,%.17g,%.17g,%.17g\n", n, h, approx, exact, err)
        end
    end
    println("Wrote results to ", outpath)
end

test(1.0)