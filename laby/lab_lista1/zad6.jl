#Tomasz Niedziałek 279754

function f(x::Float64)
    return sqrt(x^2 + 1) - 1
end

function g(x::Float64)
    return x^2 / (sqrt(x^2 + 1) + 1)
end

function test(x::Float64, n::Float64) # x brany do funcji w kolejnych potęgach malejąco do -n  
    for i in 1:n
        fv = f(x^(-1*i))
        gv = g(x^(-1*i))
        println("f(", x, "^(", (-1*i), ")) = ", fv)
        println("g(", x, "^(", (-1*i), ")) = ", gv)
        println("|f - g| = ", abs(fv - gv), "\n")
    end
end

test(8.0, 200.0) # wywołanie dla 8 jak w zadaniu, z arbitralnie dużą liczbą powtórzeń