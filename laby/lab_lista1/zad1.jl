#Tomasz Niedziałek 279754

function iter_eps(FloatX) # epsilon znajdowany iteracyjnie
    eps = one(FloatX) # jeden w danym typie
    prev = eps
    while one(FloatX) + eps != one(FloatX)
        prev = eps
        eps /= 2 * one(FloatX)
    end
    return prev
end

function iter_eta(FloatX) # eta znajdowana iteracyjnie
    eta = one(FloatX)
    prev = eta
    while eta != 0
        prev = eta
        eta /= 2 * one(FloatX)
    end
    return prev
end

function iter_max(FloatX) # tutaj mnożenie jak poprzednio, 64 nie chciał się w ten sposób policzyć (co nie jest dziwne ale to był mój pierwszy pomysł)
    max = 2 * one(FloatX)
    prev = max
    while !isinf(max)
        prev = max
        max *= 2 * one(FloatX)
    end
    max = prev
    while !isinf(max)
        prev = max
        max *= 1 + eps(FloatX)
    end
    return prev
end

function iter_max_faster(FloatX) # tutaj wypełniamy mantysę coraz mniejszymi wartościami | iteracyjne szukanie max

    max = 2 * one(FloatX)
    while !isinf(max * 2)
        max *= 2
    end

    increment = max / 2

    while increment > 0
        if !isinf(max + increment)
            max += increment
        end
        increment /= 2
    end
    
    return max
end

println("eps Float16 iter: ", iter_eps(Float16))
println("eps Float16 Julia: ", eps(Float16), "\n")

println("eps Float32 iter: ", iter_eps(Float32))
println("eps Float32 Julia: ", eps(Float32), "\n")

println("eps Float64 iter: ", iter_eps(Float64))
println("eps Float64 Julia: ", eps(Float64), "\n")

println("eta Float16 iter: ", iter_eta(Float16))
println("eta Float16 Julia: ", nextfloat(Float16(0.0)), "\n")

println("eta Float32 iter: ", iter_eta(Float32))
println("eta Float32 Julia: ", nextfloat(Float32(0.0)), "\n")

println("eta Float64 iter: ", iter_eta(Float64))
println("eta Float64 Julia: ", nextfloat(Float64(0.0)), "\n")


println("floatmin(Float32) = ", floatmin(Float32), "\n")
println("floatmin(Float64) = ", floatmin(Float64), "\n")


println("MAX Float16 iter: ", iter_max_faster(Float16))
println("MAX Float16 Julia: ", floatmax(Float16), "\n")

println("MAX Float32 iter: ", iter_max_faster(Float32))
println("MAX Float32 Julia: ", floatmax(Float32), "\n")

println("MAX Float64 iter: ", iter_max_faster(Float64))
println("MAX Float64 Julia: ", floatmax(Float64), "\n")