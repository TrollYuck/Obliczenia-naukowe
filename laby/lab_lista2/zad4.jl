#Tomasz Niedziałek 279754
using Polynomials, Printf

function calculateroots(w::Array{Float64})
    P = Polynomial(w)
    rP = roots(P)

    rp = collect(1:20)
    p = fromroots(rp)

    println("-------roots-------")
    println("{all abs()} k | z_k | P(z_k) | p(z_k) | z_k - k")
    for i in 1:20
        if (typeof(rP[i]) == Float64)
            @printf("%d & %.15f & %.15f & %.15f & %.15f\n", i, rP[i] ,abs(P(rP[i])), abs(p(rP[i])) , abs(rP[i] - i))
        else 
            println(i, " & ", rP[i]," & ", abs(P(rP[i])), " & ", abs(p(rP[i])), " & ", abs(rP[i] - i))
        end
    end
end

function calculatewilkinson(w::Array{Float64})
    P = Polynomial(w)
    rP = roots(P)

    println("-------Wilkinson-------")
    println("{all abs()} k | z_k | P(z_k) | z_k - k")
    for i in 1:20
        if (typeof(rP[i]) == Float64)
            @printf("%d & %.15f & %.15f & %.15f & %.15f\n", i, rP[i] ,abs(P(rP[i])), abs(rP[i] - i))
        else 
            println(i, " & ", rP[i]," & ", abs(P(rP[i])), " & ", " & ", abs(rP[i] - i))
        end
    end
end

function calculaterootscorrect(w::Array{Float64})

    P = Polynomial(w)

    rp = collect(1:20)
    p = fromroots(rp)

    println("-------correct-------")
    println("k | P(k) | p(k)")
    for i in 1:20
        @printf("%d & %.15f & %.15f \n", i ,abs(P(i)), p(i))
    end
end

#Wspolczynniki wielomianu z zadania 4
p=[1, -210.0, 20615.0,-1256850.0,
      53327946.0,-1672280820.0, 40171771630.0, -756111184500.0,          
      11310276995381.0, -135585182899530.0,
      1307535010540395.0,     -10142299865511450.0,
      63030812099294896.0,     -311333643161390640.0,
      1206647803780373360.0,     -3599979517947607200.0,
      8037811822645051776.0,      -12870931245150988800.0,
      13803759753640704000.0,      -8752948036761600000.0,
      2432902008176640000.0]
x = -210.0 - 2^(-23)
#Wspolczynniki do eksperymentu Wilkinsona
p1 = [1, x, 20615.0,-1256850.0,
      53327946.0,-1672280820.0, 40171771630.0, -756111184500.0,          
      11310276995381.0, -135585182899530.0,
      1307535010540395.0,     -10142299865511450.0,
      63030812099294896.0,     -311333643161390640.0,
      1206647803780373360.0,     -3599979517947607200.0,
      8037811822645051776.0,      -12870931245150988800.0,
      13803759753640704000.0,      -8752948036761600000.0,
      2432902008176640000.0]

w = p[end:-1:1]
w1 = p1[end:-1:1]

calculateroots(w)
calculaterootscorrect(w)
calculatewilkinson(w1)