#Tomasz Niedziałek 279754

function forward(x, y) # algorytm a dodawanie "w przód"
    @assert length(x) == length(y)
    S = zero(eltype(x)) #zero dla arytmetiki w której jest x
    for i in eachindex(x,y)
        S += x[i] * y[i]
    end
    return S
end

function backward(x, y) # algorytm b dodawanie "od tyłu"
    @assert length(x) == length(y)
    S = zero(eltype(x))
    for i in reverse(eachindex(x, y))
        S += x[i] * y[i]
    end
    return S
end

function large_1st(x,y) # algorytm c dodawanie "malejąco"
    @assert length(x) == length(y)
    T = eltype(x)
    p = T[x[i] * y[i] for i in eachindex(x,y)]
    pp = filter(z -> z >= 0, p)
    pn = filter(z -> z <0, p)
    sort!(pp, rev=true)
    sort!(pn, rev=true)
    sum_pp = zero(T)
    sum_pn = zero(T)
    for i in eachindex(pp)
        sum_pp += pp[i]
    end
    for i in eachindex(pn)
        sum_pn += pn[i]
    end
    return sum_pp + sum_pn
end

function small_1st(x,y) #algorytm d dodawanie "rosnąco"
    @assert length(x) == length(y)
    T = eltype(x)
    p = T[x[i] * y[i] for i in eachindex(x,y)]
    pp = filter(z -> z >= 0, p)
    pn = filter(z -> z <0, p)
    sort!(pp, rev=false)
    sort!(pn, rev=false)
    sum_pp = zero(T)
    sum_pn = zero(T)
    for i in eachindex(pp)
        sum_pp += pp[i]
    end
    for i in eachindex(pn)
        sum_pn += pn[i]
    end
    return sum_pn + sum_pp  
end


x64 = [2.718281828, -3.141592654, 1.414213562, 0.5772156649, 0.3010299957]
y64 = [1486.2497, 878366.9879, -22.37492, 4773714.647, 0.000185049]

x32 = Float32.(x64)
y32 = Float32.(y64)

println("a32 : ", forward(x32, y32))
println("a64 : ", forward(x64, y64))
println("b32 : ", backward(x32, y32))
println("b64 : ", backward(x64, y64))
println("c32 : ", large_1st(x32, y32))
println("c64 : ", large_1st(x64, y64))
println("d32 : ", small_1st(x32, y32))
println("d64 : ", small_1st(x64, y64))