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

x64_2 = [2.718281828, -3.141592654, 1.414213562, 0.577215664, 0.301029995]
y64_2 = [1486.2497, 878366.9879, -22.37492, 4773714.647, 0.000185049]

x32 = Float32.(x64)
y32 = Float32.(y64)

x32_2 = Float32.(x64_2)
y32_2 = Float32.(y64_2)

f32_1 = forward(x32, y32)
f64_1 = forward(x64, y64)
b32_1 = backward(x32, y32)
b64_1 = backward(x64, y64)
l32_1 = large_1st(x32, y32)
l64_1 = large_1st(x64, y64)
s32_1 = small_1st(x32, y32)
s64_1 = small_1st(x64, y64)

f32_2 = forward(x32_2, y32_2)
f64_2 = forward(x64_2, y64_2)
b32_2 = backward(x32_2, y32_2)
b64_2 = backward(x64_2, y64_2)
l32_2 = large_1st(x32_2, y32_2)
l64_2 = large_1st(x64_2, y64_2)
s32_2 = small_1st(x32_2, y32_2)
s64_2 = small_1st(x64_2, y64_2)

println("a32 : ", f32_1)
println("a64 : ", f64_1)
println("b32 : ", b32_1)
println("b64 : ", b64_1)
println("c32 : ", l32_1)
println("c64 : ", l64_1)
println("d32 : ", s32_1)
println("d64 : ", s64_1)

println("a32_2 : ", f32_2)
println("a64_2 : ", f64_2)
println("b32_2 : ", b32_2)
println("b64_2 : ", b64_2)
println("c32_2 : ", l32_2)
println("c64_2 : ", l64_2)
println("d32_2 : ", s32_2)
println("d64_2 : ", s64_2)

println("Diff")
println("a32_1/a32_2 : ", f32_1 / f32_2)
println("a64_1/a64_2 : ", f64_1 / f64_2)
println("b32_1/b32_2 : ", b32_1 / b32_2)
println("b64_1/b64_2 : ", b64_1 / b64_2)
println("c32_1/c32_2 : ", l32_1 / l32_2)
println("c64_1/c64_2 : ", l64_1 / l64_2)
println("d32_1/d32_2 : ", s32_1 / s32_2)
println("d64_1/d64_2 : ", s64_1 / s64_2)

println("Abs diffs")
println("a32_1 - a32_2 : ", abs(f32_1 - f32_2))
println("a64_1 - a64_2 : ", abs(f64_1 - f64_2))
println("b32_1 - b32_2 : ", abs(b32_1 - b32_2))
println("b64_1 - b64_2 : ", abs(b64_1 - b64_2))
println("c32_1 - c32_2 : ", abs(l32_1 - l32_2))
println("c64_1 - c64_2 : ", abs(l64_1 - l64_2))
println("d32_1 - d32_2 : ", abs(s32_1 - s32_2))
println("d64_1 - d64_2 : ", abs(s64_1 - s64_2))