#Tomasz Niedziałek 279754

function findx() # znajdowanie najmniejszej liczby x z przedziału [1,2] gdzie x * (1/x) != 1
    curr = one(Float64)
    count = zero(Float64)
    while (curr * (1/curr) == one(Float64) && curr < 2 * one(Float64))
        curr = nextfloat(curr)
        count += 1
    end
    return curr, count
end

function findy() # znajdowanie liczby y z przedziału [1,2] gdzie y * (1/y) != 1
    curr = 2 * one(Float64)
    count = zero(Float64)
    while (curr * (1/curr) == one(Float64))
        curr = prevfloat(curr)
        count += 1
    end
    return curr, count
end

x, c = findx()
y, c1 = findy()

println("x = ", x)
println("x * (1/x) = ", (x * (1/x)))
println("after ", c, " iterations.")

println("y = ", y)
println("y * (1/y) = ", (y * (1/y)))
println("after ", c1, " iterations.")