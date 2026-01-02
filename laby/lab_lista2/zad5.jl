#Tomasz Niedziałek 279754
function logmodel64(n::Int, p_0::Float64, r::Float64, i::Int)
    println("n = ", i, ", p_n = ", p_0)
    if n <= i
        return p_0  
    else
        p_n = p_0 + r*p_0*(1-p_0)
        return logmodel64(n, p_n, r, i + 1)
    end
end

function logmodel32(n::Int, p_0::Float32, r::Float32, i::Int)
    println("n = ", i, ", p_n = ", p_0)
    if n <= i
        return p_0  
    else
        p_n = p_0 + r*p_0*(1-p_0)
        return logmodel32(n, p_n, r, i + 1)
    end
end



function logmodelwithstop(n::Int, p_0::Float32, r::Float32, i::Int)
    println("n = ", i, ", p_n = ", p_0)
    if n <= i
        return p_0  
    elseif i == 10
        p_0 = p_0*10^3
        p_0 = trunc(p_0)
        p_0 = p_0 / 10^3
        p_n = p_0 + r*p_0*(1-p_0)
        return logmodelwithstop(n, p_n, r, i + 1)
    else
        p_n = p_0 + r*p_0*(1-p_0)
        return logmodelwithstop(n, p_n, r, i + 1)
    end
end

println("---logmodel32---")
logmodel32(40, Float32(0.01), Float32(3.0), 0)

println("---logmodelwithstop---")
logmodelwithstop(40, Float32(0.01), Float32(3.0), 0)

println("---logmodel64---")
logmodel64(40, 0.01, 3.0, 0)
