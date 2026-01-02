#Tomasz Niedziałek 279754

using Plots, Printf
export ilorazyRoznicowe, warNewton, naturalna, rysujNnfx
function ilorazyRoznicowe(x::Vector{Float64}, f::Vector{Float64}) :: Vector{Float64}
    n = length(x)
    fx = zeros(n)
    for i in 1:n
        fx[i] = f[i]
    end
    for j in 2:n
        for i in n:-1:j
            fx[i] = (fx[i] - fx[i-1]) / (x[i] - x[i-j+1])
        end
    end
    return fx
end

function warNewton(x::Vector{Float64}, fx::Vector{Float64}, t::Float64) :: Float64
    n = length(x)
    nt = fx[n]
    for i in n-1:-1:1
        nt = fx[i] + (t - x[i]) * nt
    end
    return nt
end

function naturalna(x::Vector{Float64}, fx::Vector{Float64}) :: Vector{Float64}
    n = length(x)
    a = zeros(n)
    a[n] = fx[n]
    for i in n-1:-1:1
        a[i] = fx[i] - x[i] * a[i+1]
        for j in i+1:n-1
            a[j] = a[j] - x[i] * a[j+1]
        end
    end
    return a
end

function rysujNnfx(f, a::Float64, b::Float64, n::Int; nodes::Symbol = :rownoodlegle, save_to::Union{String,Nothing}=nothing)
    dpi = 100 
    x::Vector{Float64} = zeros(n+1)
    f_val::Vector{Float64} = zeros(n+1)

    if nodes == :rownoodlegle
        h::Float64 = (b - a) / n
        for i in 1:(n+1)
            x[i] = a + (i-1)*h
            f_val[i] = f(x[i])
        end
    elseif nodes == :czebyszew
        for i in 1:(n+1)
            x[i] = 0.5 * ((b - a) * cos((2i - 1) * π / (2 * (n+1))) + (b + a))
            f_val[i] = f(x[i])
        end
    else
        throw(ArgumentError("Nieobsługiwany typ węzłów: $nodes"))
    end

    c::Vector{Float64} = ilorazyRoznicowe(x, f_val)

    num_of_points = dpi * n + 1
    gap::Float64 = (b - a) / (num_of_points - 1)

    plot_x::Vector{Float64} = zeros(num_of_points)
    plot_f::Vector{Float64} = zeros(num_of_points)
    plot_poli::Vector{Float64} = zeros(num_of_points)

    plot_x[1] = a
    plot_f[1] = f_val[1]
    plot_poli[1] = f_val[1]
    for i in 2:num_of_points
        plot_x[i] = a + (i-1)*gap
        plot_f[i] = f(plot_x[i])
        plot_poli[i] = warNewton(x, c, plot_x[i])
    end

    p = plot(plot_x, [plot_f, plot_poli], 
        label = ["funkcja" "wielomian"], 
        title = "Interpolacja funkcji f dla n = $n",
        framestyle=:origin)
    return p
end

