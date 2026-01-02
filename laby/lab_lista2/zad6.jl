#Tomasz Niedziałek 279754

using Printf, Plots

function eqation(x0::Float64, c::Float64, n::Int, i::Int)
    if n <= i
        return x0
    else
        xn = x0^2 + c 
        return eqation(xn, c, n, i + 1)
    end
end


test_cases = [
    (c = -2.0, x0 = 1.0),
    (c = -2.0, x0 = 2.0),
    (c = -2.0, x0 = 1.99999999999999), 
    (c = -1.0, x0 = 1.0),
    (c = -1.0, x0 = -1.0),
    (c = -1.0, x0 = 0.75),
    (c = -1.0, x0 = 0.25)
]

num_iterations = 40 

function run_eqation_sequence(x0::Float64, c::Float64, num_iterations::Int)
    [eqation(x0, c, k, 0) for k in 0:num_iterations]
end

function run_all(test_cases, num_iterations::Int)
    results = Vector{Vector{Float64}}(undef, length(test_cases))
    for (idx, tc) in enumerate(test_cases)
        results[idx] = run_eqation_sequence(tc.x0, tc.c, num_iterations)
    end
    return results
end

results = run_all(test_cases, num_iterations)

function print_grouped_tables(test_cases, results, num_iterations::Int)
    groups = Dict{Float64, Vector{Int}}()
    for (i, tc) in enumerate(test_cases)
        push!(get!(groups, tc.c, Int[]), i)
    end

    for c in sort(collect(keys(groups)))
        inds = groups[c]
        println("c = $(c):")
        @printf("%4s", "i")
        for idx in inds
            @printf("  %25s", @sprintf("x0=%.17f", test_cases[idx].x0))
        end
        println()
        for i in 0:num_iterations
            @printf("%4d", i)
            for idx in inds
                val = results[idx][i+1]
                @printf("  %25.17f", val) 
            end
            println()
        end
        println()
    end
end

print_grouped_tables(test_cases, results, num_iterations)




# Grafy
function plot_cobweb_diagram(c::Float64, x0::Float64, sequence::Vector{Float64})
    finite_sequence = filter(isfinite, sequence)
    
    if isempty(finite_sequence)
        min_val = -3.0
        max_val = 3.0
    else
        min_val = minimum(finite_sequence)
        max_val = maximum(finite_sequence)
    end

    if abs(max_val - min_val) < 1e-6 
        range_center = (min_val + max_val) / 2
        plot_range_min = range_center - 1.0
        plot_range_max = range_center + 1.0
    else
        padding = (max_val - min_val) * 0.2
        plot_range_min = min_val - padding
        plot_range_max = max_val + padding
    end
    
    plot_range_min = min(plot_range_min, -2.5, (1 - sqrt(complex(1-4c)))/2 |> real, (1 + sqrt(complex(1-4c)))/2 |> real)
    plot_range_max = max(plot_range_max, 2.5, (1 - sqrt(complex(1-4c)))/2 |> real, (1 + sqrt(complex(1-4c)))/2 |> real)

    x_plot_range = range(plot_range_min, stop=plot_range_max, length=500)
    y_func = x_plot_range.^2 .+ c
    
    p = plot(x_plot_range, y_func, label="f(x) = x^2 + c", color=:red, linestyle=:dash, linewidth=2,
             xlabel="x", ylabel="f(x)", 
             title="Diagram for x_{n+1} = x_n^2 + c\n(c=$c, x0=$x0)",
             legend=:topleft, grid=true, aspect_ratio=:equal, size=(700, 700))

    plot!(x_plot_range, x_plot_range, label="y = x", color=:green, linestyle=:dot, linewidth=2)

    current_x = sequence[1]
    for k in 1:(length(sequence)-1)
        next_x = sequence[k+1]
        
        if isnan(current_x) || isnan(next_x)
            break
        end

        plot!([current_x, current_x], [current_x, next_x], color=:purple, linewidth=0.8, alpha=0.7, label="", seriestype=:path)
        
        plot!([current_x, next_x], [next_x, next_x], color=:purple, linewidth=0.8, alpha=0.7, label="", seriestype=:path)
        
        current_x = next_x
    end
    
    plot!(finite_sequence[1:end-1], finite_sequence[1:end-1], seriestype=:scatter, marker=:circle, markersize=3, color=:blue, label="x_n na y=x")

    display(p)
    
    filename_base = replace("cobweb_c_$(c)_x0_$(x0)", "." => "p", "-" => "m")
    savefig(p, "$filename_base.png")
end

println("--- Generowanie grafów ---")

all_results = Vector{Vector{Float64}}(undef, length(test_cases))

for (idx, tc) in enumerate(test_cases)
    println("\nPrzetwarzanie przypadku $(idx): c = $(tc.c), x0 = $(tc.x0)")
    sequence = run_eqation_sequence(tc.x0, tc.c, num_iterations)
    all_results[idx] = sequence 

    plot_cobweb_diagram(tc.c, tc.x0, sequence)
end

println("\n--- Zakończono generowanie wykresów ---")