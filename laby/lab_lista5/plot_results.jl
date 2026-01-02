# plot_results.jl
using CSV
using DataFrames
using Plots
using Printf

# Wczytaj dane
df = CSV.read("experiment_results.csv", DataFrame)

# Wyodrębnij rozmiar z nazwy datasetu (np. "Dane10000_1_1" -> 10000)
function extract_size(dataset_name)
    m = match(r"Dane(\d+)", string(dataset_name))
    return m !== nothing ? parse(Int, m.captures[1]) : 0
end

df.Size = extract_size.(df.Dataset)

# Posortuj po rozmiarze
sort!(df, :Size)

# Unikalne metody i rozmiary
methods = unique(df.Method)
sizes = sort(unique(df.Size))

println("Znalezione metody: ", methods)
println("Znalezione rozmiary: ", sizes)
println()

# Przygotuj dane dla każdej metody
function prepare_data(df, methods, sizes, column)
    data = Dict{String, Vector{Float64}}()
    
    for method in methods
        method_str = string(method)
        values = Float64[]
        for size in sizes
            row = df[(string.(df.Method) .== method_str) .& (df.Size .== size), :]
            if nrow(row) > 0
                val = row[1, column]
                # Obsługa "N/A" lub pustych wartości
                val_str = strip(string(val))
                if val_str == "N/A" || val_str == "" || ismissing(val)
                    push!(values, NaN)
                else
                    try
                        parsed = parse(Float64, val_str)
                        # Zamień 0.0 na małą wartość dla pamięci (żeby log działał)
                        if column == :Memory_MB && parsed == 0.0
                            push!(values, 0.01)
                        else
                            push!(values, parsed)
                        end
                    catch
                        push!(values, NaN)
                    end
                end
            else
                push!(values, NaN)
            end
        end
        data[method_str] = values
    end
    
    return data
end

# Przygotuj dane dla czasu
time_data = prepare_data(df, methods, sizes, :Time_s)

# Przygotuj dane dla pamięci
memory_data = prepare_data(df, methods, sizes, :Memory_MB)

# Przygotuj dane dla błędu względnego
error_data = prepare_data(df, methods, sizes, :RelativeError)

# Kolory i style dla każdej metody
colors = [:blue, :red, :green, :purple]
markers = [:circle, :square, :diamond, :cross]
styles = [:solid, :dash, :dot, :dashdot]

# Konwertuj metody do stringów
methods_str = string.(methods)

# === WYKRES 1: CZAS WYKONANIA ===
println("Generowanie wykresu czasu wykonania...")
p1 = plot(
    title="Czas wykonania w zależności od rozmiaru problemu",
    xlabel="Rozmiar problemu (n)",
    ylabel="Czas [s]",
    xscale=:log10,
    yscale=:log10,
    legend=:topleft,
    size=(1200, 700),
    dpi=300,
    grid=true,
    framestyle=:box,
    left_margin=15Plots.mm,
    bottom_margin=10Plots.mm,
    right_margin=5Plots.mm,
    top_margin=5Plots.mm
)

for (i, method) in enumerate(methods_str)
    times = time_data[method]
    # Usuń NaN przed plotowaniem
    valid_indices = .!isnan.(times) .& (times .> 0)
    if sum(valid_indices) > 0
        plot!(p1, sizes[valid_indices], times[valid_indices],
            label=method,
            linewidth=2.5,
            marker=markers[i],
            markersize=6,
            color=colors[i],
            linestyle=styles[i]
        )
    end
end

savefig(p1, "time_comparison.png")
println("  Zapisano: time_comparison.png")

# === WYKRES 2: ZUŻYCIE PAMIĘCI ===
println("Generowanie wykresu zużycia pamięci...")
p2 = plot(
    title="Zużycie pamięci w zależności od rozmiaru problemu",
    xlabel="Rozmiar problemu (n)",
    ylabel="Pamięć [MB]",
    xscale=:log10,
    yscale=:log10,
    legend=:topleft,
    size=(1200, 700),
    dpi=300,
    grid=true,
    framestyle=:box,
    left_margin=15Plots.mm,
    bottom_margin=10Plots.mm,
    right_margin=5Plots.mm,
    top_margin=5Plots.mm
)

for (i, method) in enumerate(methods_str)
    memory = memory_data[method]
    valid_indices = .!isnan.(memory) .& (memory .> 0)
    if sum(valid_indices) > 0
        plot!(p2, sizes[valid_indices], memory[valid_indices],
            label=method,
            linewidth=2.5,
            marker=markers[i],
            markersize=6,
            color=colors[i],
            linestyle=styles[i]
        )
    end
end

savefig(p2, "memory_comparison.png")
println("  Zapisano: memory_comparison.png")

# === WYKRES 3: BŁĄD WZGLĘDNY ===
println("Generowanie wykresu błędu względnego...")
p3 = plot(
    title="Błąd względny w zależności od rozmiaru problemu",
    xlabel="Rozmiar problemu (n)",
    ylabel="Błąd względny",
    xscale=:log10,
    yscale=:log10,
    legend=:topright,
    size=(1200, 700),
    dpi=300,
    grid=true,
    framestyle=:box,
    left_margin=15Plots.mm,
    bottom_margin=10Plots.mm,
    right_margin=5Plots.mm,
    top_margin=5Plots.mm
)

for (i, method) in enumerate(methods_str)
    errors = error_data[method]
    valid_indices = .!isnan.(errors) .& (errors .> 0)
    if sum(valid_indices) > 0
        plot!(p3, sizes[valid_indices], errors[valid_indices],
            label=method,
            linewidth=2.5,
            marker=markers[i],
            markersize=6,
            color=colors[i],
            linestyle=styles[i]
        )
    end
end

savefig(p3, "error_comparison.png")
println("  Zapisano: error_comparison.png")

# === WYKRES ZŁOŻONY: WSZYSTKO NA JEDNYM ===
println("Generowanie wykresu złożonego...")
p_combined = plot(p1, p2, p3, 
    layout=(3, 1), 
    size=(1200, 1800),
    dpi=300,
    left_margin=15Plots.mm,
    bottom_margin=10Plots.mm,
    right_margin=5Plots.mm,
    top_margin=5Plots.mm
)

savefig(p_combined, "all_comparisons.png")
println("  Zapisano: all_comparisons.png")

# === TABELA PODSUMOWUJĄCA ===
println("\n" * "="^80)
println("PODSUMOWANIE WYNIKÓW")
println("="^80)

for size in sizes
    println("\nRozmiar n = $size:")
    println("-"^80)
    @printf "%-25s | %-12s | %-12s | %-15s\n" "Metoda" "Czas [s]" "Pamięć [MB]" "Błąd względny"
    println("-"^80)
    
    for method in methods_str
        subset = df[(string.(df.Method) .== method) .& (df.Size .== size), :]
        if nrow(subset) > 0
            time = subset[1, :Time_s]
            mem = subset[1, :Memory_MB]
            err = subset[1, :RelativeError]
            @printf "%-25s | %-12s | %-12s | %-15s\n" method time mem err
        end
    end
end

println("\n" * "="^80)
println("ANALIZA ZŁOŻONOŚCI")
println("="^80)

# Dla każdej metody oblicz nachylenie w skali log-log (przybliżona złożoność)
for method in methods_str
    times = time_data[method]
    valid_indices = .!isnan.(times) .& (times .> 0)
    
    if sum(valid_indices) >= 2
        log_sizes = log10.(sizes[valid_indices])
        log_times = log10.(times[valid_indices])
        
        # Prosta regresja liniowa
        n = length(log_sizes)
        mean_x = sum(log_sizes) / n
        mean_y = sum(log_times) / n
        
        numerator = sum((log_sizes .- mean_x) .* (log_times .- mean_y))
        denominator = sum((log_sizes .- mean_x).^2)
        
        slope = numerator / denominator
        
        println("\n$method:")
        @printf "  Złożoność czasowa: O(n^%.2f)\n" slope
        
        # Dla pamięci
        memory = memory_data[method]
        valid_mem = .!isnan.(memory) .& (memory .> 0)
        if sum(valid_mem) >= 2
            log_memory = log10.(memory[valid_mem])
            log_sizes_mem = log10.(sizes[valid_mem])
            
            n_mem = length(log_sizes_mem)
            mean_x_mem = sum(log_sizes_mem) / n_mem
            mean_y_mem = sum(log_memory) / n_mem
            
            num_mem = sum((log_sizes_mem .- mean_x_mem) .* (log_memory .- mean_y_mem))
            den_mem = sum((log_sizes_mem .- mean_x_mem).^2)
            
            slope_mem = num_mem / den_mem
            @printf "  Złożoność pamięciowa: O(n^%.2f)\n" slope_mem
        end
    end
end

println("\n" * "="^80)
println("Wykresy zostały wygenerowane pomyślnie!")
println("="^80)