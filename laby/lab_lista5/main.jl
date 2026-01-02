using LinearAlgebra
using Printf

if !isdefined(Main, :blocksys)
    include("blocksys.jl")
end
if !isdefined(Main, :matrixgen)
    include("blockmat.jl")
end

using .blocksys
using .matrixgen


function print_usage()
    println("Użycie: julia main.jl [opcje]")
    println("\nOpcje wejścia (Macierz A):")
    println("  -A <plik>       : Ścieżka do pliku z macierzą A")
    println("  -G              : Wygeneruj losową macierz A, wymaga -n, -l, -ck")
    println("     -n <int>     : Rozmiar macierzy n")
    println("     -l <int>     : Rozmiar bloku l")
    println("     -ck <float>  : Współczynnik uwarunkowania, domyślnie 10.0")
    
    println("\nOpcje wejścia (Wektor b):")
    println("  -b <plik>       : Ścieżka do pliku z wektorem b")
    println("  -g              : Wygeneruj b na podstawie x=(1...1) oblicza błąd")
    
    println("\nOpcje działania:")
    println("  (brak)          : Uruchamia wszystkie 4 metody i porównuje wyniki (domyślne)")
    println("  -m <metoda>     : Uruchom jedną konkretną metodę i zapisz wynik.")
    println("     Metody: gauss_no_pivot, gauss_pivot, lu_no_pivot, lu_pivot")
    println("  -test-gauss     : Porównaj warianty metody Gaussa")
    println("  -test-lu        : Porównaj warianty metody LU")
    
    println("\nOpcje wyjścia:")
    println("  -o <plik>       : Plik wynikowy (wymagane tylko przy opcji -m)")
    
    println("\nPrzykłady:")
    println("  julia main.jl -G -n 1000 -l 4 -g")
    println("  julia main.jl -A dane/A.txt -g -m gauss_pivot -o wynik.txt")
end

function parse_command_line()
    # ZMIANA: Domyślny tryb to :all, metoda to nothing
    args = Dict(
        :matrix_file => nothing,
        :vector_file => nothing,
        :output_file => nothing,
        :gen_matrix => false,
        :gen_rhs => false,
        :n => 0,
        :l => 0,
        :ck => 10.0,
        :mode => :all, # Domyślnie uruchamiamy wszystko
        :method => nothing
    )

    i = 1
    while i <= length(ARGS)
        arg = ARGS[i]
        if arg == "-A"
            args[:matrix_file] = ARGS[i+1]; i += 1
        elseif arg == "-b"
            args[:vector_file] = ARGS[i+1]; i += 1
        elseif arg == "-o"
            args[:output_file] = ARGS[i+1]; i += 1
        elseif arg == "-G"
            args[:gen_matrix] = true
        elseif arg == "-g"
            args[:gen_rhs] = true
        elseif arg == "-n"
            args[:n] = parse(Int, ARGS[i+1]); i += 1
        elseif arg == "-l"
            args[:l] = parse(Int, ARGS[i+1]); i += 1
        elseif arg == "-ck"
            args[:ck] = parse(Float64, ARGS[i+1]); i += 1
        elseif arg == "-m"
            # ZMIANA: Podanie metody przełącza tryb na pojedynczy
            args[:method] = ARGS[i+1]; args[:mode] = :single; i += 1
        elseif arg == "-test-gauss"
            args[:mode] = :test_gauss
        elseif arg == "-test-lu"
            args[:mode] = :test_lu
        elseif arg == "-h" || arg == "--help"
            print_usage()
            exit(0)
        else
            println("Nieznany argument: $arg")
            print_usage()
            exit(1)
        end
        i += 1
    end
    return args
end

function save_result(filename::String, x::Vector, relative_error=nothing)
    open(filename, "w") do f
        if !isnothing(relative_error)
            println(f, relative_error)
        end
        for val in x
            println(f, val)
        end
    end
end

function run_single_method(method_name, M, b, x_exact)
    M_work = deepcopy(M)
    
    GC.gc() 
    
    stats = @timed begin
        if method_name == :gauss_no_pivot
            solve_no_pivot(M_work, b)
        elseif method_name == :gauss_pivot
            solve_pivot(M_work, b)
        elseif method_name == :lu_no_pivot
            compute_lu_no_pivot!(M_work)
            solve_from_lu_no_pivot(M_work, b)
        elseif method_name == :lu_pivot
            perms = compute_lu_pivot!(M_work)
            solve_from_lu_pivot(M_work, perms, b)
        end
    end
    
    x_calc = stats.value
    time_sec = stats.time
    memory_mb = stats.bytes / 1024 / 1024
    
    rel_err = NaN
    if !isnothing(x_exact)
        rel_err = norm(x_calc - x_exact) / norm(x_exact)
    else
        res_vec = M * x_calc - b
        rel_err = norm(res_vec) / norm(b)
    end
    
    return time_sec, memory_mb, rel_err, x_calc
end

function main()
    if isempty(ARGS)
        print_usage()
        return
    end

    config = parse_command_line()

    if config[:mode] == :single && isnothing(config[:output_file])
        println("BŁĄD: Nie podano pliku wyjściowego (-o) wymaganego w trybie pojedynczym.")
        exit(1)
    end

    # --- Wczytywanie Macierzy ---
    matrix_file = config[:matrix_file]
    temp_matrix_file = nothing

    if config[:gen_matrix]
        if config[:n] == 0 || config[:l] == 0
            println("BŁĄD: Do generowania macierzy (-G) wymagane są -n i -l.")
            exit(1)
        end
        println("Generowanie losowej macierzy $(config[:n])x$(config[:n])...")
        temp_matrix_file = "temp_A_generated.txt"
        matrixgen.blockmat(config[:n], config[:l], config[:ck], temp_matrix_file)
        matrix_file = temp_matrix_file
    elseif isnothing(matrix_file)
        println("BŁĄD: Brak źródła macierzy (-A lub -G).")
        exit(1)
    end

    println("Wczytywanie macierzy z: $matrix_file")
    M = load_matrix_from_path(matrix_file)

    # --- Wczytywanie Wektora ---
    b = Vector{Float64}()
    x_exact = nothing

    if config[:gen_rhs]
        println("Generowanie wektora b na podstawie x = ones...")
        n_total = M.v * M.l
        x_exact = ones(n_total)
        b = M * x_exact
    elseif !isnothing(config[:vector_file])
        println("Wczytywanie wektora b z: $(config[:vector_file])")
        b = load_vector_from_path(config[:vector_file])
    else
        println("BŁĄD: Brak źródła wektora b (-b lub -g).")
        exit(1)
    end
    
    # --- Wykonanie (Testy zbiorcze lub pojedyncze) ---
    if config[:mode] in [:all, :test_gauss, :test_lu]
        
        println("\n" * "="^80) # Dłuższa linia dla nowej kolumny
        methods = []
        
        if config[:mode] == :test_gauss
            println(" TESTOWANIE: Metoda Eliminacji Gaussa")
            methods = [(:gauss_no_pivot, "Gauss bez wyboru"), (:gauss_pivot, "Gauss z wyborem")]
        elseif config[:mode] == :test_lu
            println(" TESTOWANIE: Rozkład LU + Rozwiązanie")
            methods = [(:lu_no_pivot, "LU bez wyboru"), (:lu_pivot, "LU z wyborem")]
        else 
            println(" TESTOWANIE: Wszystkie metody")
            methods = [
                (:gauss_no_pivot, "Gauss bez wyboru"), 
                (:gauss_pivot, "Gauss z wyborem"),
                (:lu_no_pivot, "LU bez wyboru"), 
                (:lu_pivot, "LU z wyborem")
            ]
        end
        
        println("="^80)
        # Zaktualizowany nagłówek i formatowanie
        @printf "%-25s | %-12s | %-12s | %-15s\n" "Wariant" "Czas [s]" "Pamięć [MB]" "Błąd względny"
        println("-"^80)

        for (sym, name) in methods
            t, mem, err, _ = run_single_method(sym, M, b, x_exact)
            # Dodano zmienną mem do wydruku
            @printf "%-25s | %-12.6f | %-12.2f | %-15.4e\n" name t mem err
        end
        println("-"^80)
        println("Pamięć i czas zawierają narzut alokacji kopii macierzy.")

    else
        # Tryb pojedynczy
        method_sym = Symbol(config[:method])
        println("Uruchamianie metody: $(config[:method])...")
        t, mem, err, x_res = run_single_method(method_sym, M, b, x_exact)
        
        println("Czas: $t s")
        println("Pamięć: $(round(mem, digits=2)) MB") # Tutaj też dodano wyświetlanie
        if !isnan(err)
            println("Błąd względny: $err")
        end
        
        println("Zapisywanie wyniku do: $(config[:output_file])")
        save_result(config[:output_file], x_res, isnan(err) ? nothing : err)
    end

    if !isnothing(temp_matrix_file)
        rm(temp_matrix_file)
    end
end

main()