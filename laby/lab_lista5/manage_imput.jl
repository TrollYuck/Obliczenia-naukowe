# Tomasz Niedziałek 279754
using LinearAlgebra

function load_matrix_from_path(filename::String) :: BlockMatrix
    open(filename) do f
        # Wczytanie nagłówka
        header = readline(f)
        parts = split(header)
        n = parse(Int, parts[1])
        l = parse(Int, parts[2])
        v = div(n, l) # liczba bloków

        # Pre-alokacja
        # As - v macierzy lxl
        As = [zeros(Float64, l, l) for _ in 1:v]
        
        # Bs - v macierzy typu BorderMatrix 
        Bs = [BorderMatrix(zeros(l), zeros(l)) for _ in 1:v]
        
        # Cs - v macierzy diagonalnych
        Cs = [Diagonal(zeros(l)) for _ in 1:v]

        M = BlockMatrix(As, Bs, Cs, l, v)

        # Wczytywanie danych
        for line in eachline(f)
            
            parts = split(line)
            i, j, val = parse(Int, parts[1]), parse(Int, parts[2]), parse(Float64, parts[3])

            # Numer rzędu blokowego (1..v)
            block_row = div(i - 1, l) + 1
            
            # indeksy lokalne
            local_i = (i - 1) % l + 1
            local_j = (j - 1) % l + 1

            min_j_A = (block_row - 1) * l + 1
            max_j_A = block_row * l

            if j >= min_j_A && j <= max_j_A
                M.A[block_row][local_i, local_j] = val
                
            elseif j > max_j_A
                if local_i == local_j
                M.C[block_row][local_i, local_i] = val
                end
                
            elseif j < min_j_A
                if local_i == 1
                    M.B[block_row].top_row[local_j] = val
                end
                if local_j == l
                    M.B[block_row].right_column[local_i] = val
                end
            end
        end

        return M
    end #file
end #load_matrix_from_path

function load_vector_from_path(filename::String) :: Vector
    open(filename) do f
        size_n = readline(f)
        n = parse(Int,size_n)
        b = zeros(Float64, n)
        index = 1
        for _ in eachline(f)
            temp = readline(f)
            b[index] = parse(Float64, temp)
            index += 1
        end
        return b
    end #file
end #load_vector_from_path

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