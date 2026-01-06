# Tomasz Niedziałek 279754

module blocksys
    export BorderMatrix, BlockMatrix, 
    load_matrix_from_path, load_vector_from_path, save_result,
    solve_no_pivot, solve_pivot, # Gauss
    compute_lu_pivot!, compute_lu_no_pivot!, # Wyliczenie LU macierzy blokowej
    solve_from_lu_no_pivot, solve_from_lu_pivot # Rozwiązanie po LU
    include("structures.jl")
    include("utils.jl")
    include("zad1.jl")
    include("zad2.jl")
    include("zad3.jl")
    include("manage_imput.jl")
end