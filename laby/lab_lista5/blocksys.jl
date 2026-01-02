# 

module blocksys
    export BorderMatrix, BlockMatrix, load_matrix_from_path, load_vector_from_path, save_result,
    lu_dense!, lu_dense_pivot!, solve_L!, solve_U!, solve_no_pivot, solve_pivot, 
    compute_W_matrix, mul_sub_B_W!, compute_lu_pivot!, compute_lu_no_pivot!,
    solve_from_lu_no_pivot, solve_from_lu_pivot
    include("structures.jl")
    include("zad1.jl")
    include("zad2.jl")
    include("zad3.jl")
    include("manage_imput.jl")
end