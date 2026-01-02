#

using Test
using LinearAlgebra

include("blocksys.jl")
include("blockmat.jl")
using .blocksys
using .matrixgen

ATOL = 1e-14

# Funkcja zamieniająca BlockMatrix na zwykłą macierz
function to_dense(M::BlockMatrix{T}) where T
    n = M.v * M.l
    l = M.l
    A_dense = zeros(T, n, n)

    for k in 1:M.v
        range_k = (k-1)*l+1 : k*l
        
        A_dense[range_k, range_k] = M.A[k]

        if k < M.v
            range_next = k*l+1 : (k+1)*l
            C_mat = M.C[k]
            A_dense[range_k, range_next] = C_mat
        end

        if k < M.v
            range_next = k*l+1 : (k+1)*l
            
            B_obj = M.B[k] 
            
            B_mat = zeros(T, l, l)
            B_mat[1, :] = B_obj.top_row
            B_mat[:, end] = B_obj.right_column
            
            A_dense[range_next, range_k] = B_mat
        end
    end
    return A_dense
end

# Generator losowych danych testowych
function generate_random_system(v::Int, l::Int)
    As = [randn(l, l) for _ in 1:v]
    for A in As
        for i in 1:l
            A[i,i] += 20.0 
        end
    end
    Bs = Vector{BorderMatrix{Float64}}(undef, v-1)
    for i in 1:v-1
        t = randn(l)
        r = randn(l)
        # Matematyczna spójność rogu (element wspólny)
        r[1] = t[end]
        Bs[i] = BorderMatrix(t, r)
    end

    Cs = [Diagonal(randn(l)) for _ in 1:v-1]

    M = BlockMatrix(As, Bs, Cs, l, v)
    
    n = v * l
    b = randn(n)
    
    return M, b
end

@testset "Testy Struktur i Pomocnicze" begin
    l = 4
    t = ones(l)
    r = ones(l)
    B = BorderMatrix(t, r)
    x = ones(l)
    y = B * x
    @test y[1] == 4.0
    @test y[2] == 1.0
    @test y[l] == 1.0
    
    @test B[1, 1] == 1.0
    @test B[2, 1] == 0.0 
end

@testset "Zadanie 1: Eliminacja Gaussa" begin
    @testset "Zadanie 1A: Solver bez wyboru elementu głównego (solve_no_pivot)" begin
        v, l = 2, 3
        M, b = generate_random_system(v, l)
        
        A_dense = to_dense(M)
        x_ref = A_dense \ b 
        
        M_test = deepcopy(M)
        
        x_calc = solve_no_pivot(M_test, b)
        
        @test x_calc ≈ x_ref atol=ATOL
    end

    @testset "Zadanie 1B: Solver z wyborem elementu głównego (solve_pivot)" begin
        v, l = 4, 4
        M, b = generate_random_system(v, l)
        
        A_dense = to_dense(M)
        x_ref = A_dense \ b
        
        M_test = deepcopy(M)
        
        x_calc = solve_pivot(M_test, b)
        
        @test x_calc ≈ x_ref atol=ATOL
        
        l = 2
        v = 2
        As = [Float64[0.0 1.0; 1.0 0.0], Float64[1.0 0.0; 0.0 1.0]] 
        t = [1.0, 1.0]; r = [1.0, 1.0]
        Bs = [BorderMatrix(t, r)] 
        Cs = [Diagonal([1.0, 1.0])]
        
        M_tricky = BlockMatrix(As, Bs, Cs, l, v)
        b_tricky = [1.0, 2.0, 3.0, 4.0]
        
        A_dense = to_dense(M_tricky)
        x_ref = A_dense \ b_tricky
        
        M_test = deepcopy(M_tricky)
        x_calc = solve_pivot(M_test, b_tricky)
        
        @test x_calc ≈ x_ref atol=ATOL
    end

    @testset "Test losowej macierzy z .matrixgen (solve_pivot)" begin
        n, l, ck = 8, 2, 5.0
        fname = "tmp_blockmat_test.txt"

        matrixgen.blockmat(n, l, ck, fname)
        M = load_matrix_from_path(fname)
        A_dense = to_dense(M)

        b = randn(n)
        x_ref = A_dense \ b

        M_test = deepcopy(M)
        x_calc = solve_pivot(M_test, b)

        @test x_calc ≈ x_ref atol=ATOL

        try
            rm(fname)
        catch
        end
    end

    @testset "Testy wydajnościowe (Duże N)" begin
        
        # N = 10 000 (l=4, v=2500)

        println("\nTestowanie N = 10 000...")
        n = 10_000
        l = 4
        v = div(n, l)

        M, _ = generate_random_system(v, l)
        x_true = randn(n)


        b = M * x_true
        
        M_test = deepcopy(M)

        x_calc_no_pivot = @time solve_no_pivot(M_test, b)
        
        rel_err_np = norm(x_calc_no_pivot - x_true) / norm(x_true)
        @test rel_err_np < ATOL

        M_test = deepcopy(M)
        x_calc_pivot = @time solve_pivot(M_test, b)
        
        rel_err_p = norm(x_calc_pivot - x_true) / norm(x_true)
        @test rel_err_p < ATOL

        # N = 500 000 (l=4, v=125 000)

        println("\nTestowanie N = 500 000...")
        n = 500_000
        l = 4
        v = div(n, l)
        
        M, _ = generate_random_system(v, l)
        x_true = randn(n)
        b = M * x_true
        
        println("  -> Solve No Pivot:")
        M_test = deepcopy(M) 
        x_calc_no_pivot = @time solve_no_pivot(M_test, b)
        
        rel_err_np = norm(x_calc_no_pivot - x_true) / norm(x_true)
        @test rel_err_np < ATOL 

        println("  -> Solve Pivot:")
        M_test = deepcopy(M)
        x_calc_pivot = @time solve_pivot(M_test, b)
        
        rel_err_p = norm(x_calc_pivot - x_true) / norm(x_true)
        @test rel_err_p < ATOL
    end

    @testset "Testy wydajnościowe" begin

        N, l = 300, 100
        
        for i in 0:4
            n = N * 10^i 
            println("\nTestowanie N=$n, l=$l...")
            v = div(n, l)

            M, _ = generate_random_system(v, l)
            x_true = randn(n)
            b = M * x_true

            M_test = deepcopy(M)
            x_np = @time solve_no_pivot(M_test, b)
            @test norm(x_np - x_true) / norm(x_true) < ATOL

            M_test = deepcopy(M)
            x_p = @time solve_pivot(M_test, b)
            @test norm(x_p - x_true) / norm(x_true) < ATOL
        end

    end
end



@testset "Zadanie 2: Faktoryzacja LU" begin
    l, v = 4, 5

    @testset "LU bez wyboru elementu" begin
        M, b = generate_random_system(v, l)
        M_copy = deepcopy(M)
        
        compute_lu_no_pivot!(M)
        
        @test M.A[1] != M_copy.A[1]
        
        x_calc = solve_from_lu_no_pivot(M, b)
        x_ref = to_dense(M_copy) \ b
        
        @test x_calc ≈ x_ref atol=ATOL
    end

    @testset "LU z częściowym wyborem elementu" begin
        M, b = generate_random_system(v, l)
        M.A[1][1,1] = 0.0
        M.A[1][2,1] = 10.0
        M_copy = deepcopy(M)
        
        perms = compute_lu_pivot!(M)
        
        @test length(perms) == v
        @test perms[1][1] != 1
        
        x_calc = solve_from_lu_pivot(M, perms, b)
        x_ref = to_dense(M_copy) \ b
        
        @test x_calc ≈ x_ref atol=ATOL
    end
end

@testset "Zadanie 3: Rozwiązywanie z LU" begin
    l, v = 4, 10
    M, b = generate_random_system(v, l)
    
    M_lu_np = deepcopy(M)
    compute_lu_no_pivot!(M_lu_np)
    
    M_lu_p = deepcopy(M)
    perms = compute_lu_pivot!(M_lu_p)
    
    @testset "Wielokrotne użycie (No Pivot)" begin
        x1 = randn(v * l)
        b1 = M * x1
        x_c1 = solve_from_lu_no_pivot(M_lu_np, b1)
        @test norm(x_c1 - x1) / norm(x1) < ATOL
        
        x2 = ones(v * l)
        b2 = M * x2
        x_c2 = solve_from_lu_no_pivot(M_lu_np, b2)
        @test norm(x_c2 - x2) / norm(x2) < ATOL
    end

    @testset "Wielokrotne użycie (Pivot)" begin
        x1 = randn(v * l)
        b1 = M * x1
        x_c1 = solve_from_lu_pivot(M_lu_p, perms, b1)
        @test norm(x_c1 - x1) / norm(x1) < ATOL
        
        x2 = ones(v * l)
        b2 = M * x2
        x_c2 = solve_from_lu_pivot(M_lu_p, perms, b2)
        @test norm(x_c2 - x2) / norm(x2) < ATOL
    end
    
    @testset "Wydajność solvera (N=100k)" begin
        n, l = 100_000, 4
        v = div(n, l)
        M_big, b = generate_random_system(v, l)
        
        compute_lu_no_pivot!(M_big)

        println("\nCzas solve_from_lu_no_pivot (N=$n):")
        @time x = solve_from_lu_no_pivot(M_big, b)
        @test length(x) == n
    end
end