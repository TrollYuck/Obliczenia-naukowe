# Tomasz Niedziałek 279754

# ZADANIE 2A: Faktoryzacja LU bez wyboru elementu
# ZŁOŻONOŚĆ CAŁKOWITA:
# - CZAS: O(v * l³) = O(n * l²) - dla stałego l: O(n)
# - PAMIĘĆ: O(1) - modyfikuje M in-place
function compute_lu_no_pivot!(M::BlockMatrix{T}) where T
    l = M.l
    v = M.v
    
    # Pętla po v blokach
    for k in 1:v
        # Faktoryzacja bloku diagonalnego
        # CZAS: O(l³), PAMIĘĆ: O(1) (in-place)
        lu_dense!(M.A[k])
        
        # Aktualizacja bloku następnego (Schur complement)
        if k < v
            # W = A_k^{-1} * C_k
            # CZAS: O(l³), PAMIĘĆ: O(l²) - tymczasowa macierz W
            W = compute_W_matrix(M.A[k], M.C[k].diag)
            
            # A_{k+1} = A_{k+1} - B_{k+1} * W
            # CZAS: O(l²) 
            # PAMIĘĆ: O(1)
            mul_sub_B_W!(M.A[k+1], M.B[k], W)
        end
    end
    # SUMA: v * [O(l³) + O(l³) + O(l²)] = O(v * l³) = O(n * l²)
end

# ZADANIE 2B: Faktoryzacja LU z częściowym wyborem elementu
# ZŁOŻONOŚĆ CAŁKOWITA:
# - CZAS: O(n * l²) - dla stałego l: O(n)
# - PAMIĘĆ: O(n) dla wektora permutacji
function compute_lu_pivot!(M::BlockMatrix{T}) where T
    l = M.l
    v = M.v
    perms = Vector{Vector{Int}}(undef, v)  # O(v*l) = O(n) pamięć
    
    for k in 1:v
        # O(l³) czas, O(l) pamięć
        perms[k] = lu_dense_pivot!(M.A[k])
        
        if k < v
            # W uwzględnia permutację p_k
            # CZAS: O(l³), PAMIĘĆ: O(l²) tymczasowa
            W = compute_W_matrix(M.A[k], M.C[k].diag, perms[k])
            
            # CZAS: O(l²), PAMIĘĆ: O(1)
            mul_sub_B_W!(M.A[k+1], M.B[k], W)
        end
    end
    return perms
end