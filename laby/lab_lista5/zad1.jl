# Tomasz Niedziałek 279754

# ZADANIE 1A: Rozwiązywanie układu Ax = b bez wyboru elementu głównego
# ZŁOŻONOŚĆ CAŁKOWITA:
# - CZAS: O(v * l³) = O(n * l²) gdzie n = v*l
#   Dla stałego l=4: O(n)
# - PAMIĘĆ: O(n) dla kopii wektora + O(v * l²) dla kopii bloków A = O(n*l) = O(n) dla stałego l
function solve_no_pivot(M_orig::BlockMatrix{T}, b::Vector{T}) where T
    x = copy(b)  # O(n) czas i pamięć
    As = deepcopy(M_orig.A)  # O(v * l²) = O(n*l) pamięć
    
    l = M_orig.l
    v = M_orig.v
    
    #   ETAP 1: Eliminacja w przód  
    # Pętla główna: v iteracji
    for k in 1:v
        # 1. Rozkład LU na bieżącym bloku diagonalnym
        # CZAS: O(l³), PAMIĘĆ: O(1) (in-place)
        lu_dense!(As[k])
        
        # 2. Aktualizacja bieżącego fragmentu wektora b (x)
        range_k = (k-1)*l+1 : k*l
        # CZAS: O(l²), PAMIĘĆ: O(1) (in-place)
        solve_L!(As[k], view(x, range_k))
        
        # 3. Jeśli to nie ostatni blok, eliminujemy blok B pod spodem
        if k < v
            # Obliczamy macierz W = A_k^{-1} * C_k
            W = compute_W_matrix(As[k], M_orig.C[k].diag)
            
            # aktualizujemy następny blok diagonalny: A_{k+1} -= B_{k+1} * W
            mul_sub_B_W!(As[k+1], M_orig.B[k], W)
            y = x[range_k]
            z = copy(y)  # O(l) pamięć
            solve_U!(As[k], z)  # O(l²) czas
            
            range_next = k*l+1 : (k+1)*l
            B_next = M_orig.B[k]
            
            # Mnożenie rzadkie B * z - O(l) czas dzięki rzadkości B
            val = dot(B_next.top_row, z)
            x[range_next[1]] -= val
            
            last_z = z[end]
            for i in 2:l
                x[range_next[i]] -= B_next.right_column[i] * last_z
            end
        end
    end
    # SUMA ETAP 1: v * [O(l³) + O(l³) + O(l²) + O(l²)] = O(v * l³) = O(n * l²)
    
    #   ETAP 2: Podstawienie wstecz  
    # Ostatni blok: O(l²) czas
    solve_U!(As[v], view(x, (v-1)*l+1 : v*l))
    
    # Pętla w górę: v-1 iteracji
    for k in v-1:-1:1
        range_k = (k-1)*l+1 : k*l
        range_next = k*l+1 : (k+1)*l
        
        # CZAS: O(l) + O(l²) + O(l) + O(l²) = O(l²)
        # PAMIĘĆ: O(l)
        C_curr = M_orig.C[k]
        current_x_next = x[range_next]
        
        correction = zeros(T, l)  # O(l) pamięć
        for i in 1:l  # O(l) czas
            correction[i] = C_curr.diag[i] * current_x_next[i]
        end
        
        solve_L!(As[k], correction)  # O(l²) czas
        x[range_k] .-= correction  # O(l) czas
        solve_U!(As[k], view(x, range_k))  # O(l²) czas
    end
    # SUMA ETAP 2: v * O(l²) = O(n * l)
    
    return x
end

# ZADANIE 1B: Rozwiązywanie układu Ax = b Z CZĘŚCIOWYM WYBOREM elementu
# ZŁOŻONOŚĆ CAŁKOWITA: identyczna jak solve_no_pivot
# - CZAS: O(n * l²) - dla stałego l=4: O(n)
# - PAMIĘĆ: O(n * l) + O(v * l) dla permutacji = O(n * l)
function solve_pivot(M_orig::BlockMatrix{T}, b::Vector{T}) where T
    x = copy(b)  # O(n) czas i pamięć
    As = deepcopy(M_orig.A)  # O(n*l) pamięć
    l = M_orig.l
    v = M_orig.v
    
    perms = Vector{Vector{Int}}(undef, v)  # O(v*l) = O(n) pamięć
    
    #   ETAP 1 
    for k in 1:v
        # O(l³) czas, O(l) pamięć dla permutacji
        perms[k] = lu_dense_pivot!(As[k])
        
        range_k = (k-1)*l+1 : k*l
        # O(l²) czas, O(l) pamięć dodatkowa dla kopii w permutacji
        solve_L!(As[k], view(x, range_k), perms[k])
        
        if k < v
            # O(l³) czas, O(l²) pamięć
            W = compute_W_matrix(As[k], M_orig.C[k].diag, perms[k])
            
            # O(l²) czas
            mul_sub_B_W!(As[k+1], M_orig.B[k], W)
            
            # O(l²) czas, O(l) pamięć
            y = copy(x[range_k])
            z = copy(y)
            solve_U!(As[k], z)
            
            range_next = k*l+1 : (k+1)*l
            B_next = M_orig.B[k]
            
            val = dot(B_next.top_row, z)
            x[range_next[1]] -= val
            
            last_z = z[end]
            for i in 2:l
                x[range_next[i]] -= B_next.right_column[i] * last_z
            end
        end
    end
    
    #   ETAP 2  
    solve_U!(As[v], view(x, (v-1)*l+1 : v*l))
    
    for k in v-1:-1:1
        range_k = (k-1)*l+1 : k*l
        range_next = k*l+1 : (k+1)*l
        
        # O(l²) czas, O(l) pamięć
        C_curr = M_orig.C[k]
        current_x_next = x[range_next]
        
        correction = zeros(T, l)
        for i in 1:l
            correction[i] = C_curr.diag[i] * current_x_next[i]
        end
        
        # Dodatkowa permutacja: O(l) czas i pamięć
        p = perms[k]
        correction_perm = zeros(T, l)
        for i in 1:l
            correction_perm[i] = correction[p[i]]
        end
        
        solve_L!(As[k], correction_perm)
        x[range_k] .-= correction_perm
        solve_U!(As[k], view(x, range_k))
    end
    
    return x
end