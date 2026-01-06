# Tomasz Niedziałek 279754

function solve_from_lu_no_pivot(M_lu::BlockMatrix{T}, b::Vector{T}) where T
    x = copy(b)  
    l = M_lu.l
    v = M_lu.v
 
    for k in 1:v
        r_curr = (k-1)*l+1 : k*l
        
        if k > 1
            r_prev = (k-2)*l+1 : (k-1)*l
            y_prev = x[r_prev]
            
            z = copy(y_prev)
            solve_U!(M_lu.A[k-1], z)
            
            B = M_lu.B[k-1]

            val = dot(B.top_row, z)
            x[r_curr[1]] -= val
            
            z_last = z[end]
            for i in 2:l
                x[r_curr[i]] -= B.right_column[i] * z_last
            end
        end
        
        solve_L!(M_lu.A[k], view(x, r_curr))
    end

    solve_U!(M_lu.A[v], view(x, (v-1)*l+1 : v*l))
    

    for k in v-1:-1:1
        r_curr = (k-1)*l+1 : k*l
        r_next = k*l+1 : (k+1)*l
        
        # CZAS: O(l) + O(l²) + O(l) + O(l²) = O(l²)
        # PAMIĘĆ: O(l)
        x_next = x[r_next]
        correction = zeros(T, l)  # O(l) pamięć
        for i in 1:l  # O(l) czas
            correction[i] = M_lu.C[k].diag[i] * x_next[i]
        end
        
        solve_L!(M_lu.A[k], correction)  # O(l²) czas
        x[r_curr] .-= correction  # O(l) czas
        solve_U!(M_lu.A[k], view(x, r_curr))  # O(l²) czas
    end
    # SUMA ETAP 2: v * O(l²) = O(n * l)
    
    return x
end

# ZADANIE 3B: Rozwiązywanie z gotowym LU (z pivotem)
# ZŁOŻONOŚĆ CAŁKOWITA:
# - CZAS: O(n * l) - dla stałego l: O(n)
# - PAMIĘĆ: O(n)
function solve_from_lu_pivot(M_lu::BlockMatrix{T}, perms::Vector{Vector{Int}}, b::Vector{T}) where T
    x = copy(b)  # O(n) czas i pamięć
    l = M_lu.l
    v = M_lu.v
    
    #   ETAP 1: Forward  
    for k in 1:v
        r_curr = (k-1)*l+1 : k*l
        
        if k > 1
            # O(l²) czas, O(l) pamięć
            r_prev = (k-2)*l+1 : (k-1)*l
            y_prev = x[r_prev]
            
            z = copy(y_prev)
            solve_U!(M_lu.A[k-1], z)
            
            B = M_lu.B[k-1]
            val = dot(B.top_row, z)
            x[r_curr[1]] -= val
            z_last = z[end]
            for i in 2:l
                x[r_curr[i]] -= B.right_column[i] * z_last
            end
        end
        
        # używamy permutacji w solve_L
        solve_L!(M_lu.A[k], view(x, r_curr), perms[k])
    end
    
    #   ETAP 2: Backward  
    solve_U!(M_lu.A[v], view(x, (v-1)*l+1 : v*l))
    
    for k in v-1:-1:1
        r_curr = (k-1)*l+1 : k*l
        r_next = k*l+1 : (k+1)*l
        
        # O(l²) czas, O(l) pamięć
        x_next = x[r_next]
        correction = zeros(T, l)
        for i in 1:l
            correction[i] = M_lu.C[k].diag[i] * x_next[i]
        end
        
        # Dodatkowa permutacja: O(l) czas i pamięć
        p = perms[k]
        correction_perm = zeros(T, l)
        for i in 1:l
            correction_perm[i] = correction[p[i]]
        end
        
        solve_L!(M_lu.A[k], correction_perm)
        x[r_curr] .-= correction_perm
        solve_U!(M_lu.A[k], view(x, r_curr))
    end
    
    return x
end
