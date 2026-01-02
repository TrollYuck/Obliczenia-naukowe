# funkcje dla bloków

# Zwykła eliminacja Gaussa dla macierzy gęstej (bez wyboru elementu)
# Modyfikuje macierz A in-place: L ląduje pod diagonalą, U na i nad diagonalą.
# ZŁOŻONOŚĆ: O(l³) czas, O(1) pamięć (in-place)
function lu_dense!(A::Matrix{T}) where T
    n = size(A, 1)
    for k in 1:n-1
        for i in k+1:n
            if abs(A[k,k]) < eps(T)
                error("Dzielie przez zero! (Zbyt mały pivot w bloku)")
            end
            factor = A[i, k] / A[k, k]
            A[i, k] = factor 
            for j in k+1:n
                A[i, j] -= factor * A[k, j]
            end
        end
    end
end

# Eliminacja Gaussa z częściowym wyborem elementu głównego
# Zwraca wektor permutacji p
# ZŁOŻONOŚĆ: O(l³) czas, O(l) pamięć (wektor permutacji)
function lu_dense_pivot!(A::Matrix{T}) :: Vector{T} where T 
    n = size(A, 1)
    p = collect(1:n)  # O(l) pamięć
    
    for k in 1:n-1
        # Szukanie max elementu w kolumnie k
        # O(l) czas
        max_val = abs(A[k,k])
        max_idx = k
        for i in k+1:n
            if abs(A[i,k]) > max_val
                max_val = abs(A[i,k])
                max_idx = i
            end
        end
        
        if max_val < eps(T)
            error("Macierz osobliwa w bloku!")
        end
        
        # Zamiana wierszy jeśli znaleziono lepszy pivot
        # O(l) czas
        if max_idx != k
            for col in 1:n
                A[k, col], A[max_idx, col] = A[max_idx, col], A[k, col]
            end
            p[k], p[max_idx] = p[max_idx], p[k]
        end
        
        # Standardowa eliminacja
        # O(l²) czas dla jednego kroku k
        for i in k+1:n
            factor = A[i, k] / A[k, k]
            A[i, k] = factor
            for j in k+1:n
                A[i, j] -= factor * A[k, j]
            end
        end
    end
    return p
end

# Rozwiązywanie układu L * x = b (Forward substitution)
# Zakłada, że L jest w dolnym trójkącie A_lu
# ZŁOŻONOŚĆ: O(l²) czas, O(l) pamięć (jeśli permutacja)
function solve_L!(A_lu::Matrix{T}, b::AbstractVector{T}, p::Vector{Int}=Int[]) where T
    n = size(A_lu, 1)
    # Jeśli podano permutację, najpierw permutujemy b
    # O(l) czas, O(l) pamięć
    if !isempty(p)
        b_copy = copy(b)
        for i in 1:n
            b[i] = b_copy[p[i]]
        end
    end
    
    # O(l²) czas
    for i in 2:n
        sum_val = zero(T)
        for j in 1:i-1
            sum_val += A_lu[i, j] * b[j]
        end
        b[i] -= sum_val
    end
end

# Rozwiązywanie układu U * x = b (Backward substitution)
# Zakłada, że U jest w górnym trójkącie A_lu
# ZŁOŻONOŚĆ: O(l²) czas, O(1) pamięć (in-place)
function solve_U!(A_lu::Matrix{T}, b::AbstractVector{T}) where T
    n = size(A_lu, 1)
    for i in n:-1:1
        sum_val = zero(T)
        for j in i+1:n
            sum_val += A_lu[i, j] * b[j]
        end
        b[i] = (b[i] - sum_val) / A_lu[i, i]
    end
end

# Obliczenie macierzy W = A^{-1} * C
# A_lu to macierz po rozkładzie LU
# ZŁOŻONOŚĆ: O(l³) czas (l razy solve O(l²)), O(l²) pamięć (macierz W)
function compute_W_matrix(A_lu::Matrix{T}, C_diag::Vector{T}, p::Vector{Int}=Int[]) where T
    l = size(A_lu, 1)
    W = zeros(T, l, l)  # O(l²) pamięć
    
    # Pętla po kolumnach: O(l) iteracji
    for col in 1:l
        rhs = zeros(T, l)  # O(l) pamięć (wielokrotnie alokowana)
        rhs[col] = C_diag[col]
        
        # A * column = rhs
        # O(l²) czas każde
        solve_L!(A_lu, rhs, p)
        solve_U!(A_lu, rhs)
        
        # Zapisujemy wynik do kolumny macierzy W
        for row in 1:l
            W[row, col] = rhs[row]
        end
    end
    return W
end

# Aktualizacja następnego bloku: A_next = A_next - B * W
# Wykorzystuje rzadką strukturę macierzy B
# ZŁOŻONOŚĆ: O(l²) czas, O(1) pamięć
function mul_sub_B_W!(A_next::Matrix{T}, B::BorderMatrix{T}, W::Matrix{T}) where T
    l = B.l
    # O(l²) czas - tylko l² operacji zamiast l³ dzięki rzadkości B
    for col in 1:l
        # Wiersz 1 macierzy B (iloczyn skalarny) - O(l) czas
        val = zero(T)
        for k in 1:l
            val += B.top_row[k] * W[k, col]
        end
        A_next[1, col] -= val
        
        # Pozostałe wiersze B (tylko ostatnia kolumna niezerowa)
        # O(l) czas
        for row in 2:l
            val = B.right_column[row] * W[l, col]
            A_next[row, col] -= val
        end
    end
end

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
    
    # === ETAP 1: Eliminacja w przód ===
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
            # a) Obliczamy macierz W = A_k^{-1} * C_k
            # CZAS: O(l³), PAMIĘĆ: O(l²)
            W = compute_W_matrix(As[k], M_orig.C[k].diag)
            
            # b) Aktualizujemy następny blok diagonalny: A_{k+1} -= B_{k+1} * W
            # CZAS: O(l²) - KLUCZOWE: wykorzystuje rzadkość B!
            # PAMIĘĆ: O(1)
            mul_sub_B_W!(As[k+1], M_orig.B[k], W)
            
            # c) Aktualizujemy następny fragment wektora b
            # CZAS: O(l²) + O(l²) + O(l) = O(l²)
            # PAMIĘĆ: O(l)
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
    
    # === ETAP 2: Podstawienie wstecz ===
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
    
    # === ETAP 1 === (analogiczny do solve_no_pivot)
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
    
    # === ETAP 2 ===
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
