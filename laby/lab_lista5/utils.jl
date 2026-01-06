# Tomasz Niedziałek 279754
# funkcje dla bloków

# Rozkład LU bloku bez wyboru pivota
# Modyfikuje macierz A in-place: L ląduje pod diagonalą, U na i nad diagonalą.
# ZŁOŻONOŚĆ: O(l³) czas, O(1) pamięć (in-place)
function lu_dense!(A::Matrix{T}) where T
    n = size(A, 1)
    for k in 1:n-1
        for i in k+1:n
            if abs(A[k,k]) < eps(T)
                error("Dzielie przez zero (Zbyt mały pivot w bloku)")
            end
            factor = A[i, k] / A[k, k]
            A[i, k] = factor 
            for j in k+1:n
                A[i, j] -= factor * A[k, j]
            end
        end
    end
end

# Rozkład LU bloku z wyborem pivota
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
    for col in 1:l
        # Wiersz 1 macierzy B (iloczyn skalarny) 
        val = zero(T)
        for k in 1:l
            val += B.top_row[k] * W[k, col]
        end
        A_next[1, col] -= val
        
        # Pozostałe wiersze B (tylko ostatnia kolumna niezerowa)
        for row in 2:l
            val = B.right_column[row] * W[l, col]
            A_next[row, col] -= val
        end
    end
end