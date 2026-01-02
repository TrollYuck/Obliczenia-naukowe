#

using LinearAlgebra

# struktura macierzy B
struct BorderMatrix{T} <: AbstractMatrix{T}
    top_row::Vector{T}
    right_column::Vector{T}
    l::Int
end

# Konstruktor zewnętrzny
function BorderMatrix(t::Vector{T}, r::Vector{T}) where T
    @assert length(t) == length(r) "Verticies must be of same length"
    return BorderMatrix(t, r, length(t))
end

# Implementacja mnożenia B * x
import Base: *
function *(B::BorderMatrix, x::AbstractVector)
    @assert length(x) == B.l "Given vector and matrix of incompatible dimentions"
    y = zeros(eltype(x), B.l)
    
    # Pierwszy wiersz
    y[1] = dot(B.top_row, x)
    
    # Pozostałe wiersze (tylko ostatnia kolumna niezerowa)
    x_last = x[end]
    for i in 2:B.l
        y[i] = B.right_column[i] * x_last
    end
    
    return y
end

# impementacja podstawowych funcji, żeby Julia widziała B jako macierz
Base.size(B::BorderMatrix) = (B.l, B.l)
function Base.getindex(B::BorderMatrix{T}, i::Int, j::Int) where T
    @boundscheck begin
        @assert i > 0 && i <= B.l "Row Index out of bounds"
        @assert j > 0 && j <= B.l "Column Index out of bounds"
    end

    if i == 1
        return B.top_row[j]
    elseif j == B.l
        return B.right_column[i]
    else 
        return zero(T)
    end
end

# macierz blokowa trójdiagonalna A
struct BlockMatrix{T}
    A::Vector{Matrix{T}}               # Macierze gęste
    B::Vector{BorderMatrix{T}}         # Macierze typu B "odwrócone L"
    C::Vector{Diagonal{T, Vector{T}}}  # Macierze diagonalne
    
    l::Int # rozmiar pojedynczego bloku
    v::Int # liczba bloków
end

# funkcja do wyswietlania BlockMatrix A w terminalu
function Base.show(io::IO, ::MIME"text/plain", M::BlockMatrix)
    n_total = M.v * M.l
    println(io, "$(M.v)×$(M.v) BlockMatrix{$(eltype(M.A[1]))}")
    println(io, "  Total size: $(n_total)×$(n_total)")
    println(io, "  Block size: $(M.l)×$(M.l)")
    
    if M.v <= 10
        println(io, "  Structure:")
        for i in 1:M.v
            row_str = "  "
            for j in 1:M.v
                if i == j
                    row_str *= "[A]"
                elseif i == j + 1
                    row_str *= "[B]"
                elseif i == j - 1
                    row_str *= "[C]"
                else
                    row_str *= " . "
                end
            end
            println(io, row_str)
        end
    end
end 

import Base: *
# Mnożenie macierzy blokowej przez wektor: y = M * x
# Złożoność: O(N), Pamięć: O(N)
function *(M::BlockMatrix{T}, x::Vector{T}) where T
    n = M.v * M.l
    if length(x) != n
        error("Wymiary się nie zgadzają! Macierz: $n x $n, Wektor: $(length(x))")
    end
    
    y = zeros(T, n)
    l = M.l
    v = M.v
    
    for k in 1:v
        range_curr = (k-1)*l+1 : k*l
        
        y_view = view(y, range_curr)
        
        # Główna przekątna (A_k * x_k)
        mul!(y_view, M.A[k], view(x, range_curr))
        
        # Pod przekątną (B_k * x_{k-1})
        if k > 1
            range_prev = (k-2)*l+1 : (k-1)*l
            y_view .+= M.B[k-1] * view(x, range_prev)
        end
        
        # Nad przekątną (C_k * x_{k+1})
        if k < v
            range_next = k*l+1 : (k+1)*l
            y_view .+= M.C[k] * view(x, range_next)
        end
    end
    
    return y
end