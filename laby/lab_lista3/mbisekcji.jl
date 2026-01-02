# Tomasz Niedziałek 279754

function mbisekcji(f, a::Float64, b::Float64, delta::Float64, epsilon::Float64) :: Tuple{Float64,Float64,Int,Int}
    u = f(a)
    v = f(b)
    e = b - a

    if sign(u) == sign(v)
        return (0.0, 0.0, 0, 1) # Błąd: Funkcja nie zmienia znaku
    end
    
    k = 0
    while true
        k += 1
        e = e/2
        c = a + e 
        w = f(c)
        
        # Kryterium zatrzymania: Osiągnięcie wymaganej dokładności przybliżenia pierwiastka (delta) 
        # LUB osiągnięcie wymaganej dokładności wartości funkcji (epsilon)
        if abs(e) <= delta || abs(w) < epsilon 
            return (c, w, k, 0) # Sukces: err = 0
        end
        
        # Aktualizacja przedziału
        if sign(w) != sign(u)
            b = c
            v = w    
        else
            a = c 
            u = w 
        end
    end
end