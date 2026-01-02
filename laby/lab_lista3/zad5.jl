# Tomasz Niedziałek 279754

include("mbisekcji.jl")

h(x) = 3*x
g(x) = (Base.MathConstants.e)^x 

f(x) = h(x) - g(x)

delta, epsilon = 10^(-4), 10^(-4)

println("-"^50)

# Punkt 1: Szukamy w przedziale [0, 1]
a1, b1 = 0.0, 1.0
(x1, val1, it1, err1) = mbisekcji(f, a1, b1, delta, epsilon)

if err1 == 0
    println("Punkt przecięcia nr 1: x ≈ $x1")
    println("Wartość różnicy funkcji: $val1")
    println("Sprawdzenie:")
    println("   3x = $(3*x1)")
    println("   e^x = $(exp(x1))")
    println("Iteracja: $it1")
else
    println("Błąd w szukaniu punktu 1 (kod błędu: $err1)")
end

println("-"^50)

# Punkt 2: Szukamy w przedziale [1, 2]
a2, b2 = 1.0, 2.0
(x2, val2, it2, err2) = mbisekcji(f, a2, b2, delta, epsilon)

if err2 == 0
    println("Punkt przecięcia nr 2: x ≈ $x2")
    println("Wartość różnicy funkcji: $val2")
    println("Sprawdzenie:")
    println("   3x = $(3*x2)")
    println("   e^x = $(exp(x2))")
    println("Iteracja: $it2")
else
    println("Błąd w szukaniu punktu 2 (kod błędu: $err2)")
end