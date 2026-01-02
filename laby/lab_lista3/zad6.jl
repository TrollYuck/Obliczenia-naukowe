# Tomasz Niedziałek 279754
using Printf

include("mbisekcji.jl")
include("mstycznych.jl")
include("msiecznych.jl")

# Funkcja 1: f1(x) = e^(1-x) - 1
# Miejsce zerowe w x = 1
f1(x) = exp(1 - x) - 1.0
pf1(x) = -exp(1 - x)  # Pochodna

# Funkcja 2: f2(x) = x * e^(-x)
# Miejsce zerowe w x = 0
f2(x) = x * exp(-x)
pf2(x) = exp(-x) * (1 - x) # Pochodna: e^-x - x*e^-x

delta = 1e-5
epsilon = 1e-5
maxit = 40

# CZĘŚĆ 1: Funkcja f1(x) = e^(1-x) - 1
println("\n" * "="^60)
println("ANALIZA FUNKCJI f1(x) = e^(1-x) - 1 (Pierwiastek w x=1)")
println("="^60)

# 1. Bisekcja: Szukamy w przedziale [0, 2], bo f(0) > 0, f(2) < 0
r_b1, v_b1, it_b1, err_b1 = mbisekcji(f1, 0.0, 2.0, delta, epsilon)

# 2. Newton: Szukamy blisko, np. x0 = 0.5
r_n1, v_n1, it_n1, err_n1 = mstycznych(f1, pf1, 0.5, delta, epsilon, maxit)

# 3. Sieczne: Punkty 0 i 2
r_s1, v_s1, it_s1, err_s1 = msiecznych(f1, 0.0, 2.0, delta, epsilon, maxit)

@printf("Metoda     | %-12s | %-12s | It | Err\n", "Wynik x", "Wartość f(x)")
println("-"^50)
@printf("Bisekcja   | %12.8f | %12.1e | %2d | %d\n", r_b1, v_b1, it_b1, err_b1)
@printf("Newton     | %12.8f | %12.1e | %2d | %d\n", r_n1, v_n1, it_n1, err_n1)
@printf("Siecznych  | %12.8f | %12.1e | %2d | %d\n", r_s1, v_s1, it_s1, err_s1)

println("\n--- Test specjalny dla f1: Newton z x0 > 1 (np. x0 = 4.0) ---")
# Pytanie z zadania: co się stanie dla x0 in (1, infinity)?
# Odp: Funkcja jest monotoniczna, wypukła/wklęsła w odpowiedni sposób, powinna zbiegać.
x0_test = 5.0
r_test, v_test, it_test, err_test = mstycznych(f1, pf1, x0_test, delta, epsilon, maxit)
@printf("Start %.1f -> Wynik: %.8f (Iteracje: %d, Err: %d)\n", x0_test, r_test, it_test, err_test)

# CZĘŚĆ 2: Funkcja f2(x) = x * e^(-x)
println("\n" * "="^60)
println("ANALIZA FUNKCJI f2(x) = x * e^(-x) (Pierwiastek w x=0)")
println("="^60)

# 1. Bisekcja: Szukamy w [-1, 0.5]. f(-1) < 0, f(0.5) > 0
r_b2, v_b2, it_b2, err_b2 = mbisekcji(f2, -1.0, 0.5, delta, epsilon)

# 2. Newton: Startujemy blisko zera, np. x0 = -0.5
r_n2, v_n2, it_n2, err_n2 = mstycznych(f2, pf2, -0.5, delta, epsilon, maxit)

# 3. Sieczne: Punkty -1 i 0.5
r_s2, v_s2, it_s2, err_s2 = msiecznych(f2, -1.0, 0.5, delta, epsilon, maxit)

@printf("Metoda     | %-12s | %-12s | It | Err\n", "Wynik x", "Wartość f(x)")
println("-"^50)
@printf("Bisekcja   | %12.8f | %12.1e | %2d | %d\n", r_b2, v_b2, it_b2, err_b2)
@printf("Newton     | %12.8f | %12.1e | %2d | %d\n", r_n2, v_n2, it_n2, err_n2)
@printf("Siecznych  | %12.8f | %12.1e | %2d | %d\n", r_s2, v_s2, it_s2, err_s2)


# TESTY SPECJALNE DLA f2 

println("\n--- Test specjalny A: Newton dla f2 z x0 > 1 (np. x0 = 2.0) ---")
# Dla x > 1 pochodna staje się ujemna, a styczna przecina oś X daleko w prawo.
x0_bad = 20.0
r_bad, v_bad, it_bad, err_bad = mstycznych(f2, pf2, x0_bad, delta, epsilon, maxit)
@printf("Start %.1f -> Wynik: %.8f (Iteracje: %d, Err: %d)\n", x0_bad, r_bad, it_bad, err_bad)

println("\n--- Test specjalny B: Newton dla f2 z x0 = 1 ---")
# W punkcie x = 1 funkcja ma ekstremum (maksimum), pochodna wynosi 0.
x0_zero_deriv = 1.0
r_zero, v_zero, it_zero, err_zero = mstycznych(f2, pf2, x0_zero_deriv, delta, epsilon, maxit)
@printf("Start %.1f -> Wynik: %.8f (Iteracje: %d, Err: %d)\n", x0_zero_deriv, r_zero, it_zero, err_zero)