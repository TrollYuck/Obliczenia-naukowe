# Tomasz Niedziałek 279754

using Printf

include("mbisekcji.jl")
include("mstycznych.jl")
include("msiecznych.jl")

f(x) = sin(x) - ((1/2)*x)^2
pf(x) = cos(x) - x/2
maxit = 20

# bisekcja
a_bi = 1.5
b_bi = 2.0
delta = (1/2)*10^(-5)
epsilon = (1/2)*10^(-5)
x_bi, fx_bi, it_bi, err_bi = mbisekcji(f, a_bi, b_bi, delta, epsilon)

# Newton
x0_new = 1.5
x_new, fx_new, it_new, err_new = mstycznych(f, pf, x0_new, delta, epsilon, maxit)

# siecznych
x0_sie = 1.0
x1_sie = 2.0
x_sie, fx_sie, it_sie, err_sie = msiecznych(f, x0_sie, x1_sie, delta, epsilon, maxit)


println("\nWyniki obliczeń dla f(x) = sin(x) - (x/2)^2")
println("Wymagana dokładność: delta=$(delta), epsilon=$(epsilon)\n")

line = "-"^78
println(line)

# %-12s - string wyrównany do lewej, szerokość 12
# %-18s - string wyrównany do lewej, szerokość 18
@printf("| %-12s | %-18s | %-18s | %-5s | %-3s |\n", 
        "Metoda", "Przybliżenie (r)", "Wartość f(r)", "Iter", "Err")

println(line)

# %.12f - float z 12 miejscami po przecinku
# %.4e  - notacja naukowa dla bardzo małych liczb
@printf("| %-12s | %18.12f | %18.4e | %5d | %3d |\n", 
        "Bisekcja", x_bi, fx_bi, it_bi, err_bi)

@printf("| %-12s | %18.12f | %18.4e | %5d | %3d |\n", 
        "Newtona", x_new, fx_new, it_new, err_new)

@printf("| %-12s | %18.12f | %18.4e | %5d | %3d |\n", 
        "Siecznych", x_sie, fx_sie, it_sie, err_sie)

println(line)

println("\nLegenda err: 0 - Sukces, 1 - Max iteracji/Brak znaku, 2 - Pochodna=0")