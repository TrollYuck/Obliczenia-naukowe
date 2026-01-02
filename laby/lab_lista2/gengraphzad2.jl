#Tomasz Niedziałek 279754
using Plots

gr()

f(x) = exp(x) * log(1 + exp(-x))

x_vals = -20:0.1:50

y_vals = f.(x_vals)

plot(x_vals, y_vals,
    label="f(x) = e^x * ln(1 + e^(-x))", 
    xlabel="x",                          
    ylabel="f(x)",                       
    title="Wykres funkcji",              
    linewidth=2,                         
    framestyle=:origin,
    legend=:bottomright                  
)

hline!([1], label="Granica y=1", linestyle=:dash, color=:red)

savefig("wykres_funkcji_zad2_2.png")

println("Wykres został pomyślnie zapisany do pliku 'wykres_funkcji_zad2.png'")