
include("functions.jl")

for n in [5, 10, 15]
    rysujNnfx(x->exp(x), 0.0, 1.0, n)
    title!("Interpolacja funkcji exp(x) dla n=$n węzłów równoodległych")
    xlabel!("x")
    ylabel!("Nn,f(x) oraz f(x)")
    savefig("wykres_exp_rownoodlegle_n$(n).png")
    println("Zapisano wykres wykres_exp_rownoodlegle_n$(n).png")


    rysujNnfx(x->sin(x)*x^2, -1.0, 1.0, n)
    title!("Interpolacja sin(x)*x^2 dla n=$n węzłów równoodległych")
    xlabel!("x")
    ylabel!("Nn,f(x) oraz f(x)")
    savefig("wykres_sinx2_rownoodlegle_n$(n).png")
    println("Zapisano wykres wykres_sinx2_rownoodlegle_n$(n).png")
end

for n in [5, 10, 15]
    # Węzły równoodległe
    rysujNnfx(x->abs(x), -1.0, 1.0, n)
    title!("Interpolacja funkcji abs(x) dla n=$n węzłów równoodległych")
    xlabel!("x")
    ylabel!("Nn,f(x) oraz f(x)")
    savefig("wykres_abs_rownoodlegle_n$(n).png")
    println("Zapisano wykres wykres_abs_rownoodlegle_n$(n).png")

    rysujNnfx(x->1/(1+x^2), -5.0, 5.0, n)
    title!("Interpolacja funkcji 1/(1+x^2) dla n=$n węzłów równoodległych")
    xlabel!("x")
    ylabel!("Nn,f(x) oraz f(x)")
    savefig("wykres_1_1plusx2_rownoodlegle_n$(n).png")
    println("Zapisano wykres wykres_1_1plusx2_rownoodlegle_n$(n).png")

    # Węzły Czebyszewa
    rysujNnfx(x->abs(x), -1.0, 1.0, n; nodes=:czebyszew)
    title!("Interpolacja funkcji abs(x) dla n=$n Czebyszewa")
    xlabel!("x")
    ylabel!("Nn,f(x) oraz f(x)")
    savefig("wykres_abs_czebyszew_n$(n).png")
    println("Zapisano wykres wykres_abs_czebyszew_n$(n).png")

    rysujNnfx(x->1/(1+x^2), -5.0, 5.0, n; nodes=:czebyszew)
    title!("Interpolacja funkcji 1/(1+x^2) dla n=$n Czebyszewa")
    xlabel!("x")
    ylabel!("Nn,f(x) oraz f(x)")
    savefig("wykres_1_1plusx2_czebyszew_n$(n).png")
    println("Zapisano wykres wykres_1_1plusx2_czebyszew_n$(n).png")
end