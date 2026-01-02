#Tomasz Niedziałek 279754

using Test

include("functions.jl")

@testset "Testy Laboratorium 4" begin

    @testset "Zadanie 1: Ilorazy Różnicowe" begin
        # Przypadek 1: Funkcja stała f(x) = 3
        # Ilorazy wyższych rzędów powinny być 0
        x = [1.0, 2.0, 3.0]
        f = [3.0, 3.0, 3.0]
        fx = ilorazyRoznicowe(x, f)
        @test fx ≈ [3.0, 0.0, 0.0] atol=1e-10

        # Przypadek 2: Funkcja liniowa f(x) = 2x + 1
        # f(-1)=-1, f(0)=1, f(1)=3
        # Oczekiwany wielomian Newtona: -1 + 2(x - (-1)) + 0(...)
        x = [-1.0, 0.0, 1.0]
        f = [-1.0, 1.0, 3.0]
        fx = ilorazyRoznicowe(x, f)
        @test fx ≈ [-1.0, 2.0, 0.0] atol=1e-10

        # Przypadek 3: Funkcja kwadratowa f(x) = x^2
        # Węzły: 0, 1, 2. Wartości: 0, 1, 4.
        # Ilorazy: 
        # rząd 0: 0, 1, 4
        # rząd 1: (1-0)/1=1, (4-1)/1=3
        # rząd 2: (3-1)/(2-0)=1
        # Wynik: [0, 1, 1] -> 0 + 1(x-0) + 1(x-0)(x-1) = x + x^2 - x = x^2
        x = [0.0, 1.0, 2.0]
        f = [0.0, 1.0, 4.0]
        fx = ilorazyRoznicowe(x, f)
        @test fx ≈ [0.0, 1.0, 1.0] atol=1e-10
    end

    @testset "Zadanie 2: Wartość wielomianu (Newton)" begin
        # Test na podstawie funkcji kwadratowej z poprzedniego testu
        # f(x) = x^2, współczynniki Newtona dla węzłów [0, 1, 2] to [0, 1, 1]
        x = [0.0, 1.0, 2.0]
        fx = [0.0, 1.0, 1.0] # obliczone wcześniej ilorazy

        # Sprawdzenie w węzłach (interpolacja musi się zgadzać)
        @test warNewton(x, fx, 0.0) ≈ 0.0 atol=1e-10
        @test warNewton(x, fx, 1.0) ≈ 1.0 atol=1e-10
        @test warNewton(x, fx, 2.0) ≈ 4.0 atol=1e-10

        # Sprawdzenie poza węzłami
        # f(3) = 3^2 = 9
        @test warNewton(x, fx, 3.0) ≈ 9.0 atol=1e-10
        # f(0.5) = 0.25
        @test warNewton(x, fx, 0.5) ≈ 0.25 atol=1e-10
    end

    @testset "Zadanie 3: Postać naturalna" begin
        # Przypadek 1: f(x) = x^2
        # Newton: 0 + 1(x) + 1(x)(x-1) = x^2
        # Naturalna: 0 + 0x + 1x^2 -> współczynniki [0, 0, 1]
        x = [0.0, 1.0, 2.0]
        fx = [0.0, 1.0, 1.0] 
        a = naturalna(x, fx)
        @test a ≈ [0.0, 0.0, 1.0] atol=1e-10

        # Przypadek 2: f(x) = 2x + 1
        # Newton dla węzłów [-1, 0, 1]: [-1, 2, 0]
        # Naturalna: 1 + 2x + 0x^2 -> współczynniki [1, 2, 0]
        x = [-1.0, 0.0, 1.0]
        fx = [-1.0, 2.0, 0.0]
        a = naturalna(x, fx)
        @test a ≈ [1.0, 2.0, 0.0] atol=1e-10
        
        # Przypadek 3: Wielomian (x-1)(x-2) = x^2 - 3x + 2
        # Węzły: 1, 2, 3. Wartości: 0, 0, 2
        # Ilorazy:
        # 0, 0, 2
        # 0, (2-0)/1=2
        # (2-0)/(3-1)=1
        # fx = [0, 0, 1]
        # Oczekiwane naturalne: 2 - 3x + 1x^2 -> [2, -3, 1]
        x = [1.0, 2.0, 3.0]
        fx = [0.0, 0.0, 1.0]
        a = naturalna(x, fx)
        @test a ≈ [2.0, -3.0, 1.0] atol=1e-10
    end

    @testset "Zadanie 4: Rysowanie" begin
        # Sprawdzamy czy funkcja się wykonuje i zwraca obiekt wykresu
        fi(x) = x^2
        p = rysujNnfx(fi, -1.0, 1.0, 5; nodes=:rownoodlegle)
        @test typeof(p) <: Plots.Plot

        p2 = rysujNnfx(fi, -1.0, 1.0, 5; nodes=:czebyszew)
        @test typeof(p2) <: Plots.Plot

        # Test błędu dla złego typu węzłów
        @test_throws ArgumentError rysujNnfx(fi, -1.0, 1.0, 5; nodes=:nieznany)
    end

end