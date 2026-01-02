# Tomasz Niedziałek 279754
using Test

include("mbisekcji.jl")
include("mstycznych.jl")
include("msiecznych.jl")

@testset "Testy Funkcji mbisekcji" begin
    
    # Test 1: Standardowa zbieżność
    @testset "Test_1_Zbieznosc_sqrt2" begin
        f = x -> x^2 - 2.0
        (r, v, it, err) = mbisekcji(f, 1.0, 2.0, 1e-7, 1e-6)
        @test err == 0
        @test isapprox(r, sqrt(2.0); atol=1e-6) # Sprawdzenie bliskości r
        @test abs(v) < 1e-6                    # Sprawdzenie bliskości f(r) do zera
    end

    # Test 2: Brak zmiany znaku (Błąd)
    @testset "Test_2_Brak_Zmiany_Znaku" begin
        f = x -> x^2 + 1.0
        (r, v, it, err) = mbisekcji(f, -1.0, 1.0, 1e-6, 1e-6)
        @test err == 1
        @test it == 0
    end
    
    # Test 3: Dominacja kryterium Delta
    @testset "Test_3_Kryterium_Delta" begin
        f = x -> x^3 - 4.0
        (r, v, it, err) = mbisekcji(f, 1.0, 2.0, 1e-3, 1e-10) # 1e-3
        # Oczekiwana liczba iteracji dla długości przedziału 1 i delty 1e-3 to ~10
        @test err == 0
        @test it <= 11 # Sprawdź, czy nie wykonano za dużo iteracji
    end

    # Test 4: Złoty Podział (równanie kwadratowe x^2 - x - 1 = 0)
    # Pierwiastek dokładny: (1 + sqrt(5)) / 2
    @testset "Test_4_Zloty_Podzial" begin
        phi_exact = (1 + sqrt(5)) / 2
        f = x -> x^2 - x - 1.0
        
        # Ustawiamy deltę o rząd wielkości mniejszą niż epsilon, 
        # ponieważ pochodna w tym punkcie wynosi ok. 2.23, więc błąd f(x) jest > błędu x.
        (r, v, it, err) = mbisekcji(f, 1.0, 2.0, 1e-7, 1e-6)
        
        @test err == 0
        @test isapprox(r, phi_exact; atol=1e-6)
        @test abs(v) < 1e-6
    end

    # Test 5: Znajdowanie liczby Pi (funkcja sinus)
    # sin(x) = 0 dla x = pi
    @testset "Test_5_Liczba_Pi" begin
        f = x -> sin(x)
        
        # Szukamy w przedziale [3.0, 4.0]. 
        # Pochodna cos(pi) = -1, więc błędy x i f(x) są do siebie zbliżone (stosunek 1:1).
        (r, v, it, err) = mbisekcji(f, 3.0, 4.0, 1e-6, 1e-6)
        
        @test err == 0
        @test isapprox(r, pi; atol=1e-6)
        @test abs(v) < 1e-6
    end
end

@testset "Testy Metody Newtona (Stycznych)" begin

    # Test 1: Pierwiastek sqrt(2) dla f(x) = x^2 - 2
    @testset "Test_1_Standard_sqrt2" begin
        f  = x -> x^2 - 2.0
        pf = x -> 2.0 * x      # Pochodna: 2x
        
        # Startujemy z 2.0. Metoda powinna zbiec bardzo szybko (3-4 iteracje).
        (r, v, it, err) = mstycznych(f, pf, 2.0, 1e-7, 1e-7, 20)
        
        @test err == 0                  # Sukces
        @test isapprox(r, sqrt(2.0); atol=1e-7)
        @test abs(v) < 1e-7
        @test it < 10                   # Newton jest szybki
    end

    # Test 2: Złoty Podział (phi) dla f(x) = x^2 - x - 1
    @testset "Test_2_Zloty_Podzial" begin
        phi_exact = (1 + sqrt(5)) / 2
        f  = x -> x^2 - x - 1.0
        pf = x -> 2.0*x - 1.0  # Pochodna: 2x - 1
        
        (r, v, it, err) = mstycznych(f, pf, 1.0, 1e-7, 1e-7, 20)
        
        @test err == 0
        @test isapprox(r, phi_exact; atol=1e-7)
    end

    # Test 3: Funkcja trygonometryczna (Liczba Pi)
    # f(x) = sin(x), pf(x) = cos(x) -> szukamy zera w okolicy 3.0 (czyli Pi)
    @testset "Test_3_Liczba_Pi" begin
        f  = x -> sin(x)
        pf = x -> cos(x)
        
        # Start z 3.0, blisko Pi (3.1415...)
        (r, v, it, err) = mstycznych(f, pf, 3.0, 1e-8, 1e-8, 20)
        
        @test err == 0
        @test isapprox(r, pi; atol=1e-8)
        @test abs(v) < 1e-8
    end

    # Test 4: Obsługa błędu - Pochodna równa zero (err = 2)
    # f(x) = x^2 - 2. W punkcie x0 = 0 pochodna wynosi 0.
    @testset "Test_4_Pochodna_Zero" begin
        f  = x -> x^2 - 2.0
        pf = x -> 2.0 * x
        
        (r, v, it, err) = mstycznych(f, pf, 0.0, 1e-7, 1e-7, 20)
        
        @test err == 2
        @test it == 0  # Powinno przerwać natychmiast
    end

    # Test 5: Obsługa błędu - Przekroczenie limitu iteracji (err = 1)
    @testset "Test_5_Max_Iteracji" begin
        f  = x -> x^2 - 2.0
        pf = x -> 2.0 * x
        
        # Start z 100.0, tylko 1 iteracja. Nie ma szans dojść do 1.41
        (r, v, it, err) = mstycznych(f, pf, 100.0, 1e-10, 1e-10, 1)
        
        @test err == 1
        # Sprawdźmy czy funkcja zwraca ostatnie przybliżenie, a nie zero
        @test r != 0.0 
    end
    
    # Test 6: Trafienie idealne na start (x0 jest pierwiastkiem)
    @testset "Test_6_Start_W_Pierwiastku" begin
        f  = x -> x^2 - 4.0
        pf = x -> 2.0 * x
        
        # Startujemy idealnie w 2.0
        (r, v, it, err) = mstycznych(f, pf, 2.0, 1e-7, 1e-7, 20)
        
        @test err == 0
        @test it == 0  # 0 iteracji, bo warunek początkowy spełniony
        @test r == 2.0
    end

end

@testset "Testy Metody Siecznych" begin

    # Test 1: Standardowa zbieżność - sqrt(2)
    # Równanie: x^2 - 2 = 0
    @testset "Test_1_Standard_sqrt2" begin
        f = x -> x^2 - 2.0
        # Punkty startowe otaczające pierwiastek (choć w siecznych nie muszą)
        x0, x1 = 1.0, 2.0 
        
        (r, v, it, err) = msiecznych(f, x0, x1, 1e-7, 1e-7, 20)
        
        @test err == 0
        @test isapprox(r, sqrt(2.0); atol=1e-7)
        @test abs(v) < 1e-7
        # Metoda siecznych jest zazwyczaj szybsza od bisekcji, ale wolniejsza od Newtona
        @test it > 0 && it < 15 
    end

    # Test 2: Złoty Podział - (1 + sqrt(5))/2
    # Równanie: x^2 - x - 1 = 0
    @testset "Test_2_Zloty_Podzial" begin
        phi_exact = (1 + sqrt(5)) / 2
        f = x -> x^2 - x - 1.0
        
        # Startujemy z punktów dodatnich
        (r, v, it, err) = msiecznych(f, 1.0, 2.0, 1e-7, 1e-7, 20)
        
        @test err == 0
        @test isapprox(r, phi_exact; atol=1e-7)
    end

    # Test 3: Funkcja trygonometryczna - Liczba Pi
    # Równanie: sin(x) = 0 w okolicy 3.
    @testset "Test_3_Liczba_Pi" begin
        f = x -> sin(x)
        # Start blisko 3.0 i 4.0
        (r, v, it, err) = msiecznych(f, 3.0, 4.0, 1e-8, 1e-8, 20)
        
        @test err == 0
        @test isapprox(r, pi; atol=1e-8)
        @test abs(v) < 1e-8
    end
     
    # Test 4: Ujemny pierwiastek
    # Równanie: x^3 + 1 = 0  ->  x = -1
    @testset "Test_4_Ujemny_Pierwiastek" begin
        f = x -> x^3 + 1.0
        # Startujemy na minusach
        x0, x1 = -2.0, -0.5
        
        (r, v, it, err) = msiecznych(f, x0, x1, 1e-7, 1e-7, 20)
        
        @test err == 0
        @test isapprox(r, -1.0; atol=1e-7)
        @test abs(v) < 1e-7
    end

     
    # Test 5: Przekroczenie limitu iteracji (err = 1)
    @testset "Test_5_Max_Iteracji" begin
        f = x -> x^2 - 2.0
        # Startujemy daleko
        x0, x1 = 10.0, 20.0
        maxit_low = 2
        
        (r, v, it, err) = msiecznych(f, x0, x1, 1e-15, 1e-15, maxit_low)
        
        @test err == 1
        @test it == maxit_low
        # Sprawdzamy czy zwrócono jakąś wartość (nie zero)
        @test r != 0.0 
    end

    # Test 6: Jeden z punktów startowych jest już pierwiastkiem
    @testset "Test_6_Start_W_Zerze" begin
        f = x -> x^2 - 4.0 # Pierwiastki: 2 i -2
        x0 = 2.0 # Idealne trafienie
        x1 = 3.0
        
        # Twoja funkcja sprawdza warunek stopu PO obliczeniu nowego kroku,
        # ale dzięki zamianie zmiennych (swap) 'a' będzie równe 2.0.
        # W pierwszej iteracji fa=0, więc powinno szybko skończyć.
        
        (r, v, it, err) = msiecznych(f, x0, x1, 1e-7, 1e-7, 20)
        
        @test err == 0
        @test r == 2.0
        @test abs(v) == 0.0
    end

end