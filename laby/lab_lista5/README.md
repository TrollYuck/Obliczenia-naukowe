Projekt implementuje efektywne algorytmy rozwiazywania układów liniowych Ax =b, gdzie macierz A jest blokowo-trójdiagonalna.

Struktura:
-   main.jl - Interface terminala.
-   blocksys.jl - Definicja modułu głownego.
-   blockmat.jl - Moduł generatora macierzy testowych profesora Zielińskiego.
-   structures.jl - Definicje struktur danych: BlockMatrix (A), BorderMatrix (B).
-   zad1.jl - Implementacja metody Eliminacji Gaussa w wersjach z pivotem oraz bez.
-   zad2.jl - Implementacja rozkładu LU macierzy metodą Eliminacji Gaussa w wersjach z pivotem oraz bez.
-   zad3.jl - Implementacja funkcji rozwiązującej układ rówań Ax = b, z uwzględnieniem rozkładu LU z zad2.jl
-   tests.jl - Testy poprawności implementacji, wyniki porównywane są z implementacją modułu LinearAlgebra
-   manage_input - Implementacja funkcji odczytu oraz zapisu do plików

Użycie:
    julia main.jl [opcje]
Dostępne opcje

Wejście danych:

    -A <plik> : Ścieżka do pliku z macierzą A,

    -b <plik> : Ścieżka do pliku z wektorem b.

    -G : Wygeneruj losową macierz A (wymaga podania -n i -l).

    -g : Wygeneruj wektor b na podstawie znanego rozwiązania x=(1,…,1)^T.

Parametry generatora (dla flagi -G):

    -n <int> : Rozmiar macierzy.

    -l <int> : Rozmiar pojedynczego bloku.

    -ck <float> : Współczynnik uwarunkowania macierzy (domyślnie 10.0).

Metody rozwiązywania (flaga -m):

    gauss_no_pivot : Eliminacja Gaussa bez wyboru elementu głównego.

    gauss_pivot : Eliminacja Gaussa z częściowym wyborem elementu głównego.

    lu_no_pivot : Rozkład LU + rozwiązanie (bez wyboru elementu głównego).

    lu_pivot : Rozkład LU + rozwiązanie (z częściowym wyborem elementu głównego).

Tryby testowe:

    -test-gauss : Uruchamia obie metody Gaussa (z pivotem i bez) i wyświetla tabelę porównawczą czasu i błędów.

    -test-lu : Uruchamia obie metody LU (z pivotem i bez) i wyświetla tabelę porównawczą.

Wyjście:

    -o <plik> : Ścieżka do pliku wynikowego