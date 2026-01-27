# Coppersmith RSA Attack

Projekt realizuje atak Coppersmitha na system RSA w przypadku małego wykładnika publicznego e. Implementacja koncentruje się na scenariuszu Stereotyped Message Attack, wykorzystując redukcję baz krat algorytmem LLL do odzyskania brakującej części wiadomości.

## Kluczowe cechy
* Własna implementacja algorytmów: Brak użycia gotowych bibliotek kryptograficznych.
* Stabilność numeryczna: Wykorzystanie biblioteki gmpy2 z precyzją 2048 bitów do obsługi ogromnych współczynników wielomianów.

## Struktura plików
* `coppersmith.py` – Główny moduł: budowanie macierzy kraty, algorytm LLL i odzyskiwanie pierwiastków.
* `poly.py` – Biblioteka pomocnicza do arytmetyki wielomianów (klasa Poly).
* `run_tests.py` – Skrypt uruchamiający testy z katalogu tests/.
* `tests/` – Katalog z danymi wejściowymi (N, e, C, M0) i oczekiwanymi wynikami (M1).

## Wymagania
Do poprawnego działania wymagany jest Python 3 oraz biblioteka gmpy2.

```bash
pip install gmpy2
```

## Instrukcja uruchomienia

### Weryfikacja poprawności (Testy)
Aby uruchomić zestaw 15 testów sprawdzających skuteczność ataku dla różnych długości kluczy należy uruchomić `run_tests.py`:
```bash
python run_tests.py
```

## Założenia techniczne
* Model ataku: M = M0 + M1, gdzie M0 (prefiks) jest znany.
* Szukanie pierwiastka: Zastosowano metodę iteracyjną w zakresie [-10^7, 10^7] ze względu na stabilność obliczeń na liczbach całkowitych przy bardzo dużych współczynnikach.
* Heurystyka: Algorytm analizuje 10 najkrótszych wektorów po redukcji bazy kraty.

