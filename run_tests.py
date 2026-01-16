import os
import json
import glob
from coppersmith import coppersmith_univariate


def run_tests(directory="tests"):
    test_files = glob.glob(os.path.join(directory, "*.*"))
    if not test_files:
        print(f"Brak plików testowych w katalogu '{directory}'")
        return

    success_count = 0
    total_count = len(test_files)

    print(f"--- ROZPOCZĘCIE TESTÓW (Suma testów: {total_count}) ---")

    for i, file_path in enumerate(test_files, 1):
        file_name = os.path.basename(file_path)
        with open(file_path, 'r') as f:
            try:
                data = json.load(f)
            except json.JSONDecodeError:
                print(f"Test {i}: {file_name} - BŁĄD: Niepoprawny format JSON")
                continue

        N = data.get('N')
        e = data.get('e')
        C = data.get('C')
        M0 = data.get('M0')
        expected_x = data.get('M1')

        bits = N.bit_length()
        if bits <= 128:
            s, t = 2, 4
        elif bits <= 512:
            s, t = 3, 5
        elif bits <= 1024:
            s, t = 4, 6
        else:
            s, t = 5, 7

        found_x = coppersmith_univariate(N, e, C, M0, s=s, t=t, delta=0.99)

        is_correct = (found_x == expected_x) if expected_x is not None else (found_x is not None)
        if is_correct:
            success_count += 1
            status = "OK"
        else:
            status = "FAIL"

        params_str = f"N={str(N)[:10]}..., e={e}, s={s}, t={t}"
        exp_str = expected_x if expected_x is not None else "N/A"
        got_str = found_x if found_x is not None else "None"

        print(f"Test {i}: {file_name}, Parametry: ({params_str}), "
              f"Oczekiwana odpowiedź: {exp_str}, Otrzymana odpowiedź: {got_str} [{status}]")

    print("-" * 50)
    print(f"PODSUMOWANIE: Poprawność: {success_count}/{total_count}")
    if total_count > 0:
        accuracy = (success_count / total_count) * 100
        print(f"Skuteczność: {accuracy:.2f}%")


if __name__ == "__main__":
    run_tests("tests")