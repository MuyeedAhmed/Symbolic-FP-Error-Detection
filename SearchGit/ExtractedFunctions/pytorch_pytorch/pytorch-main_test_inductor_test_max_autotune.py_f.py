        def f(a, b):
            # Prologue: pointwise operation on input 'a' before matmul
            a_transformed = a + 1.0
            # Matmul
            mm_result = a_transformed @ b
            # Epilogue: pointwise operation on output after matmul
            return mm_result + 2.0