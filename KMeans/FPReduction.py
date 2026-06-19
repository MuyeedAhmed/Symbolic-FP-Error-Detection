import struct
import numpy as np
import z3

class FPReduction:
    '''
    left_fold  : sequential accumulate, ((a+b)+c)+d — CPU BLAS (MKL / OpenBLAS)
    '''
    @staticmethod
    def left_fold(arr: np.ndarray) -> np.float32:
        """Sequential left-to-right accumulate: ((a0+a1)+a2)+…"""
        arr = np.asarray(arr, dtype=np.float32)
        acc = arr[0]
        for v in arr[1:]:
            acc = np.float32(acc + v)
        return acc

    @staticmethod
    def right_fold(arr: np.ndarray) -> np.float32:
        """Sequential right-to-left accumulate: …+(a2+(a1+a0))"""
        arr = np.asarray(arr, dtype=np.float32)
        acc = arr[-1]
        for v in arr[-2::-1]:
            acc = np.float32(acc + v)
        return acc
    
    '''
    tree_fold  : pairwise lane merge, (a+b)+(c+d) - GPU SIMD (cuBLAS / CUTLASS)
    '''
    @staticmethod
    def tree_fold(arr: np.ndarray) -> np.float32:
        """Pairwise tree reduction: 2-wide SIMD lane merge (GPU default)."""
        arr = np.asarray(arr, dtype=np.float32).copy()
        while len(arr) > 1:
            half  = len(arr) // 2
            pairs = np.float32(arr[:2*half:2]) + np.float32(arr[1:2*half:2])
            arr   = np.concatenate([pairs, arr[2*half:]]) if len(arr) % 2 else pairs
        return arr[0]

    '''
    chunk_fold : tiled left-fold, models configurable BLAS tile width
    '''
    @staticmethod
    def chunk_fold(arr: np.ndarray, chunk_size: int) -> np.float32:
        arr      = np.asarray(arr, dtype=np.float32)
        partials = [FPReduction.left_fold(arr[i:i+chunk_size])
                    for i in range(0, len(arr), chunk_size)]
        return FPReduction.left_fold(np.array(partials, dtype=np.float32))

    _REGISTRY = None

    @classmethod
    def get(cls, name: str):
        if cls._REGISTRY is None:
            cls._REGISTRY = {
                "left":  cls.left_fold,
                "right": cls.right_fold,
                "tree":  cls.tree_fold,
            }
            for w in [2, 4, 8, 16, 32, 64]:
                cls._REGISTRY[f"chunk_{w}"] = (lambda w=w: (lambda arr: cls.chunk_fold(arr, w)))()
        if name not in cls._REGISTRY:
            raise ValueError(f"Unknown reduction: '{name}'.  "
                             f"Available: {list(cls._REGISTRY)}")
        return cls._REGISTRY[name]

    @staticmethod
    def z3_left(terms: list, rm) -> object:
        acc = terms[0]
        for t in terms[1:]:
            acc = z3.fpAdd(rm, acc, t)
        return acc

    @staticmethod
    def z3_tree(terms: list, rm) -> object:
        t = list(terms)
        while len(t) > 1:
            pairs = [z3.fpAdd(rm, t[i], t[i+1]) for i in range(0, len(t)-1, 2)]
            if len(t) % 2:
                pairs.append(t[-1])
            t = pairs
        return t[0]

    @classmethod
    def z3_get(cls, name: str):
        """Return a Z3 symbolic fold builder by name ('left' or 'tree')."""
        if name == "left":
            return cls.z3_left
        if name == "tree":
            return cls.z3_tree
        raise ValueError(f"Z3 fold '{name}' not supported (add it to z3_get).")

    @staticmethod
    def f32_to_bits(v: float) -> int:
        return struct.unpack("I", struct.pack("f", np.float32(v)))[0]

    @staticmethod
    def bits_to_f32(n: int) -> np.float32:
        return np.float32(struct.unpack("f", struct.pack("I", n & 0xFFFF_FFFF))[0])
