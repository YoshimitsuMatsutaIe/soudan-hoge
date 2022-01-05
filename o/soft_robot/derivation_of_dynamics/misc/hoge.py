import time
import concurrent.futures

import numpy as np

M = 10000000000000
N = 1000

def func1():
    for i in range(M):
        print("func1")
        A = np.random.rand(N, N)
        B = np.random.rand(N, N)
        C = A @ B


def func2():
    for i in range(M):
        print("func2")
        A = np.random.rand(10*N, 10*N)
        B = np.random.rand(10*N, 10*N)
        C = A @ B


if __name__ == "__main__":
    executor = concurrent.futures.ProcessPoolExecutor(max_workers=2)
    executor.submit(func1)
    executor.submit(func2)