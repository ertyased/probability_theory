from decimal import *
from math import e, pi
from scipy.integrate import quad
getcontext().prec = 10
factorials = [Decimal(1)] * 10001
for i in range(1, 10001):
    factorials[i] = factorials[i - 1] * Decimal(i)


def collect_data(arr, deviation):
    most_probable = 0
    sum_probs = 0
    for i in range(len(arr)):
        if abs(len(arr) / 2 - i) < deviation:
            sum_probs += arr[i]
        if arr[i] > arr[most_probable]:
            most_probable = i
    return [sum_probs, most_probable, arr[most_probable]]


def test_calc(n, p):
    global factorials
    n_non_dec = int(n)
    q = 1 - p
    probabilities = [0] * (n_non_dec + 1)
    probabilities[0] = q**n
    for i in range(1, n_non_dec):
        probabilities[i] = probabilities[i - 1] * ((n - i) / (i + 1)) * p / q
    deviation = (n * p * q) ** 2
    return collect_data(probabilities, deviation)


def test_Puasson(n, p):
    n_non_dec = int(n)
    l = n * p
    const = Decimal(e)**(-l)
    probabilites = [0] * (n_non_dec + 1)
    for i in range(n_non_dec):
        probabilites[i] = const * l ** i / factorials[i]
    deviation = (n * p * (1 - p)) ** 2
    return collect_data(probabilites, deviation)


def test_local_Moivre_Laplace(n, p):
    n_non_dec = int(n)
    q = 1 - p
    probabilites = [0] * (n_non_dec + 1)
    for i in range(n_non_dec):
        x_i = (i - n * p) / (n * p * q)**Decimal(0.5)
        probabilites[i] = 1 / (2 * Decimal(pi) * n * p * q)**Decimal(0.5) * Decimal(e) ** (- x_i ** 2 / 2)

    deviation = (n * p * (1 - p)) ** 2
    return collect_data(probabilites, deviation)


def f(x):
    return e**(-x**2/2)


def phi(x):
    const = 1 / (2 * pi)**(0.5)
    return const * quad(f, 0, x)[0]

def test_integral_Moivre_Laplace(n, p):
    n = int(n)
    p = float(p)
    q = (1 - p)
    deviation = (n * p * q) ** 0.5
    left = n / 2 - deviation
    right = n / 2 + deviation
    x1 = (left - n * p) / deviation
    x2 = (right - n * p) / deviation
    left = phi(x1)
    right = phi(x2)
    return right - left


n = [Decimal(100), Decimal(1000), Decimal(10000)]
p = [Decimal("0.001"), Decimal("0.01"), Decimal("0.1"), Decimal("0.25"), Decimal("0.5")]
for i in n:
    for j in p:
        print(f"Testing for n = {i} and p = {j}:")
        test1 = test_calc(i, j)
        print(f"\tResults for rough calculations:\n\t\tProbability: {test1[0]}\n\t\tMost probable: {test1[1]}\n\t\tMost probable probability:{test1[2]}")
        test1 = test_Puasson(i, j)
        print(f"\tResults for Puasson calculations:\n\t\tProbability: {test1[0]}\n\t\tMost probable: {test1[1]}\n\t\tMost probable probability:{test1[2]}")
        test1 = test_local_Moivre_Laplace(i, j)
        print(f"\tResults for local Moivre-Laplace calculations:\n\t\tProbability: {test1[0]}\n\t\tMost probable: {test1[1]}\n\t\tMost probable probability:{test1[2]}")
        test1 = test_integral_Moivre_Laplace(i, j)
        print(f"\tResults for integral Moivre-Laplace calculations:\n\t\tProbability: {test1}")


