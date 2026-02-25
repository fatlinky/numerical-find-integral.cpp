import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

def f(x):
    with np.errstate(divide='ignore', invalid='ignore'):
        result = x**2 * np.exp(x) * np.sin(1/x)
        if isinstance(result, np.ndarray):
            result[np.isinf(result) | np.isnan(result)] = np.nan
        elif np.isnan(result) or np.isinf(result):
            result = np.nan
    return result

def trapezoidal_rule(a, b, n):
    h = (b - a) / n
    x = np.linspace(a, b, n+1)
    y = f(x)
    if np.any(np.isnan(y)):
        raise ValueError("Function is undefined at some point in the interval.")
    area = (h/2) * (y[0] + 2 * np.sum(y[1:-1]) + y[-1])
    return area

def simpsons_rule(a, b, n):
    if n % 2 == 1:
        n += 1  # Simpson's rule needs even number of intervals
    h = (b - a) / n
    x = np.linspace(a, b, n+1)
    y = f(x)
    if np.any(np.isnan(y)):
        raise ValueError("Function is undefined at some point in the interval.")
    area = (h/3) * (y[0] + 2*np.sum(y[2:n:2]) + 4*np.sum(y[1:n:2]) + y[n])
    return area

def plot_function(a, b, n, method):
    x_dense = np.linspace(a, b, 1000)
    y_dense = f(x_dense)

    plt.figure(figsize=(10, 6))
    plt.plot(x_dense, y_dense, 'b', label='f(x)')
    plt.fill_between(x_dense, 0, y_dense, where=~np.isnan(y_dense), color='skyblue', alpha=0.4)

    # Subdivision points
    x = np.linspace(a, b, n + 1)
    y = f(x)

    # Draw vertical lines at subdivision points
    for xi in x:
        plt.axvline(x=xi, color='gray', linestyle='--', linewidth=0.6)

    # For trapezoidal rule, draw trapezoids
    if method == 'trapezoidal':
        for i in range(n):
            xs = [x[i], x[i], x[i+1], x[i+1]]
            ys = [0, y[i], y[i+1], 0]
            plt.fill(xs, ys, color='green', alpha=0.4, edgecolor='orangered', linewidth=1)

    plt.title(f'Function & Area Under the Curve ({method.title()} Rule)')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.legend()
    plt.grid(True)
    plt.show()


def main():
    a = float(input("Enter lower limit a: "))
    b = float(input("Enter upper limit b: "))
    n = int(input("Enter number of subdivisions n: "))
    method = input("Choose method (trapezoidal/simpson): ").strip().lower()

    if a <= 0 <= b:
        print("Error: Function is undefined at x = 0. Please choose an interval that does not include x = 0.")
        return

    try:
        if method == 'trapezoidal':
            area = trapezoidal_rule(a, b, n)
        elif method == 'simpson':
            area = simpsons_rule(a, b, n)
        else:
            print("Invalid method selected.")
            return

        exact_area, _ = quad(f, a, b, limit=100, points=[0])  # tell quad 0 is a known discontinuity

        error = abs(exact_area - area)

        print(f"\nApproximate Area ({method.title()} Rule) = {area:.6f}")
        print(f"Exact Area = {exact_area:.6f}")
        print(f"Absolute Error = {error:.6f}")
        print(f"Relative Error = {error/exact_area:.6%}")

        plot_function(a, b, n,method)

    except ValueError as e:
        print(f"Error: {e}")

if __name__ == "__main__":

    main()