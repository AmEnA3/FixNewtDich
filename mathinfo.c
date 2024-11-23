#include <stdio.h>
#include <math.h>

// Define the function and its derivative (for Newton's method)
double f(double x) {
    return x*x - 2;  // Example: f(x) = x^2 - 2 (root is sqrt(2))
}

double f_prime(double x) {
    return 2*x;  // Derivative of f(x) = x^2 - 2 is f'(x) = 2x
}

// Fixed-Point Method
double fixed_point(double (*g)(double), double x0, double tol, int max_iter) {
    double x = x0;
    for (int i = 0; i < max_iter; i++) {
        double x_new = g(x);
        if (fabs(x_new - x) < tol) {
            return x_new;
        }
        x = x_new;
    }
    return x; // Return the result after max_iter if tolerance is not met
}

// Newton's Method
double newton(double (*f)(double), double (*f_prime)(double), double x0, double tol, int max_iter) {
    double x = x0;
    for (int i = 0; i < max_iter; i++) {
        double fx = f(x);
        double fpx = f_prime(x);
        if (fabs(fx) < tol) {
            return x;
        }
        x = x - fx / fpx;
    }
    return x; // Return after max_iter if tolerance is not met
}

// Bisection Method
double bisection(double (*f)(double), double a, double b, double tol, int max_iter) {
    if (f(a) * f(b) > 0) {
        printf("f(a) and f(b) must have opposite signs!\n");
        return -1;  // Error: No root found in the interval
    }

    double c;
    for (int i = 0; i < max_iter; i++) {
        c = (a + b) / 2;
        if (fabs(f(c)) < tol) {
            return c;  // Root found
        }
        if (f(c) * f(a) < 0) {
            b = c;  // Root is in the left half
        } else {
            a = c;  // Root is in the right half
        }
    }
    return c;  // Return after max_iter if tolerance is not met
}

// Example of g(x) for Fixed-Point Method (rearranged f(x) = 0)
double g(double x) {
    return sqrt(2);  // Example for fixed-point method
}

int main() {
    double tol = 1e-6;  // Tolerance
    int max_iter = 100; // Maximum number of iterations
    double x0 = 1.0;    // Initial guess for Newton's and Fixed-Point methods
    double a = 0.0, b = 2.0;  // Initial interval for Bisection Method

    // Fixed-Point Method
    double root_fp = fixed_point(g, x0, tol, max_iter);
    printf("Fixed-Point Method Root: %f\n", root_fp);

    // Newton's Method
    double root_newton = newton(f, f_prime, x0, tol, max_iter);
    printf("Newton's Method Root: %f\n", root_newton);

    // Bisection Method
    double root_bisect = bisection(f, a, b, tol, max_iter);
    printf("Bisection Method Root: %f\n", root_bisect);

    return 0;
}
