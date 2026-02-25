
#include <iostream>
#include <cmath>
#include <vector>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <limits>
#include <string>
#include <algorithm>
#include <tuple>

// Define the function to integrate
double customFunc(double val) {
    if (val == 0) {
        return 0;
    }
    return val * val * exp(val) * sin(1 / val);
}

// Trapezoidal Rule Implementation
std::tuple<double, double, std::vector<double>, std::vector<double>>
calculateAreaUsingTrapezoidal(double (*func)(double), double lower, double upper, int intervals) {
    auto start = std::chrono::high_resolution_clock::now();

    if (lower == 0) lower = 1e-10;

    double stepSize = (upper - lower) / intervals;
    std::vector<double> xs(intervals + 1);
    std::vector<double> ys(intervals + 1);

    for (int i = 0; i <= intervals; i++) {
        xs[i] = lower + i * stepSize;
        ys[i] = func(xs[i]);
    }

    double totalArea = 0.5 * ys[0];
    for (int i = 1; i < intervals; i++) {
        totalArea += ys[i];
    }
    totalArea += 0.5 * ys[intervals];
    totalArea *= stepSize;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> timeElapsed = end - start;

    return std::make_tuple(totalArea, timeElapsed.count(), xs, ys);
}

// Simpson's Rule Implementation
std::tuple<double, double, std::vector<double>, std::vector<double>>
calculateAreaUsingSimpsons(double (*func)(double), double lower, double upper, int intervals) {
    auto start = std::chrono::high_resolution_clock::now();

    if (lower == 0) lower = 1e-10;
    if (intervals % 2 != 0) intervals++;

    double stepSize = (upper - lower) / intervals;
    std::vector<double> xs(intervals + 1);
    std::vector<double> ys(intervals + 1);

    for (int i = 0; i <= intervals; i++) {
        xs[i] = lower + i * stepSize;
        ys[i] = func(xs[i]);
    }

    double totalArea = ys[0];
    for (int i = 1; i < intervals; i++) {
        totalArea += (i % 2 == 0) ? 2 * ys[i] : 4 * ys[i];
    }
    totalArea += ys[intervals];
    totalArea *= stepSize / 3;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> timeElapsed = end - start;

    return std::make_tuple(totalArea, timeElapsed.count(), xs, ys);
}

// Generate plot using gnuplot
void createPlot(double (*func)(double), double lower, double upper, int intervals,
    const std::string& methodName, double result,
    const std::vector<double>& xs, const std::vector<double>& ys) {
    std::ofstream dataFile("data.txt");
    for (size_t i = 0; i < xs.size(); i++) {
        dataFile << xs[i] << " " << ys[i] << std::endl;
    }
    dataFile.close();

    std::ofstream smoothDataFile("smooth_data.txt");
    int points = 1000;
    double increment = (upper - lower) / (points - 1);

    for (int i = 0; i < points; i++) {
        double xVal = lower + i * increment;
        if (std::abs(xVal) < 1e-10) xVal = 1e-10;
        smoothDataFile << xVal << " " << func(xVal) << std::endl;
    }
    smoothDataFile.close();

    std::ofstream gnuplotScript("plot_data.gp");
    gnuplotScript << "set terminal png size 1000,600\n"
        << "set output '" << methodName << "_plot.png'\n"
        << "set title '" << methodName << " Area = " << std::fixed << std::setprecision(10) << result << ", intervals = " << intervals << "'\n"
        << "set xlabel 'x'\nset ylabel 'f(x)'\nset grid\n"
        << "set style fill transparent solid 0.5 noborder\n"
        << "plot 'smooth_data.txt' with lines title 'f(x) = x^2e^x sin(1/x)' lw 2, \\\n"
        << "     'data.txt' with points pt 7 ps 0.5 lc rgb 'red' title 'Sample Points', \\\n"
        << "     'data.txt' using 1:2:(0) with filledcurves title 'Area'\n";
    gnuplotScript.close();

    system("gnuplot plot_data.gp");
    std::cout << "Plot saved as '" << methodName << "_plot.png'\n";
}

// Error Estimation
template<typename T>
std::tuple<std::vector<T>, std::vector<T>>
calculateError(double (*func)(double), double lower, double upper, const std::vector<int>& intervalsList, const std::string& method) {
    std::vector<T> results;

    for (int intervals : intervalsList) {
        double area, time;
        std::vector<double> x, y;
        if (method == "trapezoidal")
            std::tie(area, time, x, y) = calculateAreaUsingTrapezoidal(func, lower, upper, intervals);
        else
            std::tie(area, time, x, y) = calculateAreaUsingSimpsons(func, lower, upper, intervals);

        results.push_back(area);
    }

    std::vector<T> errors;
    for (size_t i = 0; i < results.size() - 1; ++i) {
        errors.push_back(std::abs(results[i + 1] - results[i]));
    }

    return { results, errors };
}

int main() {
    std::cout << "Numerical Integration Calculator\n";
    std::cout << "Function: ∫₀¹ x²e^x sin(1/x) dx\n\n";

    try {
        double lowerBound = 0.0, upperBound = 1.0;
        int numIntervals = 100;
        std::string methodChoice = "3";

        std::string input;
        std::cout << "Enter lower bound [default 0]: ";
        std::getline(std::cin, input);
        if (!input.empty()) lowerBound = std::stod(input);

        std::cout << "Enter upper bound [default 1]: ";
        std::getline(std::cin, input);
        if (!input.empty()) upperBound = std::stod(input);

        std::cout << "Enter number of intervals [default 100]: ";
        std::getline(std::cin, input);
        if (!input.empty()) numIntervals = std::stoi(input);

        std::cout << "Choose method: 1-Trapezoidal, 2-Simpson's, 3-Both [default 3]: ";
        std::getline(std::cin, input);
        if (!input.empty()) methodChoice = input;

        double trapArea = 0, trapTime = 0, simpsonArea = 0, simpsonTime = 0;
        std::vector<double> trapX, trapY, simpsonX, simpsonY;

        if (methodChoice == "1" || methodChoice == "3") {
            std::tie(trapArea, trapTime, trapX, trapY) = calculateAreaUsingTrapezoidal(customFunc, lowerBound, upperBound, numIntervals);
            std::cout << "\n=== Trapezoidal Rule ===\n";
            std::cout << "  Area       : " << std::fixed << std::setprecision(10) << trapArea << "\n";
            std::cout << "  Time Taken : " << std::fixed << std::setprecision(6) << trapTime << " seconds\n";
            createPlot(customFunc, lowerBound, upperBound, numIntervals, "Trapezoidal_Rule", trapArea, trapX, trapY);
        }

        if (methodChoice == "2" || methodChoice == "3") {
            std::tie(simpsonArea, simpsonTime, simpsonX, simpsonY) = calculateAreaUsingSimpsons(customFunc, lowerBound, upperBound, numIntervals);
            std::cout << "\n=== Simpson's Rule ===\n";
            std::cout << "  Area       : " << std::fixed << std::setprecision(10) << simpsonArea << "\n";
            std::cout << "  Time Taken : " << std::fixed << std::setprecision(6) << simpsonTime << " seconds\n";
            createPlot(customFunc, lowerBound, upperBound, numIntervals, "Simpsons_Rule", simpsonArea, simpsonX, simpsonY);
        }

        std::vector<int> intervalsList = { numIntervals, 2 * numIntervals, 4 * numIntervals, 8 * numIntervals };

        std::cout << "\n=== Error Estimation ===\n";

        if (methodChoice == "1" || methodChoice == "3") {
            auto [trapResults, trapErrors] = calculateError<double>(customFunc, lowerBound, upperBound, intervalsList, "trapezoidal");
            std::cout << "\n--- Trapezoidal Rule Convergence ---\n";
            for (size_t i = 0; i < trapResults.size(); ++i)
                std::cout << "  Intervals = " << intervalsList[i] << ": Result = " << std::fixed << std::setprecision(10) << trapResults[i] << "\n";
            for (size_t i = 0; i < trapErrors.size(); ++i)
                std::cout << "  Error between " << intervalsList[i] << " and " << intervalsList[i + 1] << " = "
                << std::scientific << trapErrors[i] << "\n";
        }

        if (methodChoice == "2" || methodChoice == "3") {
            auto [simpResults, simpErrors] = calculateError<double>(customFunc, lowerBound, upperBound, intervalsList, "simpson");
            std::cout << "\n--- Simpson's Rule Convergence ---\n";
            for (size_t i = 0; i < simpResults.size(); ++i)
                std::cout << "  Intervals = " << intervalsList[i] << ": Result = " << std::fixed << std::setprecision(10) << simpResults[i] << "\n";
            for (size_t i = 0; i < simpErrors.size(); ++i)
                std::cout << "  Error between " << intervalsList[i] << " and " << intervalsList[i + 1] << " = "
                << std::scientific << simpErrors[i] << "\n";
        }

        if (methodChoice == "3") {
            std::cout << "\n=== Method Comparison ===\n";
            std::cout << "  Trapezoidal Rule : Area = " << std::fixed << std::setprecision(10) << trapArea << ", Time = " << trapTime << "s\n";
            std::cout << "  Simpson's Rule   : Area = " << simpsonArea << ", Time = " << simpsonTime << "s\n";
            double diff = std::abs(trapArea - simpsonArea);
            std::cout << "  Difference        : " << std::scientific << diff << "\n";
            if (simpsonTime < trapTime)
                std::cout << "  Simpson's Rule is " << std::fixed << std::setprecision(2) << trapTime / simpsonTime << "x faster\n";
            else
                std::cout << "  Trapezoidal Rule is " << std::fixed << std::setprecision(2) << simpsonTime / trapTime << "x faster\n";
        }

    }
    catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
