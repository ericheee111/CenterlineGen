// PolynomialRegressionTest.cpp

#ifdef TEST_POLYNOMIAL_REGRESSION

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include "PolynomialRegression.h"

const int IMG_WIDTH = 512;
const int IMG_HEIGHT = 512;

std::pair<std::vector<double>, std::vector<double>> generatePolynomialData(const std::vector<double>& coeffs, size_t numPoints, double noiseStdDev = 1.0) {
    std::vector<double> x(numPoints);
    std::vector<double> y(numPoints);
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 10.0);
    std::normal_distribution<double> noise(0.0, noiseStdDev);

    for (size_t i = 0; i < numPoints; ++i) {
        x[i] = distribution(generator);
        y[i] = 0.0;
        for (size_t j = 0; j < coeffs.size(); ++j) {
            y[i] += coeffs[j] * std::pow(x[i], j);
        }
        y[i] += noise(generator);
    }
    return std::make_pair( x, y );
}

void plotDataAndFit(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& coeffs, const char* filename) {
    std::vector<unsigned char> image(IMG_WIDTH * IMG_HEIGHT * 3, 0);

    auto getY = [&](double xVal) {
        double result = 0;
        for (size_t i = 0; i < coeffs.size(); ++i) {
            result += coeffs[i] * std::pow(xVal, i);
        }
        return result;
        };

    // Plot original data points (blue)
    for (size_t i = 0; i < x.size(); ++i) {
        int imgX = static_cast<int>((x[i] / *std::max_element(x.begin(), x.end())) * IMG_WIDTH);
        int imgY = static_cast<int>((1.0 - (y[i] / *std::max_element(y.begin(), y.end()))) * IMG_HEIGHT);
        if (imgX >= 0 && imgX < IMG_WIDTH && imgY >= 0 && imgY < IMG_HEIGHT) {
            image[3 * (imgY * IMG_WIDTH + imgX) + 0] = 255;
            image[3 * (imgY * IMG_WIDTH + imgX) + 1] = 255;
            image[3 * (imgY * IMG_WIDTH + imgX) + 2] = 255;
        }
    }

    // Plot fitted curve (red)
    for (int imgX = 0; imgX < IMG_WIDTH; ++imgX) {
        double xVal = (imgX / static_cast<double>(IMG_WIDTH)) * *std::max_element(x.begin(), x.end());
        double yVal = getY(xVal);
        int imgY = static_cast<int>((1.0 - (yVal / *std::max_element(y.begin(), y.end()))) * IMG_HEIGHT);
        if (imgY >= 0 && imgY < IMG_HEIGHT) {
            image[3 * (imgY * IMG_WIDTH + imgX) + 0] = 255;
            image[3 * (imgY * IMG_WIDTH + imgX) + 1] = 0;
            image[3 * (imgY * IMG_WIDTH + imgX) + 2] = 0;
        }
    }

    stbi_write_png(filename, IMG_WIDTH, IMG_HEIGHT, 3, image.data(), IMG_WIDTH * 3);
}

bool areVectorsEqual(const std::vector<double>& a, const std::vector<double>& b, double epsilon = 1e-6) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (std::abs(a[i] - b[i]) > epsilon) return false;
    }
    return true;
}

void testPolynomialRegression(const std::vector<double>& coeffs, int order, const char* plotFilename) {
    auto pair = generatePolynomialData(coeffs, 100, 5.0); // Generate 100 data points with noise
    auto x = pair.first;
    auto y = pair.second;
    PolynomialRegression<double> poly;
    std::vector<double> fittedCoeffs;

    if (poly.fitIt(x, y, order, fittedCoeffs)) {
        std::cout << "Fitted coefficients for order " << order << ": ";
        for (const auto& coeff : fittedCoeffs) {
            std::cout << coeff << " ";
        }
        std::cout << std::endl;

        plotDataAndFit(x, y, fittedCoeffs, plotFilename);
    }
    else {
        std::cout << "Polynomial fitting failed for order " << order << "." << std::endl;
    }
}

int main() {
    // Quadratic polynomial y = 3x^2 + 2x + 1
    std::vector<double> quadraticCoeffs = { 1, 2, 3 };
    testPolynomialRegression(quadraticCoeffs, 2, "quadratic.png");

    // Cubic polynomial y = x^3 - 6x^2 + 11x - 6
    std::vector<double> cubicCoeffs = { -6, 11, -6, 1 };
    testPolynomialRegression(cubicCoeffs, 3, "cubic.png");

    // Quartic polynomial y = 4x^4 - 3x^3 + 2x^2 - x + 5
    std::vector<double> quarticCoeffs = { 5, -1, 2, -3, 4 };
    testPolynomialRegression(quarticCoeffs, 4, "quartic.png");

    return 0;
}

#endif