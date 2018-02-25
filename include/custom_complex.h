//
// Created by nephtys on 23.01.18.
//

#ifndef SCARLET_STUDY_CUSTOM_COMPLEX_H
#define SCARLET_STUDY_CUSTOM_COMPLEX_H

#include <math.h>

typedef struct complex_number {
    double real;
    double imag;
} complex_number_t;



complex_number_t complexify(double r, double i) {
    return (complex_number_t) {
            .real = r, .imag = i
    };
}

complex_number_t cmul(complex_number_t number, double scalar) {
    return (complex_number_t) {
            .real = number.real * scalar, .imag = number.imag * scalar
    };
}
complex_number_t cdiv(complex_number_t number, double scalar) {
    return (complex_number_t) {
            .real = number.real / scalar, .imag = number.imag / scalar
    };
}
complex_number_t cinverse(complex_number_t number) {
    //for input x: cinverse(x) = 1/x
    const double a = number.real;
    const double b = number.imag;
    const double a2_b2 = (a*a) + (b*b);
    return (complex_number_t) {
            .real = (a / a2_b2) , .imag = ((-1.0) * b / a2_b2)
    };
}
complex_number_t cmulc(complex_number_t a, complex_number_t b) {
    return (complex_number_t) {
            .real = (a.real * b.real) - (a.imag * b.imag), .imag = (a.real * b.imag) + (a.imag * b.real)
    };
}
complex_number_t cadd(complex_number_t a, complex_number_t b) {
    return (complex_number_t) {
            .real = a.real + b.real, .imag = a.imag + b.imag
    };
}
complex_number_t csub(complex_number_t a, complex_number_t b) {
    return (complex_number_t) {
            .real = a.real - b.real, .imag = a.imag - b.imag
    };
}
double complex_abs(complex_number_t number) {
    return sqrt( (number.real * number.real) + (number.imag * number.imag)) ;
}
double complex_real(complex_number_t number) {
    return number.real;
}
double complex_imag(complex_number_t number) {
    return number.imag;
}
complex_number_t complex_conj(complex_number_t number) {
    return (complex_number_t) {
            .real = number.real, .imag = (-1.0) * number.imag
    };
}
complex_number_t complex_exp(complex_number_t number) {
    const double x = number.real;
    const double y = number.imag;
    //Euler's formula for complex exponentials: e^z with z = x + (y*i) <=> e^x  (cos(y) + i sin(y))
    const double r = exp(x) * cos(y);
    const double i = exp(x) * sin(y);

    complex_number_t obj = {
            .real = r, .imag = i
    };
    return obj;
}


#endif //SCARLET_STUDY_CUSTOM_COMPLEX_H
