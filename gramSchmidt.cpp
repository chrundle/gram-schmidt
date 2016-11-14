#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "linearalgebra.h"


void gramSchmidt (double ** a, double ** v, double ** q, int length) {
    int i, j, k;
    double r;

    for(i = 0; i < length; i++) {
        vec_copy(a[i], v[i], length);
    }

    for(i = 0; i < length; i++) {

        r = norm(v[i], length);
        scalar_mult(v[i], r, length, q[i]);

        for(j = i+1; j < length; j++) {
            r = dot_product(q[i], v[j], length);
            scalar_sub(q[i], r, length, v[j]);
        }
    }
}


int main () {
    int i, j, n;
    double x;

    n = 6;

    double ** a = new double*[n];
    double ** v = new double*[n];
    double ** q = new double*[n];

    for(i = 0; i < n; i++) {
        a[i] = new double[n];
        v[i] = new double[n];
        q[i] = new double[n];
    }

    for(i = 0; i < n; i++) {
        for(j = i; j < n; j++) {
            a[i][j] = i + j + 1;
        }
    }

    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            if(a[j][i] >= 0) {
                std::cout << " ";
            }
            std::cout << a[j][i] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;

    gramSchmidt(a, v, q, n);

    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            if(q[j][i] >= 0) {
                std::cout << " ";
            }
            std::cout << q[j][i] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;

    for(i = 0; i < n; i++) {
        for(j = i; j < n; j++) {
            x = dot_product(q[i], q[j], n);
            printf("q[%i] * q[%i] = %lg\n", i, j, x);
        }
    }

    delete[] a, v, q;

    return 0;
}
