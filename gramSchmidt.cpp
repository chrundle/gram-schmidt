#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "../linearalgebra.h"


/* ----------------------- gramSchmidt ----------------------- */
/*  Given a matrix A of dimension n by n, this algorithm 
    computes a QR decomposition of A, where Q is a unitary 
    n by n matrix and R is a n by n upper triangular matrix
    and A = QR.    
    
    Input variables:
        a: pointer to array of arrays, the ith array of
            which should correspond to the ith column of the 
            matrix A. During the algorithm, the columns of Q 
            will replace the columns of A.
        r: pointer to array of arrays in which the ith 
            column of the upper triangular matrix R will be 
            stored in the ith subarray of r.
        n: number of rows and columns in A.

    Features: This implementation has time complexity O(n^3)
    and requires O(1) additional memory. 

    Remarks: Due to the nature of the problem, if A is nearly
    rank-deficient then the resulting columns of Q may not
    exhibit the orthogonality property.                        */

void gramSchmidt (double ** a, double ** r, int n) {
    int i, j;

    for(i = 0; i < n; i++) {
        r[i][i] = norm(a[i], n);                  // r_ii = ||a_i||
        scalar_div(a[i], r[i][i], n, a[i]);      // a_i = a_i/r_ii
        for(j = i+1; j < n; j++) {
            r[i][j] = dot_product(a[i], a[j], n); // r_ij = a_i*a_j
            scalar_sub(a[i], r[i][j], n, a[j]);   // a_j -= r_ij q_i
        }
    }
}


int main () {
    int i, j, n;
    double x;

    /* let user set the dimension of matrix A */
    std::cout << "Enter the dimension n (where A is a n by n matrix): ";
    std::cin >> n;

    /* allocate memory for the matrices A and R */
    double ** a = new double*[n];
    double ** r = new double*[n];
    for(i = 0; i < n; i++) {
        a[i] = new double[n];
        r[i] = new double[n];
    }

    /* initialize the values in matrix A */
    for(i = 0; i < n; i++) {
        for(j = i; j < n; j++) {
            a[i][j] = j - i + 1; // this choice of values was arbitrary
        }
    }

    /* print the matrix A before calling gramSchmidt */
    std::cout << "A = " << std::endl;
    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {

            printf("%9.7lg ", a[j][i]);
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    gramSchmidt(a, r, n); // execute gramSchmidt to determine QR factorization

    /* print the matrix Q resulting from gramSchmidt */
    std::cout << "Q = " << std::endl;
    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            if(a[j][i] >= 0) {
                std::cout << " ";
            }
            printf("%9.7lg ", a[j][i]);
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    /* print the matrix R resulting from gramSchmidt */
    std::cout << "R = " << std::endl;
    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            printf("%9.7lg ", r[i][j]);
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    /* print numerical evidence that columns of Q are orthonormal */
    printf("Numerical verification that {q_1, ..., q_%i} is an "
           "orthonormal set:\n", n);
    for(i = 0; i < n; i++) {
        for(j = i; j < n; j++) {
            x = dot_product(a[i], a[j], n);
            printf("q_%i * q_%i = %lg\n", i + 1, j + 1, x);
        }
    }
    
    /* free memory */
    for(i = 0; i < n; i++) {
        delete[] a[i];
        delete[] r[i];
    }
    delete[] a;  
    delete[] r;

    return 0;       // exit main
}
