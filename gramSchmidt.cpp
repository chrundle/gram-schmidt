#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "../linearalgebra.h"


/* ----------------------- gramSchmidt ----------------------- */
/*  Given a matrix A of dimension m by n, this algorithm 
    computes a QR decomposition of A, where Q is a unitary 
    m by n matrix and R is a n by n upper triangular matrix
    and A = QR.    
    
    Input variables:
        a   : pointer to array of arrays, the ith array of
                which should correspond to the ith column of the 
                matrix A. During the algorithm, the columns of Q 
                will replace the columns of A.
        r   : pointer to array of arrays in which the ith 
                column of the upper triangular matrix R will be 
                stored in the ith subarray of r.
        m   : number of columns in A.
        n   : number of rows in A.
        thin: TRUE  => thin QR factorization computed
              FALSE => full QR factorization computed

    Features: This implementation has time complexity O(m n^2)
    and requires O(1) additional memory.

    Remarks: Due to the nature of the problem, if A is nearly
    rank-deficient then the resulting columns of Q may not
    exhibit the orthogonality property.                        */

void gramSchmidt (double ** a, double ** r, int m, int n, bool full) {
    int i, j;
    double anorm, tol = 10e-7;

    for(i = 0; i < n; i++) {
        r[i][i] = norm(a[i], m);                  // r_ii = ||a_i||

        if(r[i][i] > tol) {
            scalar_div(a[i], r[i][i], m, a[i]);   // a_i = a_i/r_ii
        }
        else if(i == 0) { // set a[0] = [1 0 0 ... 0]^T
            a[i][0] = 1;
            for(j = 1; j < m; j++) {
                a[i][j] = 0;
            }
        }
        else{ // need to choose a_i orthogonal to < a_1, ... a_{i-1} >
            for(j = 0; j < m; j++) {
                a[i][j] = -a[0][i] * a[0][j];
            }
            a[i][i] += 1;

            for(j = 1; j < i; j++) {
                scalar_sub(a[j], a[j][i], m, a[i]);
            }

            anorm = norm(a[i], m);
            scalar_div(a[i], anorm, m, a[i]);
        }

        for(j = i+1; j < n; j++) {
            r[j][i] = dot_product(a[i], a[j], m); // r_ij = a_i*a_j
            scalar_sub(a[i], r[j][i], m, a[j]);   // a_j -= r_ij a_i
        }
    }

    /* if full QR factorization requested, we choose remaining 
       columns of Q so that the m columns of Q form an 
       orthonormal set                                          */
    if(full) {
        for(; i < m; i++) {
            for(j = 0; j < m; j++) {
                    a[i][j] = -a[0][i] * a[0][j];
                }
                a[i][i] += 1;
    
                for(j = 1; j < i; j++) {
                    scalar_sub(a[j], a[j][i], m, a[i]);
                }
    
                anorm = norm(a[i], m);
                scalar_div(a[i], anorm, m, a[i]);
        }
    }
}


int main () {
    int i, j, n, m, q_n, r_m;
    bool full;
    double x;

    /* let user set the dimension of matrix A */
    std::cout << "Enter the dimension m (where A is a m by n matrix): ";
    std::cin >> m;
    std::cout << "Enter the dimension n (where A is a m by n matrix): ";
    std::cin >> n;

    if(m != n) {
        /* check if m < n */
        if(m < n) {
            printf("For a successful factorization, this implementation "
                   "requires n <= m.\nTerminating program.\n");
            return 0;
        }
        /* let user choose either full or thin QR factorization */
        std::cout << "Enter either 0 to compute a thin QR factorization"
                  << std::endl;
        std::cout << "          or 1 to compute a full QR factorization: ";
        std::cin >> full;
    }
    else { // else m == n so full and thin QR factorization are identical */
        full = 1;
    }

    /* set dimensions of matrices Q and R based on full or thin QR */
    if(full) { // Q is m by m and R is m by n
        q_n = m;
        r_m = m;
    }
    else { // Q is m by n and R is n by n
        q_n = n;
        r_m = n;
    }

    /* allocate memory for the matrices A and R */
    double ** a = new double*[q_n];
    double ** r = new double*[n];
    for(i = 0; i < n; i++) {
        a[i] = new double[m];
        r[i] = new double[r_m];
    }
    for(; i < q_n; i++) {
        a[i] = new double[m];
    }

    /* initialize the values in matrix A (only n columns regardless of
       thin QR or full QR) */
    for(i = 0; i < n; i++) {
        for(j = i; j < m; j++) {
            a[i][j] = j - i + 1; // this choice of values was arbitrary
        }
    }

    /* print the matrix A before calling gramSchmidt */
    std::cout << "A = " << std::endl;
    for(i = 0; i < m; i++) {
        for(j = 0; j < n; j++) {
            printf("%9.6lg ", a[j][i]);
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    /* execute gramSchmidt to compute QR factorization */
    gramSchmidt(a, r, m, n, full);

    /* print the matrix Q resulting from gramSchmidt */
    std::cout << "Q = " << std::endl;
    for(i = 0; i < m; i++) {
        for(j = 0; j < q_n; j++) {
            if(a[j][i] >= 0) {
                std::cout << " ";
            }
            printf("%9.6lg ", a[j][i]);
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    /* print the matrix R resulting from gramSchmidt */
    std::cout << "R = " << std::endl;
    for(i = 0; i < r_m; i++) {
        for(j = 0; j < n; j++) {
            printf("%9.6lg ", r[j][i]);
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    /* print numerical evidence that columns of Q are orthonormal */
    printf("Numerical verification that {q_1, ..., q_%i} is an "
           "orthonormal set:\n", q_n);
    for(i = 0; i < q_n; i++) {
        for(j = i; j < q_n; j++) {
            x = dot_product(a[i], a[j], m);
            printf("q_%i * q_%i = %lg\n", i + 1, j + 1, x);
        }
    }
    
    /* free memory */
    for(i = 0; i < n; i++) {
        delete[] a[i];
        delete[] r[i];
    }
    for(; i < q_n; i++) {
        delete[] a[i];
    }
    delete[] a;  
    delete[] r;

    return 0;       // exit main
}
