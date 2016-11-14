double norm (double * x, int length) {
    int i;
    double a, sum;

    for(i = 0; i < length; i++) {
        a = x[i];
        sum += a * a;
    }

    return sqrt(sum);
}

void vec_copy (double * x, double * y, int length) {
    int i;

    for(i = 0; i < length; i++) {
        y[i] = x[i];
    }
}

void scalar_mult (double * x, double r, int length, double * y) {
    int i;

    for(i = 0; i < length; i++) {
        y[i] = x[i]/r;
    }
}

void scalar_sub (double * x, double r, int length, double * y) {
    int i;

    for(i = 0; i < length; i++) {
        y[i] -= r * x[i];
    }
}

double dot_product (double * x, double * y, int length) {
    int i;
    double a, sum;

    for(i = 0; i < length; i++) {
        sum += x[i] * y[i];
    }

    return sum;
}
