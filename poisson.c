int64_t random_poisson(bitgen_t *bitgen_state, double lam) {
    if (lam >= 10.0) {
        // For large lambda, use Hörmann's Transformed Rejection Method
        return random_poisson_ptrs(bitgen_state, lam);
    } else if (lam == 0.0) {
        return 0;
    } else {
        // For small lambda, use the sequential multiplication method
        return random_poisson_mult(bitgen_state, lam);
    }
}


int64_t random_poisson_mult(bitgen_t *bitgen_state, double lam) {
    int64_t X;
    double prod, U, enlam;

    enlam = exp(-lam);
    X = 0;
    prod = 1.0;

    while (1) {
        U = next_double(bitgen_state); // Generates a uniform random float in [0, 1)
        prod *= U;
        if (prod <= enlam) {
            return X;
        }
        X++;
    }
}

int64_t random_poisson_ptrs(bitgen_t *bitgen_state, double lam) {
    double slam, loglam, a, b, invalpha, vr, us;
    int64_t k;

    slam = sqrt(lam);
    loglam = log(lam);
    b = 0.931 + 2.53 * slam;
    a = -0.059 + 0.02483 * b;
    invalpha = 1.1239 + 1.1328 / (b - 3.4);
    vr = 0.9277 - 0.62241 / slam;

    while (1) {
        double U = next_double(bitgen_state) - 0.5;
        double V = next_double(bitgen_state);
        double u = 0.5 - fabs(U);

        if ((U < 0.0) && (u < 0.0)) continue;

        us = 0.5 - U * U / u;
        k = (int64_t)floor((2.0 * a / u + b) * U + lam + 0.43);

        // Immediate acceptance via the "squeeze" boundary conditions
        if ((us >= 0.0) && (V <= vr)) {
            return k;
        }
        if ((k < 0) || ((us < 0.0) && (V >= us))) {
            continue;
        }

        // Mathematical acceptance step (safeguard log evaluation)
        if (log(V * invalpha / (a / (u * u) + b)) <=
            k * loglam - lam - log_factorial(k)) {
            return k;
        }
    }
}
