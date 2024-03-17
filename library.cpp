#include "library.h"
#include <iostream>
#include <array>

SnowFlake::SnowFlake(int maxIterations, int alpha, int beta, int gamma, int n)
{
    _max_iterations = maxIterations;
    _alpha = alpha;
    _beta = beta;
    _gamma = gamma;
    _n = n;
    _s = new double[n*n];
    _c = new double[n*n];
    _r = new double[n*n];
    _nr = new double[n*n];
    _av = new double[n*n];

}

SnowFlake* CreateSnowFlake(int maxIterations, int alpha, int beta, int gamma, int n)
{
    const auto snowFlake =  new SnowFlake(maxIterations, alpha, beta, gamma, n);
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
            snowFlake->_s[i*n+j] = beta;
    const int mid  = (n-1)/2;
    snowFlake->_s[mid*n+mid] = 1;
    return snowFlake;
}

bool Iterate(const SnowFlake* snowFlake) {
    int n = snowFlake->_n;
    double gamma = snowFlake->_gamma;
    double alpha = snowFlake->_alpha;
    double* r = snowFlake->_r;
    double* nr = snowFlake->_nr;
    double* av = snowFlake->_av;
    double* s = snowFlake->_s;
    double* c = snowFlake->_c;

    // Initializing r, nr, and av to zero
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            r[i*n+j] = 0;
            nr[i*n+j] = 0;
            av[i*n+j] = 0;
            c[i*n+j] = 0;
        }
    }

    // Modifying c based on conditions
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (s[i*n+j] >= 1) {
                c[i*n+j] = 1;
                c[(i + 1) % n * n+j] = 1;
                c[(i - 1 + n) % n * n +j] = 1;
                c[i*n+(j - 1 + n) % n] = 1;
                c[i*n+(j + 1) % n] = 1;
                if (i % 2 == 0) {
                    c[(i + 1) % n * n + (j - 1 + n) % n] = 1;
                    c[(i - 1 + n) % n * n + (j - 1 + n) % n] = 1;
                } else {
                    c[(i + 1) % n * n + (j + 1) % n] = 1;
                    c[(i - 1 + n) % n * n + (j + 1) % n] = 1;
                }
            }
        }
    }

    // Calculating r and nr based on conditions
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (c[i*n + j]>0) {
                r[i*n+j] = s[i*n+j] + gamma;
            } else {
                nr[i*n+j] = s[i*n+j];
            }
        }
    }

    // Calculating av based on conditions
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i % 2 == 0) {
                av[i*n+j] = nr[i*n+j] * (1 - alpha * 0.5) +
                            nr[i*n+(j - 1 + n) % n] * (alpha / 12) +
                            nr[i*n+(j + 1) % n] * (alpha / 12) +
                            nr[(i - 1 + n) % n *n + j] * (alpha / 12) +
                            nr[(i + 1) % n * n + j] * (alpha / 12) +
                            nr[(i + 1) % n *n + (j - 1 + n) % n] * (alpha / 12) +
                            nr[(i - 1 + n) % n *n + (j - 1 + n) % n] * (alpha / 12);
            } else {
                av[i*n+j] = nr[i*n+j] * (1 - alpha * 0.5) +
                            nr[i*n+(j - 1 + n) % n] * (alpha / 12) +
                            nr[i*n+(j + 1) % n] * (alpha / 12) +
                            nr[(i - 1 + n) % n * n + j] * (alpha / 12) +
                            nr[(i + 1) % n * n + j] * (alpha / 12) +
                            nr[(i + 1) % n * n + (j + 1) % n] * (alpha / 12) +
                            nr[(i - 1 + n) % n * n + (j + 1) % n] * (alpha / 12);
            }
        }
    }

    // Assigning s based on r and av
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            s[i*n+j] = r[i*n+j] + av[i*n+j];
        }
    }


    return true;
}


double* GetGrid(const SnowFlake* snowFlake) {
    return snowFlake->_s;
}

