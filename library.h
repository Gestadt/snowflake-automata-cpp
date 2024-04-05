#ifndef SNOWFLAKE_AUTOMATA_CPP_LIBRARY_H
#define SNOWFLAKE_AUTOMATA_CPP_LIBRARY_H

class SnowFlake {
private:

public:
    SnowFlake(int maxIterations, double alpha, double beta, double gamma, int n);
    int _max_iterations;
    double _alpha{};
    double _beta{};
    double _gamma{};
    int _n{};
    double* _s; // linear storage of nxn grid
    double* _c; // temp
    double* _r; // receptive
    double* _nr; // not receptive
    double* _av; // average
    int _iteration;
};

extern "C"
{
    __attribute__((visibility("default"))) SnowFlake* CreateSnowFlake(int maxIterations, double alpha, double beta, double gamma, int n);
    __attribute__((visibility("default"))) bool Iterate(SnowFlake* snowFlake);
    __attribute__((visibility("default"))) double* GetGrid(const SnowFlake* snowFlake);
}

#endif //SNOWFLAKE_AUTOMATA_CPP_LIBRARY_H
