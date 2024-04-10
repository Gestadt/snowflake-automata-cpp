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
    SnowFlake* CreateSnowFlake(int maxIterations, double alpha, double beta, double gamma, int n);
    bool Iterate(SnowFlake* snowFlake);
    double* GetGrid(const SnowFlake* snowFlake);
    int Test(int t);
}

#endif //SNOWFLAKE_AUTOMATA_CPP_LIBRARY_H
