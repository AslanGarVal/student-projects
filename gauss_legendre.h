#ifndef _gauss_legendre_h
#define _gauss_legendre_h

#include <vector>

class GaussLegendreIntegrator{
    private:

    int order;
    std::vector<double> nodes;
    std::vector<double> weights;

    public:

    GaussLegendreIntegrator(int j);
    double integrate(double a, double b, double function(double));
    double special_integrate(double a, double b, double function(int, double, double),double x, int n);
};

#endif