#include <iostream>
#include <cmath>
#include <gauss_legendre.h> // compilar usando "g++ 'nombre.cpp' gauss_legendre.cpp -I '/direccion/del/.h'"

double g(double x){
    return 1/sqrt(1-x*x);
}

double f(double x){
    return 1/x;
}

int main(){
    for (int i = 1; i <= 10; i++)
    {
        GaussLegendreIntegrator I = GaussLegendreIntegrator(i);
        double integral_g = I.integrate(-1.0, 1.0, g), integral_f = I.integrate(1.0, 2.0, f);
        std::cout << "integral de 1/x = " << integral_f << ", integral de 1/sqrt(1-x^2) = " << integral_g << " en n = " << i << std::endl;
    }
    
    return 0;
}