#include <iostream>
#include <vector>
#include <cmath>
#include <gauss_legendre.h> // compilar usando "g++ 'nombre.cpp' gauss_legendre.cpp -I '/direccion/del/.h'"

double dy(double v){ //v=dy/dt
    return v;
}

double dv(double x, double v){
    return -2.0*x*v;
}

double f(double x){
    return 1/sqrt(M_PI)*exp(-x*x);
}

GaussLegendreIntegrator I(10);

int main(){

    /************************************* Solucion con IVP *****************************************************/
    std::cout << " " << std::endl;
    std::cout << "/******** Solución con RK4 *******/" << std::endl;
    std::cout << " " << std::endl;

    double h = 0.01; //step-size

    double y0 = 0.0, v0 = 2/sqrt(M_PI), x0 = 0;  //condiciones iniciales

    //generar donde guardaremos el vector solucion
    std::vector<double> x, v, y;
    x.push_back(x0);
    v.push_back(v0);
    y.push_back(y0);

    //ahora a generar las y, en realidad es sist de 2 ecs.
    for(int i = 0; i <= 600; i++){
        double k_1 = dy(v0);
        double k_2 = dy(v0 + h/2*k_1);
        double k_3 = dy(v0 + h/2*k_2);
        double k_4 = dy(v0 + h*k_3);

        double l_1 = dv(x0,v0);
        double l_2 = dv(x0+h/2, v0+h/2*l_1);
        double l_3 = dv(x0+h/2, v0+h/2*l_2);
        double l_4 = dv(x0+h, v0+h*l_3);

        y0 += h/6*(k_1+2*k_2+2*k_3+k_4);
        v0 += h/6*(l_1+2*l_2+2*l_3+l_4);
        x0 += h;

        x.push_back(x0);
        y.push_back(y0);
        v.push_back(v0);
    }

    //imprime los valores de erf en los vals. pedidos (usamos la antiparidad de la función para ahorrarnos unos calculos)
    std::cout << "erf(3) = " << y[300] << " = -erf(-3)" << std::endl; 
    std::cout << "erf(2) = " << y[200] << " = -erf(-2)" << std::endl; 
    std::cout << "erf(1) = " << y[100] << " = -erf(-1)" << std::endl; 
    std::cout << "erf(0) = " << y[0] << std::endl; 

    /************************************** Solución con cuadratura Gauss-Legendre de orden 10 ***********************************/
    std::cout << " " << std::endl;
    std::cout << "/******** Solución con cuadratura Gauss-Legendre de orden 10 *******/" << std::endl;
    std::cout << " " << std::endl;

    for (int i = 0; i < 4; i++)
    {
        std::cout << "erf(" << i << ") = " << I.integrate(-(double) i, (double) i, f) << " = -erf(-" << i << ")" << std::endl;
    }
    
    return 0;
}