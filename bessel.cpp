#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <gauss_legendre.h> // compilar usando "g++ 'nombre.cpp' gauss_legendre.cpp -I '/direccion/del/.h'"

GaussLegendreIntegrator I(20); //Si usamos n=10, no es muy preciso :/, usar otro metodo es todavía más tardado

class Bessel{//Aquí generamos la función de Bessel de orden n
    private:

    int order;

    public:

    Bessel(int j): order(j){} //Constructor

    static double aux_func(int n, double t, double x){//funcion para la representacion integral
        Bessel j = Bessel(n);
        double p = cos(j.order*t-x*sin(t));
        return 1/M_PI*p;
    }
    
    double eval(double x){ //evalúa usando la representación integral
        return I.special_integrate(0, M_PI, aux_func, x, this->order);
    }
};

class Intervalo{
    private:

    float a,b;
    public:

    Intervalo(float c, float d): a(c), b(d){}
    void resize_left(float x){
        this->a = x;
    }
    void resize_right(float x){
        this->b = x;
    }
    float mid_point(){
        return (a+b)/2;
    }
    float get_start(){
        return a;
    }
    float get_end(){
        return b;
    }
};

double bisection(Bessel P, Intervalo I, double tol){//Implementación del método de bisección
    double x0 = I.mid_point(), a = I.get_start(), b = I.get_end();
    int iterations = 0;
    while (std::abs(I.get_end()-I.get_start()) > tol && iterations <= 50)
    {
        if (P.eval(x0) == 0.0)
        {
            break;
        } else
        {
            if (P.eval(a)*P.eval(x0) < 0){
                I.resize_right(x0);
                x0 = I.mid_point();
            }
            else if (P.eval(x0)*P.eval(b) < 0){
                I.resize_left(x0);
                x0 = I.mid_point();
            }
        }
        iterations ++;
    }
    return x0;
}



int main(){
    std::vector<Bessel> J;
    std::vector<double> zeros;
    for (int i = 0; i < 6; i++)
    {
        J.push_back(Bessel(i));
    }
    double error = 0.000001;
    std::cout.precision(7);

    std::ofstream output;
    output.open("BesselRoots.txt");
    
    /* Comencemos con J(0), por inspección gráfica, sus ceros están contenidos en subintervalos de largo 1, empezando en 2. 
    Calcularemos los primeros 9, pues usaremos la misma propiedad de interlacing para calcular el resto. 
    Empieza con más ceros de los necesarios para generar las siguientes particiones (cada iteración obtiene un cero menos) */
    std::vector<Intervalo> Partition;
    for (int i = 1; i < 13; i++)
    {
        Partition.push_back(Intervalo(2.5*((float)i-1),2.5*((float)i)));
    }
    output << "n = 0:" << std::endl; 
    for (Intervalo &I: Partition)
    {
        if (std::abs(J[0].eval(bisection(J[0], I, error))) < error)
        {
            zeros.push_back(bisection(J[0], I, error));
        }   
    }
    for (int i = 0; i < zeros.size()-1; i++)
    {
        output << zeros[i] << ",";
    } 
    output << zeros[zeros.size()-1] << std::endl;
    
    for (int i = 1; i < 5; i++) 
    {
        Partition.clear();
        for (int k = 1; k < zeros.size(); k++)
        {
            Partition.push_back(Intervalo(zeros[k-1], zeros[k]));
        }
        zeros.clear();
        for (Intervalo &I: Partition)
        {
            if (std::abs(J[i].eval(bisection(J[i], I, error))) < error)
            {
                zeros.push_back(bisection(J[i], I, error));
            }   
        }
        output << "n = " << i << ":" << std::endl; 
        for (int i = 0; i < zeros.size()-1; i++)
        {
            output << zeros[i] << ",";
        }    
        output << zeros[zeros.size()-1] << std::endl; 
    }
    return 0;
}