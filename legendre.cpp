#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

std::vector<double> set_coeffs(int m){ //funcion para obtener los coeficientes, se invoca en el constructor
    if (m == 0) 
    {               
        return {1};
    }
    if (m == 1)
    {
        return {0, 1};
    }

    std::vector<double> coeffs(m+1);

    std::vector<double> v = set_coeffs(m-1);
    std::vector<double> u = set_coeffs(m-2);

    double a = (2.0*m-1.0)/m;
    double b = (m-1.0)/m;

    int first = 1;
    
    if ( m % 2 == 0 )
    {
        coeffs[0] = -b * u[0];
        first = 2;
    }
    for (int i = first; i < m-1; i += 2)
    {
        coeffs[i] = (a * v[i-1] - b * u[i]);
    }
    coeffs[m] = a * v[m - 1];

    return coeffs;
}


class LegendrePolynomial{
    private:

    int n;
    std::vector <double> coeffs;

    public:

    LegendrePolynomial(int j): n(j){
        for(double c: set_coeffs(j)){
            this-> coeffs.push_back(c);
        }
    }

    void get_coeffs(){ // Para imprimir los coeficientes, invoquese este metodo
        for(auto &c : coeffs){
            std::cout << c << ", " ;
        }
        std::cout << " " << std::endl;
    }

    double eval(double x){ //evalua en x
        double p = 0;
        for (int i = 0; i < this-> n+1; i++)
        {
            double c = this-> coeffs[i];
            p += c*pow(x,i);
        }
        
        return p;
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

double derivative(int n, double x){ //Calcula la derivada del P(n)
    double f = n/(x*x-1);
    LegendrePolynomial P(n), Q(n-1);
    return f*(x*P.eval(x)-Q.eval(x));
}

double get_weight(int n, double x){//Obtiene los pesos para la cuadratura Gauss-Legendre
    double p = derivative(n, x);
    return 2/((1-x*x)*p*p);
}

double bisection(LegendrePolynomial P, Intervalo I, double tol){//Implementación del método de bisección
    double x0 = I.mid_point(), a = I.get_start(), b = I.get_end();

    while (std::abs(I.get_end()-I.get_start()) > tol)
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
    }
    return x0;
}

double fixed_point(double x){
    LegendrePolynomial P(3);
    return 2.0/3.0*(P.eval(x) + 3.0/2.0*x);
}

int main(){
    int n; //Esta parte da los coeficientes del polinomio "n"
    std::cin >> n;
    LegendrePolynomial P(n);
    P.get_coeffs(); 
 
    /****************************************** CEROS **********************************************************************/

    std::ofstream output, output2;
    output.open("roots.txt");
    output2.open("weights.txt");

    std::vector<LegendrePolynomial> Pols;
    for (int i = 0; i <= 20; i++)
    {
        LegendrePolynomial p = LegendrePolynomial(i);
        Pols.push_back(p);
    }
    
    output << "n = 0:\nno hay raíces" << std::endl; //caso n=0
    output2 << "n = 0:\nno hay nodos" << std::endl;
    output << "n = 1:\n0.0" << std::endl; //caso n=1
    output2 << "n = 1:\n2.0" << std::endl;

    double precision = 0.000001;
    std::cout.precision(7);

    // Aplica Método de Newton para P2

    double x0 = 0.5;
    while (std::abs(Pols[2].eval(x0)) > precision)
    {
        double q = Pols[2].eval(x0)/derivative(2,x0);
        x0 -= q;
    }

    output << "n = 2:\n" << "-" << x0 << "," << x0 << std::endl;
    output2 << "n = 2:\n" << get_weight(2,-x0) << "," << get_weight(2,x0) << std::endl;
    
    //Aplicaremos punto fijo para P3

    double root4 = 0.6, prev = root4;
    while (1)
    {
        prev = fixed_point(root4);
        if (std::abs(root4-prev) < precision)
        {
            output << "n=3:\n";
            break;
        }
        
        root4 = prev;
    }

    Intervalo I(0.1,0.9);
    output << "-" << bisection(Pols[3], I, precision) << "," << prev << "," << bisection(Pols[3], I, precision) << std::endl;
    output2 << "n = 3:" << std::endl;
    output2 << get_weight(3,-bisection(Pols[3], I, precision)) << "," << get_weight(3, 0) << "," << get_weight(3,bisection(Pols[3], I, precision)) << std::endl;
    std::vector<double> zeros = {-bisection(Pols[3], I, precision),0, bisection(Pols[3],I,precision)};
    std::vector<double> weights = {get_weight(3,zeros[0]), get_weight(3,zeros[1]), get_weight(3,zeros[2])};

    /*
      Para el resto usaremos la propiedad de interlacing: los ceros de P(n) generan una particion del intervalo (-1,1) donde hay
      un único cero de P(n+1) en cada subintervalo de la particion. Aplicamos bisección en cada subintervalo.
    */
    for (int n = 4; n <= 20; n++)
    { 
        std::vector<Intervalo> Partition;

        Intervalo p_i = Intervalo(-1.0, zeros[0]);
        Partition.push_back(p_i);
        for (int i = 0; i < n-2; i++)
        {
            Intervalo p = Intervalo(zeros[i], zeros[i+1]);
            Partition.push_back(p);
        }
        Intervalo p_f = Intervalo(zeros[n-2], 1.0);
        Partition.push_back(p_f);

        zeros.clear();
        weights.clear();

        for(Intervalo &inter: Partition){
            zeros.push_back(bisection(Pols[n], inter, precision));
            weights.push_back(get_weight(n,bisection(Pols[n], inter, precision)));
        }
       
        output << "n = " << n << ":" << std::endl; 
        output2 << "n = " << n << ":" << std::endl; 
        for (int i = 0; i < n-1; i++)
        {
            output << zeros[i] << ",";
            output2 << weights[i] << ",";
        }
        output << zeros[n-1] << std::endl;
        output2 << weights[n-1] << std::endl;
    }

    output.close();
    output2.close();

    return 0; //La parte de la cuadratura y las integrales viene aparte en gauss_legendre.cpp, gauss_legendre.h, integrales.cpp
}