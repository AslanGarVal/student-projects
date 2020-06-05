#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include"matplotlibcpp.h"

namespace plt = matplotlibcpp;

double mu = 1.0;//parametros del sistema

double dx(double v){ //v=dx/dt
    return v;
}

double dv(double x, double v){
    return mu*(1-x*x)*v-x;
}

int main(){
    double x0[3] = {0.0, 1.0, -2.161}; //condiciones iniciales
    double v0[3] = {2.16108, 0, 2.53}; 
    double t0 = 0;

    double dt = 0.001; //step-size

    std::vector<double> t = {t0};
    std::vector<double> x = {x0[0]};
    std::vector<double> v = {v0[0]};
    

    for (int j = 0; j <= 2; j++)
    {
        x.push_back(x0[j]); 
        v.push_back(v0[j]); 
        t.push_back(t0);

        plt::figure();

        for (int i = 0; i <= 10000; i++)
        {
            double k_1 = dx(v0[j]);
            double k_2 = dx(v0[j] + dt/2*k_1);
            double k_3 = dx(v0[j] + dt/2*k_2);
            double k_4 = dx(v0[j] + dt*k_3);

            double l_1 = dv(x0[j],v0[j]);
            double l_2 = dv(x0[j]+dt/2, v0[j]+dt/2*l_1);
            double l_3 = dv(x0[j]+dt/2, v0[j]+dt/2*l_2);
            double l_4 = dv(x0[j]+dt, v0[j]+dt*l_3);

            x0[j] += dt/6*(k_1+2*k_2+2*k_3+k_4);
            v0[j] += dt/6*(l_1+2*l_2+2*l_3+l_4);
            t0 += dt;

            t.push_back(t0);
            x.push_back(x0[j]);
            v.push_back(v0[j]);
        }
        std::string subtitle = "vanderpol", index = std::to_string(j+1), format = ".png";
        std::string title = subtitle+index+format;
        plt::named_plot("Velocidad",t,v);   
        plt::named_plot("Posicion",t,x);
        plt::xlabel("Tiempo");
        plt::title(index+" set de condiciones iniciales");
        plt::legend();
        plt::save(title);
        t0 = 0; 
        t.clear(); x.clear(); v.clear();
    }
    return 0;
}