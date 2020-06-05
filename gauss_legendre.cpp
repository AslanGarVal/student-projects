#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string> 
#include <utility>
#include <gauss_legendre.h>

std::pair<std::vector<double>, std::vector<double>> set_nodes(int m){ //obtiene las raíces y los pesos de los archivos "roots" y "weights"
    std::ifstream input_n, input_w;
    input_n.open("roots.txt");
    input_w.open("weights.txt");
    std::vector<double> nodes, weights;
    std::pair <std::vector<double>, std::vector<double>> parcito;

    int line_counter = 1, line_counter_n = 1;
    std::string line_w, line_n;

    if (input_w.is_open()){
        while (std::getline(input_w, line_w))
        {
            if (line_counter == 2*(m+1))
            {
                while (line_w.find(",") != std::string::npos)
                {
                    std::string weight = line_w.substr(0,line_w.find(","));
                    weights.push_back(std::stof(weight));
                    line_w.erase(0, weight.length()+1);
                }
                weights.push_back(std::stof(line_w));
            }   
            line_counter ++;
        }
    }

    if (input_n.is_open()){
        while (std::getline(input_n, line_n))
        {
            if (line_counter_n == 2*(m+1))
            {
                while (line_n.find(",") != std::string::npos)
                {
                    std::string node = line_n.substr(0,line_n.find(","));
                    nodes.push_back(std::stof(node));
                    line_n.erase(0, node.length()+1);
                }
                nodes.push_back(std::stof(line_n));
            }   
            line_counter_n ++;
        }
    }
    parcito.first = nodes;
    parcito.second = weights;
    return parcito;
}

GaussLegendreIntegrator::GaussLegendreIntegrator(int j): order(j){//Constructor del integrador
        for(double &n: set_nodes(j).first){
            this->nodes.push_back(n);
        }
        for(double &w: set_nodes(j).second){
            this->weights.push_back(w);
        }
}

double GaussLegendreIntegrator::integrate(double a, double b, double function(double)){//La parte que integra
        double integral = 0;
        for (int i = 0; i < this-> order; i++)
        {
            double s = this->weights[i]*function((b-a)/2*this->nodes[i]+(a+b)/2);
            integral += s;
        }

        return (b-a)/2*integral; 
    }

double GaussLegendreIntegrator::special_integrate(double a, double b, double function(int, double, double), double x, int n){//Arreglo feo para integrar Bessel con G-L
        double integral = 0;
        for (int i = 0; i < this-> order; i++)
        {
            double s = this->weights[i]*function(n, (b-a)/2*this->nodes[i]+(a+b)/2, x);
            integral += s;
        }

        return (b-a)/2*integral; 
    }

