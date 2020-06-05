#ifndef INCLUDE_MATRIX
#define INCLUDE_MATRIX
#include<iostream>
#include<vector>
#include<string>
#include<cstring>
#include<utility>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<iterator>

#include<random> //To get the mersenne twister engine
#include<chrono> //To get time

////////////// Generador aleatorio compartido //////////////////
class RandGen{
public:
    static RandGen& get(){
        if(_instance == nullptr){
            _instance = new RandGen();
        }
        return *_instance;
    }
    std::mt19937 gen;
private:
    static RandGen* _instance;
    RandGen(){
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        gen = std::mt19937(seed);
    };
    ~RandGen(){};
};

RandGen* RandGen::_instance = nullptr;
/////////////////////////////////////////////////////////////////

class Matrix{
public:
    // Primero ponemos los constructores
    // Recibiendo m y n por aparte
    Matrix(unsigned short m, unsigned short n);
    // El trivial, construye la matriz de 1x1 con 0 como valor
    Matrix(): Matrix(1,1){};
    // Recibiendo m y n en un std::pair, reutilizamos el constructor de arriba
    Matrix(std::pair<unsigned short, unsigned short> shape): Matrix(shape.first, shape.second){};
    // Constructor de lista, el más poderoso
    Matrix(std::initializer_list< std::initializer_list<double> > vals);

    // Métodos numéricos que sólo aplican a matrices cuadradas.
    // Método que intenta resolver la ecuación Ax=b, siendo A la instancia que llama
    Matrix simple_solve(Matrix &b);
    // Del par que debe regresar, la entrada "first" debe ser la matriz G, y "second" A^{(n-1)}
    std::pair<Matrix, Matrix> GEM();
    // Del par que debe regresar, la entrada "first" debe ser la matriz L, y "second" U (a.k.a. A^{(n-1)})
    std::pair<Matrix, Matrix> LU();

    // Método para obtener una copia de la matriz pero traspuesta
    Matrix T();

    // Método para obtener una referencia a una entrada de la matriz
    double & get(int i, int j);

    // Convierte la matriz en una cadena de caracteres
    std::string print();

    // Guarda la matriz en un archivo
    void save(std::string file_name);

    // Genera una matriz identidad
    static Matrix identity(int n);

    // Genera una matriz diagonal a partir de una matriz columna
    static Matrix diag(Matrix &d);

    // Genera un par de matrices A y b que satisfacen la ecuación Ax = b
    static std::pair<Matrix,Matrix> from_solution(Matrix &x);

    // Carga una matriz guardada como un CSV.
    static Matrix load_from_file(std::string file_path);

    // Entrega las dimensiones de la matriz en un par.
    std::pair<unsigned short, unsigned short> shape() const { return std::pair<unsigned short, unsigned short>(this->m, this->n);}
protected:
    // Vector de vectores (renglones) de la matriz
    std::vector<std::vector<double>> values;
    // Dimensiones de la matriz
    unsigned short m,n;
};

// Estos operadores son los que afectan a la matriz por la izquierda //
std::ostream& operator<<(std::ostream& os, Matrix obj){
    os << obj.print();
    return os;
}
///////////////////////////////////////////////////////////////////////

// Cambios que se hicieron
// -> Correción del print, usaba dos veces el mismo índice
// -> simple_solve usaba b en lugar de solution

Matrix::Matrix(unsigned short m, unsigned short n): m(m), n(n){
    for(int i = 0; i < m; i++){
        this->values.push_back(std::vector<double>(n));
    }
}

// Constructor de lista, el más poderoso
Matrix::Matrix(std::initializer_list< std::initializer_list<double> > vals){
    this->m = vals.size();
    this->n = vals.begin()->size();
    for(std::vector<double> row: vals){
        if (row.size() != this->n){
            throw std::string("Uno de los renglones tiene un largo distinto");
        }
        this->values.push_back(row);
    }
}

/////////// Ahora vamos con los métodos de la matriz cuadrada ///////////
Matrix Matrix::simple_solve(Matrix &b){
    if(this->m != this->n) throw std::string("La matriz debe ser cuadrada!");
    bool is_upper_triangular = true;
    for(int i = 1; i < this->m; i++){
        for(int j = 0; j < i-1; j++){
            if(this->get(i,j) != 0){
                is_upper_triangular = false;
                break;
            }
        }
        if(!is_upper_triangular) {
            break;
        }
    }

    bool is_lower_triangular = true;
    for(int i = 0; i < this->m-1; i++){
        for(int j = i+1; j < this->n; j++){
            if(this->get(i,j) != 0){
                is_lower_triangular = false;
                break;
            }
        }
        if(!is_lower_triangular) {
            break;
        }
    }
    Matrix solution(this->m, 1);
    if(is_upper_triangular && is_lower_triangular){
        for(int i = 0; i < this->m; i++){
            solution.get(i,0) = b.get(i,0)/this->get(i,i);
        }
    } else if(is_upper_triangular) {
        for(int i = this->m-1; i>=0; i--){
            double sumita = 0;
            for(int j = i+1; j < this->m; j++){
                sumita+=this->get(i,j)*solution.get(j,0);
            }
            solution.get(i,0) = (1.0/this->get(i,i))*(b.get(i,0) - sumita);
        }
    } else if(is_lower_triangular){
        for(int i = 0; i<this->m; i++){
            double sumita = 0;
            for(int j = 0; j < i; j++){
                sumita+=this->get(i,j)*solution.get(j,0);
            }
            solution.get(i,0) = (1.0/this->get(i,i))*(b.get(i,0) - sumita);
        }
    } else {
        throw std::string("Esa matriz es demasiado complicada");
    }
    return solution;
}

std::pair<Matrix, Matrix> Matrix::GEM(void){
    if(this->m != this->n) throw std::string("La matriz debe ser cuadrada!");
    Matrix G = Matrix::identity(this->m), A = *this;
    for (int k = 1; k < A.n; k++)
        {
            //implementación del pivoteo
            float max = std::abs(A.get(k-1,k-1)); 
            int max_row = k-1;
            for (int m = k-1; m < A.n; m++)
            {
                if (std::abs(A.get(m,k-1)) > max)
                {
                    max = std::abs(A.get(m,k-1));
                    max_row = m;
   
                    for (int n = 0; n < A.n; n++)
                    {
                        double tmp = A.get(max_row,n);
                        double tmpG = G.get(max_row,n);
                        A.get(max_row,n) = A.get(k-1,n);                        
                        G.get(max_row,n) = G.get(k-1,n);
                        A.get(k-1,n) = tmp;
                        G.get(k-1,n) = tmpG;    
                    } 
                }
            } 

            //implementación de la eliminación gaussiana
            for (int i = k; i < G.n; i++) 
            {
                double a = A.get(i,k-1)/A.get(k-1,k-1);
                for (int j = 0; j < G.n; j++)
                {
                    G.get(i,j) -= a*G.get(k-1,j);
                    A.get(i,j) -= a*A.get(k-1,j); 
                }     
            }
        } 
    // Hagan magia aquí
    std::pair<Matrix, Matrix> result(G,A);
    return result;
}

std::pair<Matrix, Matrix> Matrix::LU(void){
    if(this->m != this->n) throw std::string("La matriz debe ser cuadrada!");
    Matrix L = Matrix::identity(this-> m), A = *this;

        for (int k = 1; k < A.n; k++)
        {
            //implementación de la factorización LU (sin pivoteo!!)
            for (int i = k; i < L.n; i++) 
            {
                double a = A.get(i,k-1)/A.get(k-1,k-1);
                L.get(i,k-1) = a;
                for (int j = 0; j < L.n; j++)
                {
                    A.get(i,j) -= a*A.get(k-1,j); 
                }      
            }
        } 
    // Hagan magia aquí
    std::pair<Matrix, Matrix> result(L,A);
    return result;
}

Matrix Matrix::T(){
    Matrix result(this->n, this->m);
    for(int i = 0; i < this->m; i++){
        for(int j = 0; j < this->n ; j++){
            result.get(j,i) = this->get(i,j);
        }
    }
    return result;
}

// Método para obtener una referencia a una entrada de la matriz
double & Matrix::get(int i, int j){
    if (i >= this->m || j >= this->n) throw std::string("Los índices exceden las dimensiones de la matriz.");
    return this->values[i][j];
}

// Convierte la matriz en una cadena de caracteres
std::string Matrix::print() {
    std::string result;
    for(int i = 0; i< this->m; i++){
        for(int j = 0; j < this->n; j++){
            result+= std::to_string(this->get(i,j)) + ",\t";
        }
        if(i<this->m-1) result += '\n';
    }
    return result;
}

void Matrix::save(std::string file_name){
    std::ofstream f(file_name);
    if(!f.good()) throw std::string("No se pudo guardar la matriz");
    std::stringstream ss;
    for(int i = 0; i < this->m; i++){
        for(int j = 0; j < this->n; j++){
            ss << this->get(i,j) << (j < this->n-1 ? "," : (i < this->m-1 ? "\n": ""));
        }
    }
    std::copy(
        std::istreambuf_iterator<char>(ss),
        std::istreambuf_iterator<char>(),
        std::ostreambuf_iterator<char>(f)
    );
    f.close();
}

Matrix Matrix::identity(int n){
    Matrix result(n,n);
    for(int i = 0; i < n; i++) result.get(i,i) = 1;
    return result;
};

Matrix Matrix::diag(Matrix &d){
    if(d.n != 1) throw std::string("La diagonal debe ser una matriz columna");
    Matrix result(d.m, d.m);
    for(int i = 0; i < d.m; i++){
        result.get(i,i) = d.get(i,0);
    }
    return result;
}

std::pair<Matrix,Matrix> Matrix::from_solution(Matrix &x){
    if (x.n != 1) throw std::string("La solución debe ser una matriz columna!");
    std::pair<unsigned short, unsigned short> shape = x.shape();
    std::pair<Matrix, Matrix> result = std::pair<Matrix,Matrix>(
        Matrix(shape.first, shape.first),
        Matrix(shape.first, 1)
    );

    std::uniform_int_distribution<> distr(-17, 17);
    // El generador de aleatorios
    RandGen &rg = RandGen::get();

    for(int i = 0; i < shape.first; i++){
        double b = 0;
        for(int j = 0; j < shape.first; j++){
            result.first.get(i,j) = distr(rg.gen);
            b+= result.first.get(i,j) * x.get(j,0);
        }
        result.second.get(i,0) = b;
    }
    return result;
}

std::vector<std::string> inline split(const std::string &source, const char *delimiter = ",", bool keepEmpty = false){
    std::vector<std::string> results;
    size_t prev = 0;
    size_t next = 0;
    while ((next = source.find_first_of(delimiter, prev)) != std::string::npos){
        if (keepEmpty || (next - prev != 0)) results.push_back(source.substr(prev, next - prev));
        prev = next + 1;
    }
    if (prev < source.size()) results.push_back(source.substr(prev));
    return results;
}

// Carga una matriz guardada como un CSV.
Matrix Matrix::load_from_file(std::string file_path){
    std::ifstream f(file_path);
    if(!f.good()){
        throw std::string("No se encontró el archivo...");
    }
    std::string line;
    unsigned short m=0,n=0;
    bool first_line = true;
    while(std::getline(f, line)){
        if(first_line){
            first_line = false;
            std::vector<std::string> tokens = split(line);
            n = tokens.size();
        }
        m++;
    }
    f.clear();
    f.seekg(0, std::ios::beg);
    Matrix result(m,n);
    unsigned short current_line_idx = 0;
    while(std::getline(f, line)){
        std::vector<std::string> tokens = split(line);
        if(tokens.size() != n) throw std::string("Error en la línea ")+std::to_string(current_line_idx+1);
        for(int i = 0; i < tokens.size(); i++){
            result.get(current_line_idx, i) = std::stod(tokens[i]);
        }
        current_line_idx++;
    }
    f.close();
    return result;
}

// Sobrecarga de suma entre dos matrices
Matrix operator+(Matrix a, Matrix b) {
    std::pair<unsigned short, unsigned short> a_shape = a.shape(), b_shape = b.shape();
    if(a_shape != b_shape) throw std::string("Las dimensiones no son correctas.");
    Matrix result(a_shape);
    for(int i = 0; i < a_shape.first; i++){
        for(int j = 0; j < a_shape.second; j++){
            result.get(i,j) = a.get(i,j)+b.get(i,j);
        }
    }
    return result;
}

// Sobrecarga de resta de dos matrices
Matrix operator-(Matrix a, Matrix b) {
    std::pair<unsigned short, unsigned short> a_shape = a.shape(), b_shape = b.shape();
    if(a_shape != b_shape) throw std::string("Las dimensiones no son correctas.");
    Matrix result = a;
    for(int i = 0; i < a_shape.first; i++){
        for(int j = 0; j < a_shape.second; j++){
            result.get(i,j) -= b.get(i,j);
        }
    }
    return result;
}

// Sobrecarga del operador de producto para dos matrices
Matrix operator*(Matrix a, Matrix b){
    std::pair<unsigned short, unsigned short> a_shape = a.shape(), b_shape = b.shape();
    if(a_shape.second != b_shape.first) throw std::string("Las dimensiones para el producto no son correctas.");
    Matrix result(a_shape.first, b_shape.second);
    for(int i = 0; i < a_shape.first; i++){
        for(int j = 0; j < b_shape.second; j++){
            for(int k = 0; k < a_shape.second; k++){
                result.get(i,j) += a.get(i,k)*b.get(k,j);
            }
        }
    }
    return result;
}

// Sobrecarga del operador de producto para una matriz (izquierda) y un escalar (derecha)
Matrix operator*(Matrix A, double v){
    std::pair<unsigned short, unsigned short> a_shape = A.shape();
    Matrix result(a_shape);
    for(int i = 0; i < a_shape.first; i++){
        for(int j = 0; j < a_shape.second; j++){
            result.get(i,j) = v*A.get(i,j);
        }
    }
    return result;
}

Matrix operator*(double v, Matrix A){
    return A*v;
}

#endif //INCLUDE_MATRIX