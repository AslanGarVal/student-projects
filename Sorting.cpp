#include "SortMachines.h"

int main() {
    // Iniciamos la máquina, con un BubbleSort
    SortMachine *my_sort_machine = new HeapSort();

    std::vector<int> x(10);
    for(int i = 0; i < x.size(); i++){
        x[i] = i;
    }
    // Alteramos el orden de x
    my_sort_machine->shuffle(x);
    
    // Lo imprimimos nadamás porque sí
    SortMachine::print(x);

    // En principio, ordenamos
    my_sort_machine->sort(x);
    
    // Imprimimos el resultado
    SortMachine::print(x);
    return 0;
}