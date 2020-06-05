#ifndef INCLUDE_SORT_MACHINES_H
#define INCLUDE_SORT_MACHINES_H

#include<vector>
#include<iostream>
#include<algorithm>
#include<random> //To get the mersenne twister engine
#include<chrono> //To get time

////////////////////
// ESTO IGNORENLO //
////////////////////
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
/////////////////////////////////
// DEJEN DE IGNORAR DESDE AQUI //
/////////////////////////////////
class SortMachine{
public:
    SortMachine(){};

    virtual void sort(std::vector<int> &target) = 0;

    // Imprime el vector target de una manera bonita
    static void print(std::vector<int> &target){
        std::cout << target[0];
        for(int i = 1; i < target.size(); i++){
            std::cout << ", " << target[i];
        }
        std::cout << "\n";
    }

    // Hace aleatorio el vector target
    void shuffle(std::vector<int> &target){
        std::shuffle(target.begin(), target.end(), RandGen::get().gen);
    }
    std::string name = "No Name Sort";
};

/************************************** *
 * A partir de aquí comienza su trabajo *
 * *************************************/

class BubbleSort: public SortMachine{
public:
    BubbleSort(){this->name = "Bubble Sort";}
    void sort(std::vector<int> &target){
        for (int i = 0; i < target.size()-1; i++)
        {
            for (int j = 0; j < target.size()-1-i; j++)
            {
                if (target[j] > target[j+1])
                {
                    int x = target[j];
                    target[j] = target[j+1];
                    target[j+1] = x;
                }
            }
        }
        

    }
};
class InsertionSort: public SortMachine{
public:
    InsertionSort(){this->name = "Insertion Sort";}
    void sort(std::vector<int> &target){
        int i = 1;
        while (i < target.size())
        {
            int x = target[i];
            int j = i-1;

            while (j >= 0 && target[j] > x)
            {
                target[j+1] = target[j];
                j = j-1;
            }
            target[j+1] = x;
            i ++;  
        }   
    }
};
class MergeSort: public SortMachine{
public:
    MergeSort(){this->name = "Merge Sort";}
    void Merge(std::vector<int> &target, int p, int q, int r){
        int n1 = q-p+1, n2 = r-q;
        int i, j, k;
        std::vector<int> left(n1), right(n2);

        for (i = 0; i < n1; i++)
        {
            left[i] = target[p+i];
        }
        for (j = 0; j < n2; j++)
        {
            right[j] = target[q+j+1];
        }
        i = 0; j = 0;

        for (k = p; i < n1 && j < n2; k++)
        {
            if (left[i] < right[j])
            {
                target[k] = left[i++];
            }
            else
            {
                target[k] = right[j++];
            }
        }    

        while (i < n1)
        {
            target[k++] = left[i++];
        }
        while (j < n2)
        {
            target[k++] = right[j++];
        }              
    }
    void MergeSort1(std::vector<int> &target, int p, int r){
        int q;
        if (p < r)
        {
            q = (p+r)/2;
            MergeSort1(target, p, q);
            MergeSort1(target, q+1, r);
            Merge(target, p, q ,r);
        }
        
    }

    void sort(std::vector<int> &target){
        MergeSort1(target, 0, target.size());
    }
};
class HeapSort: public SortMachine{
public:
    HeapSort(){this->name = "Heap Sort";}
    void heapify(std::vector<int> &target, int n, int i){
        int largest =  i;
        int left = 2*i+1;
        int right = 2*i+2;

        if (left < n && target[left] > target[largest])
        {
            largest = left;
        }

        if (right < n && target[right] > target[largest])
        {
            largest = right;
        }

        if (largest != i)
        {
            int x = target[i];
            target[i] = target[largest];
            target[largest] = x;

            heapify(target, n, largest);
        }
        
    }
    void sort(std::vector<int> &target){
        int n = target.size();
        for (int i = n/2 - 1; i >= 0; i--)
        {
            heapify(target,n,i);
        }

        for (int i = n-1; i >= 0; i--)
        {
            int x = target[0];
            target[0] = target[i];
            target[i] = x;

            heapify(target, i, 0);
        }   
    }
};

#endif //INCLUDE_SORT_MACHINES_H