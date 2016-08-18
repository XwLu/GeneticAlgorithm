#include <iostream>
#include "function.h"

using namespace std;

int POP[100][100],G,best_generation=0,best_individual[100],pop_new[100][100] ;
double best_fitness,fitness_value[100],fitness_table[100],fitness_avg[1000];

int main()
{
    int pop_size, chromo_size,  generation_size;
    double cross_rate, mutate_rate;
    bool elitism;

    elitism = true;//ѡ��Ӣ����
    pop_size = 20;//��Ⱥ��С
    chromo_size = 16;//Ⱦɫ���С
    generation_size = 1000;//��������
    cross_rate = 0.9;// �������
    mutate_rate = 0.2;//�������

    GeneticAlgorithm(pop_size, chromo_size, generation_size, cross_rate, mutate_rate, elitism );

    return 0;
}
