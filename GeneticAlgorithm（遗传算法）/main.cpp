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

    elitism = true;//选择精英操作
    pop_size = 20;//种群大小
    chromo_size = 16;//染色体大小
    generation_size = 1000;//迭代次数
    cross_rate = 0.9;// 交叉概率
    mutate_rate = 0.2;//变异概率

    GeneticAlgorithm(pop_size, chromo_size, generation_size, cross_rate, mutate_rate, elitism );

    return 0;
}
