extern int POP[100][100],G,best_generation,best_individual[100] ,pop_new[100][100];
extern double best_fitness,fitness_value[100],fitness_table[100],fitness_avg[1000];

void GeneticAlgorithm(int pop_size, int chromo_size, int generation_size, double cross_rate, double mutate_rate, bool elitism );
