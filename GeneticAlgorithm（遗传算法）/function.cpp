#include "function.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

using namespace std;
double get_rand(int range)
{
    static int i=1;

    if(i)
    {
       srand((unsigned)time(NULL));    //��ʱ�����ֲ��������
       i=0;
    }
    return rand()%range;
}
/*�������뺯��*/
int myround(double num)
{
    double x;
    x=num-(int)num;
    if (x>=0.5) return (int)num+1;
    else return (int)num;
}
/*��ʼ����Ⱥ
pop_size����Ⱥ��С
chromo_size��Ⱦɫ�峤��*/
void initilize(int pop_size,int chromo_size)
{
  int i,j;
  double x;
  int num;
  for(i=0;i<pop_size;i++)
  {
   for(j=0;j<chromo_size;j++)
    {
    num = get_rand(10000);
    x=num/10000.0;
    POP[i][j]=myround(x);
    }
  }
}
/*��Ӧ�ȼ��㺯��*/
void fitness(int pop_size, int chromo_size)
{
    int i,j;
    for(i=0;i<pop_size;i++)
        fitness_value[i]=0;
    for(i=0;i<pop_size;i++)
      {
        for(j=0;j<chromo_size;j++)
        {
         if(POP[i][j] == 1)
          fitness_value[i] = fitness_value[i] + pow(2.0,(double)j);
        }
      }
    for(i=0;i<pop_size;i++)
    {
            fitness_value[i] = -1+fitness_value[i] * (3-(-1))/(pow(2.0,(double)chromo_size-1));//����δ֪����ȡֵ
            fitness_value[i] = 6-pow(fitness_value[i]-2,2.0);      //���ݱ��ʽ��ļ�ֵ��
    }

}
/*������*/
void myrank(int pop_size, int chromo_size)
{
    int i,j,k,mymin,temp1[100];
    double mytemp;
    for(i=0;i<pop_size;i++)
        fitness_table[i] = 0;
    mymin=0;
    mytemp=0;
    temp1[chromo_size-1]=0;
    for (i=0;i<pop_size;i++)
    {
        mymin = i;
        for(j = i+1;j<pop_size;j++)
        {
            if (fitness_value[j]<fitness_value[mymin])
            mymin = j;
        }

       if(mymin!=i)
        {
        mytemp = fitness_value[i];
        fitness_value[i] = fitness_value[mymin];
        fitness_value[mymin] = mytemp;
        for(k=0;k<chromo_size;k++)
        {
           temp1[k] = POP[i][k];
           POP[i][j] = POP[mymin][k];
           POP[mymin][k] = temp1[k];
        }
        }
    }
    for(i=0;i<pop_size;i++)
    {
        if(i==0)
            fitness_table[i] = fitness_table[i] + fitness_value[i];
        else
            fitness_table[i] = fitness_table[i-1] + fitness_value[i];
    }
    fitness_avg[G] = fitness_table[pop_size-1]/(double)pop_size;

    if(fitness_value[pop_size-1]>best_fitness)
    {
           best_fitness = fitness_value[pop_size-1];
           best_generation = G;
           for (j=0;j<chromo_size;j++)
          best_individual[j] = POP[pop_size-1][j];
    }
}

void selection(int pop_size, int chromo_size, bool elitism)
{
    double r,x;
    int i,j,num,first,last,mid,idx,p;
    for(i=0;i<pop_size;i++)
    {
        num=get_rand(10000);
        x=num/10000.0;
        r = x * fitness_table[pop_size-1];
        first = 0;
        last = pop_size-1;
        mid = myround((double)(first+last)/2.0);
        idx = -1;
        while (first <= last&&idx == -1)
        {
            if(r > fitness_table[mid])
                first = mid;
            else if (r < fitness_table[mid])
                last = mid;
            else {idx = mid;break;}

            mid = myround((double)(last+first)/2.0);
            if(last - first == 1)
                idx = last;
        }
        for (j=0;j<chromo_size;j++)
            pop_new[i][j] = POP[idx][j];
    }
         if (elitism)
            p = pop_size - 1;
         else
            p = pop_size;
         for (i=0;i<p;i++)
            for(j=0;j<chromo_size;j++)
            POP[i][j] = pop_new[i][j];
}
/*���亯��*/
void crossover(int pop_size,int chromo_size, double cross_rate)
{
    int num,cross_pos,j,i,temp;
    double x,y;
    for(i=0;i<pop_size;i=i+2)
    {
        num=get_rand(10000);
        x=num/10000.0;
        if( x < cross_rate)
        {
            num = get_rand(10000);
            y = num/10000.0;
            cross_pos = myround(y * chromo_size);//�����������λ��
            if(cross_pos == 0)
                continue;
            for(j=cross_pos;j<chromo_size;j++)
            {
             temp = POP[i][j];
             POP[i][j] = POP[i+1][j];
             POP[i+1][j] = temp;
            }
        }
    }
}
/*���캯��*/
void mutation(int pop_size, int chromo_size, double mutate_rate)
{
    int i,num,mutate_pos;
    double x,y;
    for(i = 0; i<pop_size;i++)
    {
        num=get_rand(10000);
        x=num/10000.0;
        if(x < mutate_rate)
        {
            num = get_rand(10000);
            y = num/10000.0;
            mutate_pos = myround(y * chromo_size);
            POP[i][mutate_pos] = 1 - POP[i][mutate_pos];
        }
    }
}
/*�Ŵ����̺���*/
void GeneticAlgorithm(int pop_size, int chromo_size, int generation_size, double cross_rate, double mutate_rate, bool elitism )
{
    int i,j,*m,p;
    double n,q;
    for(i=0;i<generation_size;i++)
    fitness_avg[i] = 0;

    fitness_value[pop_size-1] = 0;
    best_fitness = 0;
    best_generation = 0;
    initilize(pop_size,chromo_size); //��ʼ��
    for (G = 1; G<=generation_size;G++)
    {
        fitness(pop_size, chromo_size);//��Ӧ�ȼ���
        myrank(pop_size,chromo_size);//�Ը��尴��Ӧ�ȴ�С����
        selection(pop_size,chromo_size,elitism);//ѡ�����
        crossover(pop_size,chromo_size,cross_rate);//�������
        mutation(pop_size,chromo_size,mutate_rate);//�������
    }
    m = best_individual; //�����Ѹ���
    n = best_fitness;// ��������Ӧ��
    p = best_generation;//�����Ѹ�����ִ�

    /*�Բ�ͬ�ĺ�����һ����Ҫ��д*/
    q=0;
    for(j=0;j<chromo_size;j++)
        if (best_individual[j] == 1)
        q = q + pow(2.0,(double)j);
    q = -1+q * (3.0-(-1.0))/pow(2.0,chromo_size-1); //δ֪���Ľ��ֵ

    cout<<"�����Ѹ���";
    for(i=0;i<chromo_size;i++)
    cout<<*(m+i);
    cout<<endl;
    cout<<"��������Ӧ��"<<n<<endl;
    cout<<"�����Ѹ�����ִ�"<<p<<endl;
    cout<<"������ֵ"<<q<<endl;

}
