
#include <omp.h>
#include <vector> 
#include <iostream>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */

using namespace std;

int main()
{

vector<int> vec;
for (int i = 0; i < 100000000; i++)
{
 vec.push_back(i);
}
for (auto i = vec.begin(); i < vec.end(); i++)
{
  (*i)++;

}

clock_t timer_t;
timer_t = clock();

omp_set_num_threads(4);
//#pragma omp parallel
//{
#pragma omp parallel for
for (auto i = vec.begin(); i < vec.end(); i++)
{
//#pragma omp single nowait
//{
 (*i)++;
//}
}

//}

timer_t = clock() - timer_t;
printf ("It took me %d clicks (%f seconds).\n",timer_t,((float)timer_t)/CLOCKS_PER_SEC);



vector<int> vec1;
clock_t timer_t1;
timer_t1 = clock();

for (auto i = vec.begin(); i < vec.end(); i++)
{
  (*i)++;

}

timer_t1 = clock() - timer_t1;
printf ("It took me %d clicks (%f seconds).\n",timer_t1,((float)timer_t1)/CLOCKS_PER_SEC);


return 0;
}
