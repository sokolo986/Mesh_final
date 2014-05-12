
#include <omp.h>
#include <vector> 
#include <iostream>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */

using namespace std;

int main()
{
int sz = 100000000;
vector<int> vec;
vec.reserve(sz);

for (int i = 0; i < sz; i++)
{
 vec.push_back(i);
}

for (auto i = vec.begin(); i < vec.end(); i++)
{
  (*i)++;

}

clock_t timer_t;
for (int j=0; j <= 9; ++j){
	
	
	
	//omp_set_num_threads(j);
	timer_t = clock();
	//#pragma omp parallel
	//{
	#pragma omp parallel for num_threads(j)
		
		for (auto i = vec.begin(); i < vec.end(); i++)
		{
			//if (i==vec.begin())
			//	printf ("J= %d and thread %d\n",j,omp_get_thread_num());
			//if (40000000==vec.end()-i)
			//	printf ("J= %d and thread %d\n\n",j,omp_get_thread_num());
		//#pragma omp single nowait
		//{
		 (*i)++;
		
		//}
		}
	
	timer_t = clock() - timer_t;
	printf ("%d threads took %ld clicks (%f seconds).\n",j,timer_t,((float)timer_t)/CLOCKS_PER_SEC);
}


vector<int> vec1;
clock_t timer_t1;
timer_t1 = clock();

for (auto i = vec.begin(); i < vec.end(); i++)
{
  (*i)++;
}

timer_t1 = clock() - timer_t1;
printf ("Standard time took me %ld clicks (%f seconds).\n",timer_t1,((float)timer_t1)/CLOCKS_PER_SEC);
vec.clear();

return 0;
}
