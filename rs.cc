// rs.cc
// clang++ -std=c++11 rs.c -o rs or g++ -std=c++11 rs.c -o rs
//_RUNNING_STATISTICS_MT_: clang++ -std=c++11 rs.cc -o rs -lpthread or g++ -std=c++11 rs.cc -o rs -lpthread
#include <iostream>
#include <iomanip>
#include <vector>
#include <random>
#include <tuple>

#include "RunningStatistics.hh"

#ifdef _RUNNING_STATISTICS_MT_
template <typename Int, typename Real>
std::shared_timed_mutex RunningStatistics<Int, Real>::mutex;
#endif

void run(RunningStatistics<long, double>& instance)
{
   std::random_device rd{};
   std::mt19937 gen{rd()};
   std::normal_distribution<> d{0,5};
   for(int n=0; n<1000000; ++n)
   {
      instance(d(gen));
      if (n && (n % 1000) == 0)
      {
         long nn;
         double mean, var, skew, kurt;
         std::tie(nn, mean, var, skew, kurt) = instance();
         double dev = sqrt(var);
         std::cout << std::fixed << std::setprecision(4) << n << ": " << mean  << ", " << dev << std::endl;
      }
   }
}

int main(int argc, char **argv)
{
   RunningStatistics<long, double> instance;
#if defined(_RUNNING_STATISTICS_MT_) || defined(_RUNNING_STATISTICS_MT_LOCKFREE_)
   std::vector<std::thread> threads;
   for (unsigned n=0; n<std::thread::hardware_concurrency()+1; n++)
      threads.emplace_back(run, std::ref(instance));
   for (std::thread& t : threads) t.join();
#else
   run(instance);
#endif
}

