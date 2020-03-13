/**
A templatized single header implementation of Welfords algorithm for online
mean and variance.
Adapted for modern C++ from https://www.johndcook.com/blog/skewness_kurtosis/
which is based on Knuth Volume 2 (Seminumerical Algorithms) which in turn
describes Welford's [Welford (1962).Note on a method for calculating corrected
sums of squares and products] algorithm.

Requires C++14 for multithreaded version which uses reader-writer locks.

    Copyright (c) 2018: Donald Munro.
    Permission is hereby granted, free of charge, to any person obtaining a copy of this software
    and associated documentation files (the "Software"), to deal in the Software without restriction,
    including without limitation the rights to use, copy, modify, merge, publish, distribute,
    sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    The above copyright notice and this permission notice shall be included in all copies or
    substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef _RUNNING_STATISTICS_H_
#define _RUNNING_STATISTICS_H_

/*
 * WARNING: Its often better to prefer per thread instances and
 * combine them with operator+ or += or use fibers and a single non-MT instance.
 */
//#define _RUNNING_STATISTICS_MT_

#include <cmath>
#ifdef _RUNNING_STATISTICS_MT_
#include <thread>
#include <mutex>
#include <shared_mutex> //reader-writer locks
#endif
#include <tuple>

template <typename Int =long, typename Real=double>
class RunningStatistics
//=====================
{
private:
   Int n;
   Real M1, M2, M3, M4;
#ifdef _RUNNING_STATISTICS_MT_
   static std::shared_timed_mutex mutex;
#endif

public:
   RunningStatistics() : n(0), M1(0), M2(0), M3(0), M4(0) {}

   void clear() { n = 0; M1 = M2 = M3 = M4 = static_cast<Real>(0); }

   void operator()(const Real x)
   //-----------------------------
   {
#ifdef _RUNNING_STATISTICS_MT_
      std::unique_lock<std::shared_timed_mutex> lock(mutex);
#endif
      Real delta, delta_n, delta_n2, term1;
      Int n1 = n;
      n++;
      delta = x - M1;
      delta_n = delta / n;
      delta_n2 = delta_n * delta_n;
      term1 = delta * delta_n * n1;
      M1 += delta_n;
      M4 += term1 * delta_n2 * (n*n - 3*n + 3) + 6 * delta_n2 * M2 - 4 * delta_n * M3;
      M3 += term1 * delta_n * (n - 2) - 3 * delta_n * M2;
      M2 += term1;
   }

   friend RunningStatistics operator+(const RunningStatistics a, const RunningStatistics b)
   //--------------------------------------------------------------------------------------
   {
#ifdef _RUNNING_STATISTICS_MT_
      std::unique_lock<std::shared_timed_mutex> lock(mutex);
#endif
      RunningStatistics combined;
      combined.n = a.n + b.n;

      Real delta = b.M1 - a.M1;
      Real delta2 = delta*delta;
      Real delta3 = delta*delta2;
      Real delta4 = delta2*delta2;

      combined.M1 = (a.n*a.M1 + b.n*b.M1) / combined.n;

      combined.M2 = a.M2 + b.M2 +
                    delta2 * a.n * b.n / combined.n;

      combined.M3 = a.M3 + b.M3 +
                    delta3 * a.n * b.n * (a.n - b.n)/(combined.n*combined.n);
      combined.M3 += 3.0*delta * (a.n*b.M2 - b.n*a.M2) / combined.n;

      combined.M4 = a.M4 + b.M4 + delta4*a.n*b.n * (a.n*a.n - a.n*b.n + b.n*b.n) /
                                  (combined.n*combined.n*combined.n);
      combined.M4 += 6.0*delta2 * (a.n*a.n*b.M2 + b.n*b.n*a.M2)/(combined.n*combined.n) +
                     4.0*delta*(a.n*b.M3 - b.n*a.M3) / combined.n;

      return combined;
   }

   RunningStatistics& operator+=(const RunningStatistics& rhs)
   //---------------------------------------------------------
   {
#ifdef _RUNNING_STATISTICS_MT_
      std::unique_lock<std::shared_timed_mutex> lock(mutex);
#endif
      RunningStatistics combined = *this + rhs;
      *this = combined;
      return *this;
   }

   // return count, mean, variance, skew, kurtosis
   std::tuple<Int, Real, Real, Real, Real> operator()() const
   //--------------------------------------------------
   {
#ifdef _RUNNING_STATISTICS_MT_
      std::shared_lock<std::shared_timed_mutex> lock(mutex);
#endif
      return std::make_tuple(n, M1, M2/(n-1.0), sqrt(double(n)) * M3/ pow(M2, 1.5),
                             double(n)*M4 / (M2*M2) - 3.0);
   }

   // return count, mean, variance
   std::tuple<Int, Real, Real> mean_var()
   //--------------------------------------------------
   {
#ifdef _RUNNING_STATISTICS_MT_
      std::shared_lock<std::shared_timed_mutex> lock(mutex);
#endif
      return std::make_tuple(n, M1, M2/(n-1.0));
   }

   Int size() const
   //--------------
   {
#ifdef _RUNNING_STATISTICS_MT_
      std::shared_lock<std::shared_timed_mutex> lock(mutex);
#endif
      return n;
   }

   Real mean() const
   //---------------
   {
#ifdef _RUNNING_STATISTICS_MT_
      std::shared_lock<std::shared_timed_mutex> lock(mutex);
#endif
      return M1;
   }

   Real variance() const
   //-------------------
   {
#ifdef _RUNNING_STATISTICS_MT_
      std::shared_lock<std::shared_timed_mutex> lock(mutex);
#endif
      return M2/(n-1.0);
   }

   Real deviation() const
   //---------------------
   {
#ifdef _RUNNING_STATISTICS_MT_
      std::shared_lock<std::shared_timed_mutex> lock(mutex);
#endif
      return sqrt(variance());
   }

   Real skewness() const
   //-------------------
   {
#ifdef _RUNNING_STATISTICS_MT_
      std::shared_lock<std::shared_timed_mutex> lock(mutex);
#endif
      return sqrt(double(n)) * M3/ pow(M2, 1.5);
   }

   Real kurtosis() const
   //-------------------
   {
#ifdef _RUNNING_STATISTICS_MT_
      std::shared_lock<std::shared_timed_mutex> lock(mutex);
#endif
      return double(n)*M4 / (M2*M2) - 3.0;
   }
};
#endif // _RUNNING_STATISTICS_H_
