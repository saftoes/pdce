// (C) 2018 University of Bristol. See License.txt


#include "Tools/time-func.h"
#include "Exceptions/Exceptions.h"


long long timeval_diff(struct timeval *start_time, struct timeval *end_time)
{ struct timeval temp_diff;
  struct timeval *difference;
  difference=&temp_diff;
  difference->tv_sec =end_time->tv_sec -start_time->tv_sec ;
  difference->tv_usec=end_time->tv_usec-start_time->tv_usec;
  while(difference->tv_usec<0)
  { difference->tv_usec+=1000000;
    difference->tv_sec -=1;
  }
  return 1000000LL*difference->tv_sec+difference->tv_usec;
}

double timeval_diff_in_seconds(struct timeval *start_time, struct timeval *end_time)
{
    return double(timeval_diff(start_time, end_time)) / 1e6;
}


long long timespec_diff(struct timespec *start_time, struct timespec *end_time)
{
  long long sec =end_time->tv_sec -start_time->tv_sec ;
  long long nsec=end_time->tv_nsec-start_time->tv_nsec;
  while(nsec<0)
  { nsec+=1000000000;
    sec -=1;
  }
  return 1000000000*sec+nsec;
}


double convert_ns_to_seconds(long long x)
{
  return double(x) / 1e9;
}


double Timer::elapsed()
{
  long long res = elapsed_time;
  if (running)
    res += elapsed_since_last_start();
  return convert_ns_to_seconds(res);
}

double Timer::idle()
{
  if (running)
    throw Processor_Error("Timer running.");
  else
    return convert_ns_to_seconds(elapsed_since_last_start());
}
