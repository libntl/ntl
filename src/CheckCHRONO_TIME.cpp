#include <NTL/ctools.h>

#if (NTL_CXX_STANDARD >= 2011)

#include <chrono>

double _ntl_GetWallTime( )
{
   auto current_time = std::chrono::steady_clock::now();
   auto duration_in_seconds = std::chrono::duration<double>(current_time.time_since_epoch());
   return duration_in_seconds.count();
}

#endif


int main()
{
   double t = _ntl_GetWallTime( );
   return 0;
}
