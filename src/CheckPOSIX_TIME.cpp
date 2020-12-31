
#include <unistd.h>

#if defined(_POSIX_VERSION) && defined(_POSIX_TIMERS) && (_POSIX_TIMERS > 0)
#include <ctime>  

#if defined(CLOCK_MONOTONIC)
#define Clk_ID CLOCK_MONOTONIC
#elif defined(CLOCK_REALTIME)
#define Clk_ID CLOCK_REALTIME
#elif defined(CLOCK_HIGHRES)
#define Clk_ID CLOCK_HIGHRES
#endif

#endif


#if (defined(Clk_ID))
// POSIX clock_gettime()

int retval;

double _ntl_GetWallTime( )
{
   using namespace std;
   timespec ts;
   if ((retval=clock_gettime(Clk_ID, &ts)))
      return -1;
   else
      return double(ts.tv_sec) + double(ts.tv_nsec) / 1000000000.0;
}
#endif

int main()
{
   double t = _ntl_GetWallTime();
   return retval;
}
