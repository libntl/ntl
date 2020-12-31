
#if defined(__MACH__) && defined(__APPLE__)
#include <mach/mach.h>
#include <mach/mach_time.h>

static inline double InitTimeConvert()
{
   mach_timebase_info_data_t timeBase;
   (void)mach_timebase_info( &timeBase );
   return (double)timeBase.numer / (double)timeBase.denom / 1000000000.0;
}

double _ntl_GetWallTime( )
{
   static double timeConvert = InitTimeConvert();
   // even in a multi-threaded environment, this will
   // be safely initialized, according to C++11 standard

   return double(mach_absolute_time()) * timeConvert;
}
#endif


int main()
{
   double t = _ntl_GetWallTime( );
   return 0;
}
