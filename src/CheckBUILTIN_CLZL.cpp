#include <NTL/mach_desc.h>
#include <cstdlib>

using namespace std;


long CountLeadingZeros(unsigned long x)
{
   return __builtin_clzl(x);
}

int main()
{
   unsigned long x = atoi("3");
   if (CountLeadingZeros(x) == NTL_BITS_PER_LONG-2)
      return 0;
   else
      return -1;
}



