

#include <NTL/config.h>

#include <cstdlib>
#include <cstdio>

using namespace std;

unsigned long _ntl_GetPID();


int main()
{
   printf("%lu\n", _ntl_GetPID());
   return 0;
}
