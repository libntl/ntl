#include <NTL/ctools.h>

#include <cstdlib>
#include <iostream>

#include <NTL/simde_pclmul.h>



#if (NTL_BITS_PER_LONG != 64)
#error "PCLMUL not supported"
#endif

void
mul1 (unsigned long *c, unsigned long a, unsigned long b)
{
   pclmul_mul1(c, a, b);
}
   


using namespace std;


int main()
{

// Example with 8 bits per word
//        11110000
//           x 110
//----------------
//       111100000
//      1111000000
//----------------
//      1000100000

   unsigned long a = ((unsigned long) _ntl_nofold(15)) << (NTL_BITS_PER_LONG-4);
   unsigned long b = _ntl_nofold(6);
   unsigned long c[2];

   mul1(c, a, b);

   unsigned long c0 = ((unsigned long) _ntl_nofold(1)) << (NTL_BITS_PER_LONG-3);
   unsigned long c1 = _ntl_nofold(2);


   if (c[0] == c0 && c[1] == c1) 
      return 0;
   else
      return -1;
}
