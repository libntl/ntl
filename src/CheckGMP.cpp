
#include <NTL/config.h>


#ifdef NTL_GMP_LIP
#include <gmp.h>
#endif


int main() { 
#ifdef NTL_GMP_LIP
   mpz_t pie;
   mpz_init_set_str (pie, "3141592653589793238462643383279502884", 10);
#endif

   return 0; 
}







