#ifdef NTL_FP_CONTRACT_OFF
#pragma fp_contract(off)
#endif

void touch_double(double* x) {}

double val_double(double x) { return x; }

double power2(long k)
{
   long i;
   double res;

   res = 1;

   for (i = 1; i <= k; i++)
      res = res * 2;

   return res;
}

double fma_test(double a, double b, double c)
{
   double t1 = a*b;
   double t2 = t1 + c;
   return t2;
}


