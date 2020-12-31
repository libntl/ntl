
void touch_double(double* x);
double val_double(double x);
double power2(long k);
double fma_test(double a, double b, double c);


long FMADetected(long dp)
{
   double x = power2(0) + power2(dp-1);
   double y = power2(0) + power2(dp-1);

   touch_double(&x);
   touch_double(&y);

   double z = x*y;
   touch_double(&z);
   z = -z;
   touch_double(&z);

   double lo = fma_test(x, y, z);
   return lo != 0;
}

long DoublePrecision()
{
   double eps, one, res;
   long k;

   one = val_double(1.0);
   eps = val_double(1.0);

   k = 0;

   do {
      double tmp;

      k++;
      eps *= 1.0/2.0;
      tmp = 1.0 + eps;
      touch_double(&tmp);
      res = tmp - one;
   } while (res == eps);

   return k;
}


int main()
{
   long dp = DoublePrecision();
   long fma = FMADetected(dp);

   return int(fma);
}
