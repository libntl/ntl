
#include <NTL/ZZ_pX.h>
#include <NTL/FFT.h>

#include <cstdio>

NTL_CLIENT


double clean_data(double *t)
{
   double x, y, z;
   long i, ix, iy, n;

   x = t[0]; ix = 0;
   y = t[0]; iy = 0;

   for (i = 1; i < 5; i++) {
      if (t[i] < x) {
         x = t[i];
         ix = i;
      }
      if (t[i] > y) {
         y = t[i];
         iy = i;
      }
   }

   z = 0; n = 0;
   for (i = 0; i < 5; i++) {
      if (i != ix && i != iy) z+= t[i], n++;
   }

   z = z/n;  

   return z;
}


void print_flag()
{

#if defined(NTL_FFT_LAZYMUL)
printf("FFT_LAZYMUL ");
#endif

#if defined(NTL_SPMM_ULL)
printf("SPMM_ULL ");
#endif

#if defined(NTL_AVOID_BRANCHING)
printf("AVOID_BRANCHING ");
#endif

#if defined(NTL_FFT_BIGTAB)
printf("FFT_BIGTAB ");
#endif

printf("\n");

}


int main()
{

   SetSeed(ZZ(0));


   long n, k;

   n = 200;
   k = 10*NTL_ZZ_NBITS;

   ZZ p;

   RandomLen(p, k);
   if (!IsOdd(p)) p++;


   ZZ_p::init(p);         // initialization

   ZZ_pX f, g, h, r1, r2, r3;

   random(g, n);    // g = random polynomial of degree < n
   random(h, n);    // h =             "   "
   random(f, n);    // f =             "   "

   SetCoeff(f, n);  // Sets coefficient of X^n to 1
   

   // For doing arithmetic mod f quickly, one must pre-compute
   // some information.

   ZZ_pXModulus F;
   build(F, f);

   PlainMul(r1, g, h);  // this uses classical arithmetic
   PlainRem(r1, r1, f);

   MulMod(r2, g, h, F);  // this uses the FFT

   MulMod(r3, g, h, f);  // uses FFT, but slower

   // compare the results...

   if (r1 != r2) {
      printf("999999999999999 ");
      print_flag();
      return 0;
   }
   else if (r1 != r3) {
      printf("999999999999999 ");
      print_flag();
      return 0;
   }

   double t;
   long i, j;
   long iter;

   const int nprimes = 30;
   const long L = 12; 
   const long N = 1L << L;
   long r;
   

   for (r = 0; r < nprimes; r++) UseFFTPrime(r);

   vec_long A1[nprimes], A2[nprimes];
   vec_long B1[nprimes], B2[nprimes];

   for (r = 0; r < nprimes; r++) {
      A1[r].SetLength(N);
      A2[r].SetLength(N);
      B1[r].SetLength(N);
      B2[r].SetLength(N);

      for (i = 0; i < N; i++) {
         A1[r][i] = RandomBnd(GetFFTPrime(r));
         A2[r][i] = RandomBnd(GetFFTPrime(r));
      }
   }

   for (r = 0; r < nprimes; r++) {
      long *A1p = A1[r].elts();
      long *A2p = A2[r].elts();
      long *B1p = B1[r].elts();
      long *B2p = B2[r].elts();
      long q = GetFFTPrime(r);
      mulmod_t qinv = GetFFTPrimeInv(r);
 
      FFTFwd(B1p, A1p, L, r);
      FFTFwd(B2p, A2p, L, r);
      for (i = 0; i < N; i++) B1p[i] = NormalizedMulMod(B1p[i], B2p[i], q, qinv);
      FFTRev1(B1p, B1p, L, r);
   }

   iter = 1;

   do {
     t = GetTime();
     for (j = 0; j < iter; j++) {
        for (r = 0; r < nprimes; r++) {
           long *A1p = A1[r].elts();
           long *A2p = A2[r].elts();
           long *B1p = B1[r].elts();
           long *B2p = B2[r].elts();
           long q = GetFFTPrime(r);
           mulmod_t qinv = GetFFTPrimeInv(r);

           FFTFwd(B1p, A1p, L, r);
           FFTFwd(B2p, A2p, L, r);
           for (i = 0; i < N; i++) B1p[i] = NormalizedMulMod(B1p[i], B2p[i], q, qinv);
           FFTRev1(B1p, B1p, L, r);
        }
     }
     t = GetTime() - t;
     iter = 2*iter;
   } while(t < 1);

   iter = iter/2;

   iter = long((3/t)*iter) + 1;


   double tvec[5];
   long w;

   for (w = 0; w < 5; w++) {
     t = GetTime();
     for (j = 0; j < iter; j++) {
        for (r = 0; r < nprimes; r++) {
           long *A1p = A1[r].elts();
           long *A2p = A2[r].elts();
           long *B1p = B1[r].elts();
           long *B2p = B2[r].elts();
           long q = GetFFTPrime(r);
           mulmod_t qinv = GetFFTPrimeInv(r);

           FFTFwd(B1p, A1p, L, r);
           FFTFwd(B2p, A2p, L, r);
           for (i = 0; i < N; i++) B1p[i] = NormalizedMulMod(B1p[i], B2p[i], q, qinv);
           FFTRev1(B1p, B1p, L, r);
        }
     }
     t = GetTime() - t;
     tvec[w] = t;
   }

   t = clean_data(tvec);

   t = floor((t/iter)*1e13);

   if (t < 0 || t >= 1e15)
      printf("999999999999999 ");
   else
      printf("%015.0f ", t);

   printf(" [%ld] ", iter);

   print_flag();

   return 0;
}
