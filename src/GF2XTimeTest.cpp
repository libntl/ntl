
#include <NTL/GF2X.h>

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


#ifdef NTL_GF2X_ALTCODE
printf("NTL_GF2X_ALTCODE ");
#endif


#ifdef NTL_GF2X_ALTCODE1
printf("NTL_GF2X_ALTCODE1 ");
#endif


#ifdef NTL_GF2X_NOINLINE
printf("NTL_GF2X_NOINLINE ");
#endif


printf("\n");

}

int main()
{
   long n, i, j, iter, s, k;
   double t;

   SetSeed(ZZ(0));


   for (i = 0; i < 10000; i++) {
      GF2X a, b, c, d;
      long da = RandomBnd(5*NTL_BITS_PER_LONG);
      long db = RandomBnd(5*NTL_BITS_PER_LONG);
      long dc = RandomBnd(5*NTL_BITS_PER_LONG);
      long dd = RandomBnd(5*NTL_BITS_PER_LONG);
      random(a, da);  random(b, db);  random(c, dc);  random(d, dd);

      if ((a + b)*(c + d) != c*a + d*a + c*b + d*b) {
	 printf("999999999999999 ");
	 print_flag();
	 return 0;
      }
   }
   

   n = 16;
   s = 56;

   GF2X *a = new GF2X[s];
   GF2X *b = new GF2X[s];

   GF2X c;

   for (k = 0; k < s; k++) {
      random(a[k], (n + (k % 7))*NTL_BITS_PER_LONG);
      random(b[k], (n + (k % 8))*NTL_BITS_PER_LONG);
   }

   for (k = 0; k < s; k++) mul(c, a[k], b[k]);


   iter = 1;

   do {
     t = GetTime();
     for (i = 0; i < iter; i++) {
        for (j = 0; j < 1; j++) for (k = 0; k < s; k++) mul(c, a[k], b[k]);
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
     for (i = 0; i < iter; i++) {
        for (j = 0; j < 1; j++) for (k = 0; k < s; k++) mul(c, a[k], b[k]);
     }
     t = GetTime() - t;
     tvec[w] = t;
   }


   t = clean_data(tvec);

   t = floor((t/iter)*1e14);

   if (t < 0 || t >= 1e15)
      printf("999999999999999 ");
   else
      printf("%015.0f ", t);

   printf(" [%ld] ", iter);

   print_flag();

   return 0;
}
   

