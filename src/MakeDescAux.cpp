
#include <cstdlib>
using namespace std;

int val_int(int x) 
{ volatile int _x = x; return _x; }

unsigned int val_uint(unsigned int x) 
{ volatile unsigned int _x = x;  return _x; }
 
long val_long(long x) 
{ volatile long _x = x; return _x; }

unsigned long val_ulong(unsigned long x) 
{ volatile unsigned long _x = x; return _x; }
 
size_t val_size_t(size_t x) 
{  volatile size_t _x = x; return _x; }

double val_double(double x) 
{  volatile double _x = x; return _x; }

long double val_ldouble(double x) 
{ volatile double _x = x; return _x; }



 
void touch_int(int* x) 
{ *x = val_int(*x); }

void touch_uint(unsigned int* x) 
{ *x = val_uint(*x); }
 
void touch_long(long* x) 
{ *x = val_long(*x); }

void touch_ulong(unsigned long* x) 
{ *x = val_ulong(*x); }

void touch_size_t(size_t* x) 
{ *x = val_size_t(*x); }
 
void touch_double(double* x) 
{ *x = val_double(*x); }

void touch_ldouble(long double* x) 
{ *x = val_ldouble(*x); }





double fma_test(double a, double b, double c)
{
   double t1 = a*b;
   double t2 = t1 + c;
   return t2;
}

double reassoc_test(double a, double b, double c, double d)
{
   double t1 = a*c + a*d;
   double t2 = b*c + b*d;
   return t1 + t2;
   // an optimizing compiler that reassociates will almost
   // surely compute this as (a+b)*(c+d).
}

double power2(long k)
{
   long i;
   double res;

   res = 1;

   for (i = 1; i <= k; i++)
      res = res * 2;

   return res;
}
