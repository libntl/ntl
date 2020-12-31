#include <NTL/GF2EX.h>
#include <NTL/GF2XFactoring.h>

namespace NTL {

void HalfGCD(GF2EX&,GF2EX&);
void PlainRem(GF2EX& r, const GF2EX& a, const GF2EX& b, GF2XVec& x);

}

NTL_CLIENT


#define TIME_IT(t, action) \
do { \
   double _t0, _t1; \
   long _iter = 1; \
   long _cnt = 0; \
   do { \
      _t0 = GetTime(); \
      for (long _i = 0; _i < _iter; _i++) { action; _cnt++; } \
      _t1 = GetTime(); \
   } while ( _t1 - _t0 < 2 && (_iter *= 2)); \
   t = (_t1 - _t0)/_iter; \
} while(0)



void TestGCD1(long bnd, const GF2EX& a, const GF2EX& b)
{

   long n = deg(a) + 1;
   GF2EX u(INIT_SIZE, n), v(INIT_SIZE, n);
   GF2XVec tmp(n, 2*GF2E::WordLength());

   u = a;
   v = b;
   while (deg(v) > bnd) {
      PlainRem(u, u, v, tmp);
      swap(u, v);
   } 

}


long test(long k)
{
   GF2X P;

   BuildIrred(P, k);
   GF2EPush push(P);

   for (long n = 42; ; n = long(n*1.4)) {
      cerr << ",";
      GF2EX d, a, b, u, v;
      random(a, n);
      SetCoeff(a, n);
      random(b, n);
      double t1, t2;
      TIME_IT(t1, u=a; v=b; HalfGCD(u, v));
      TIME_IT(t2, TestGCD1(deg(v), a, b));
      double t = t1/t2;
      if (t <= 1) return n;
   }
}

int main()
{
#if 1
   cerr << "0.5 " << test(32) << "\n";
   for (long i = 1; i <= 4; i+=1) {
      cerr << i << " " << test(64*i) << "\n";
   }
   for (long i = 8; i <= 16; i+=4) {
      cerr << i << " " << test(64*i) << "\n";
   }
#endif
   for (long i = 24; i <= 48; i+=8) {
      cerr << i << " " << test(64*i) << "\n";
   }
}


