#include <NTL/GF2EX.h>
#include <NTL/GF2XFactoring.h>

namespace NTL {

void DivRemPlain(GF2EX& q, GF2EX& r, const GF2EX& a, const GF2EX& b, bool plain);

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


long test(long k)
{
   GF2X P;

   BuildIrred(P, k);
   GF2EPush push(P);

   for (long n = 25; ; n+=25) {
      cerr << ",";
      GF2EX a, b, q, r;
      random(a, 2*n);
      random(b, n);
      double t1, t2;
      TIME_IT(t1, DivRemPlain(q, r, a, b, false));
      TIME_IT(t2, DivRemPlain(q, r, a, b, true));
      double t = t1/t2;
      if (t <= 0.95) return n;
   }
}

int main()
{
   cerr << "0.5 " << test(32) << "\n";
   for (long i = 1; i <= 50; i++) {
      cerr << i << " " << test(64*i) << "\n";
   }

   for (long i = 75; i <= 200 ; i+=25) {
      cerr << i << " " << test(64*i) << "\n";
   }
}


