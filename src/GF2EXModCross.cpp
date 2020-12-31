#include <NTL/GF2EX.h>
#include <NTL/GF2XFactoring.h>

namespace NTL {

void BuildPlain(GF2EXModulus& F, const GF2EX& f, bool plain);

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

   for (long n = 5; ; n+=5) {
      cerr << ",";
      GF2EX a, r, f;
      random(a, 2*n-1);
      random(f, n);
      SetCoeff(f, n);
      GF2EXModulus F1, F2;
      BuildPlain(F1, f, false);
      BuildPlain(F2, f, true);
      double t1, t2;
      TIME_IT(t1, rem(r, a, F1));
      TIME_IT(t2, rem(r, a, F2));
      double t = t1/t2;
      if (t <= 0.95) return n;
   }
}

int main()
{
   cerr << "0.5 " << test(32) << "\n";
   for (long i = 1; i <= 40 ; i++) {
      cerr << i << " " << test(64*i) << "\n";
   }
}


