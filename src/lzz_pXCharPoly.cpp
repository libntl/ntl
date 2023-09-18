#include <NTL/mat_poly_lzz_p.h>


NTL_START_IMPL

static
void HessCharPoly(zz_pX& g, const zz_pX& a, const zz_pX& f)
{
   long n = deg(f);
   if (n <= 0 || deg(a) >= n)
      LogicError("HessCharPoly: bad args");

   mat_zz_p M;
   M.SetDims(n, n);

   long i, j;

   zz_pX t;
   t = a;

   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) 
         M[i][j] = coeff(t, j);

      if (i < n-1) 
         MulByXMod(t, t, f);
   }

   CharPoly(g, M);
}

// "well-known" algorithm according to https://arxiv.org/pdf/1109.4323.pdf p.11
void CharPolyMod(zz_pX& g, const zz_pX& a, const zz_pX& ff)
{
   zz_pX f = ff;
   MakeMonic(f);
   long n = deg(f);

   if (n <= 0 || deg(a) >= n)
      LogicError("CharPolyMod: bad args");

   if (IsZero(a)) {
      clear(g);
      SetCoeff(g, n);
      return;
   }

   if (zz_p::modulus() < n+1) {
      MinPolyMod(g, a, f);
      if (deg(g) < n)
         HessCharPoly(g, a, f);
      return;
   }

   // compute tr(x^0 mod f), tr(x^1 mod f), ..., tr(x^(n-1) mod f) as rev(f')/rev(f)
   // see for instance https://cr.yp.to/papers/newton.pdf
   vec_zz_p T(INIT_SIZE, n);
   {
      zz_pX rev_f = reverse(f);
      zz_pX rev_df = reverse(diff(f));
      zz_pX inv_rev_f = InvTrunc(rev_f, n);
      zz_pX quot = MulTrunc(rev_df, inv_rev_f, n);
      VectorCopy(T, quot, n);
   }

   // compute tr(a^0 mod f), tr(a^1 mod f), ..., tr(a^n mod f)
   vec_zz_p tr = ProjectPowers(T, n+1, a, ff);

   // recover characteristic polynomial using exp(∫g'/g) = g;
   // "well-known" according to https://dl.acm.org/doi/pdf/10.1145/1145768.1145814 Prop. 9
   zz_pX trpoly;
   trpoly.SetLength(n+1);
   trpoly[0] = 0;
   for (long i = 1; i < n+1; ++i)
      trpoly[i] = -tr[i] / i;
   trpoly.normalize();

   ExpTrunc(g, trpoly, n+1);
   reverse(g, g, n);
}

NTL_END_IMPL
