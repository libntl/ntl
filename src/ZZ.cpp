
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/Lazy.h>
#include <NTL/fileio.h>
#include <NTL/SmartPtr.h>

#include <NTL/BasicThreadPool.h>



#if defined(NTL_HAVE_AVX2)
#include <immintrin.h>
#elif defined(NTL_HAVE_SSSE3)
#include <emmintrin.h>
#include <tmmintrin.h>
#endif

#if defined(NTL_HAVE_KMA)
#include <NTL/linux_s390x.h>
#endif



NTL_START_IMPL





const ZZ& ZZ::zero()
{
   
   static const ZZ z; // GLOBAL (relies on C++11 thread-safe init)
   return z;
}


const ZZ& ZZ_expo(long e)
{
   NTL_TLS_LOCAL(ZZ, expo_helper);

   conv(expo_helper, e);
   return expo_helper;
}



void AddMod(ZZ& x, const ZZ& a, long b, const ZZ& n)
{
   NTL_ZZRegister(B);
   conv(B, b);
   AddMod(x, a, B, n);
}


void SubMod(ZZ& x, const ZZ& a, long b, const ZZ& n)
{
   NTL_ZZRegister(B);
   conv(B, b);
   SubMod(x, a, B, n);
}

void SubMod(ZZ& x, long a, const ZZ& b, const ZZ& n)
{
   NTL_ZZRegister(A);
   conv(A, a);
   SubMod(x, A, b, n);
}



// ****** input and output


static NTL_CHEAP_THREAD_LOCAL long iodigits = 0;
static NTL_CHEAP_THREAD_LOCAL long ioradix = 0;
// iodigits is the greatest integer such that 10^{iodigits} < NTL_WSP_BOUND
// ioradix = 10^{iodigits}

static void InitZZIO()
{
   long x;

   x = (NTL_WSP_BOUND-1)/10;
   iodigits = 0;
   ioradix = 1;

   while (x) {
      x = x / 10;
      iodigits++;
      ioradix = ioradix * 10;
   }

   if (iodigits <= 0) TerminalError("problem with I/O");
}


istream& operator>>(istream& s, ZZ& x)
{
   long c;
   long cval;
   long sign;
   long ndigits;
   long acc;
   NTL_ZZRegister(a);

   if (!s) NTL_INPUT_ERROR(s, "bad ZZ input");

   if (!iodigits) InitZZIO();

   a = 0;

   SkipWhiteSpace(s);
   c = s.peek();

   if (c == '-') {
      sign = -1;
      s.get();
      c = s.peek();
   }
   else
      sign = 1;

   cval = CharToIntVal(c);

   if (cval < 0 || cval > 9) NTL_INPUT_ERROR(s, "bad ZZ input");

   ndigits = 0;
   acc = 0;
   while (cval >= 0 && cval <= 9) {
      acc = acc*10 + cval;
      ndigits++;

      if (ndigits == iodigits) {
         mul(a, a, ioradix);
         add(a, a, acc);
         ndigits = 0;
         acc = 0;
      }

      s.get();
      c = s.peek();
      cval = CharToIntVal(c);
   }

   if (ndigits != 0) {
      long mpy = 1;
      while (ndigits > 0) {
         mpy = mpy * 10;
         ndigits--;
      }

      mul(a, a, mpy);
      add(a, a, acc);
   }

   if (sign == -1)
      negate(a, a);

   x = a;
   return s;
}


// The class _ZZ_local_stack should be defined in an empty namespace,
// but since I don't want to rely on namespaces, we just give it a funny 
// name to avoid accidental name clashes.

struct _ZZ_local_stack {
   long top;
   Vec<long> data;

   _ZZ_local_stack() { top = -1; }

   long pop() { return data[top--]; }
   long empty() { return (top == -1); }
   void push(long x);
};

void _ZZ_local_stack::push(long x)
{
   if (top+1 >= data.length()) 
      data.SetLength(max(32, long(1.414*data.length())));

   top++;
   data[top] = x;
}


static
void PrintDigits(ostream& s, long d, long justify)
{
   NTL_TLS_LOCAL_INIT(Vec<char>, buf, (INIT_SIZE, iodigits));

   long i = 0;

   while (d) {
      buf[i] = IntValToChar(d % 10);
      d = d / 10;
      i++;
   }

   if (justify) {
      long j = iodigits - i;
      while (j > 0) {
         s << "0";
         j--;
      }
   }

   while (i > 0) {
      i--;
      s << buf[i];
   }
}
      

   

ostream& operator<<(ostream& s, const ZZ& a)
{
   ZZ b;
   _ZZ_local_stack S;
   long r;
   long k;

   if (!iodigits) InitZZIO();

   b = a;

   k = sign(b);

   if (k == 0) {
      s << "0";
      return s;
   }

   if (k < 0) {
      s << "-";
      negate(b, b);
   }

   do {
      r = DivRem(b, b, ioradix);
      S.push(r);
   } while (!IsZero(b));

   r = S.pop();
   PrintDigits(s, r, 0);

   while (!S.empty()) {
      r = S.pop();
      PrintDigits(s, r, 1);
   }
      
   return s;
}



long GCD(long a, long b)
{
   long u, v, t, x;

   if (a < 0) {
      if (a < -NTL_MAX_LONG) ResourceError("GCD: integer overflow");
      a = -a;
   }

   if (b < 0) {
      if (b < -NTL_MAX_LONG) ResourceError("GCD: integer overflow");
      b = -b;
   }


   if (b==0)
      x = a;
   else {
      u = a;
      v = b;
      do {
         t = u % v;
         u = v; 
         v = t;
      } while (v != 0);

      x = u;
   }

   return x;
}

         

void XGCD(long& d, long& s, long& t, long a, long b)
{
   long  u, v, u0, v0, u1, v1, u2, v2, q, r;

   long aneg = 0, bneg = 0;

   if (a < 0) {
      if (a < -NTL_MAX_LONG) ResourceError("XGCD: integer overflow");
      a = -a;
      aneg = 1;
   }

   if (b < 0) {
      if (b < -NTL_MAX_LONG) ResourceError("XGCD: integer overflow");
      b = -b;
      bneg = 1;
   }

   u1=1; v1=0;
   u2=0; v2=1;
   u = a; v = b;

   while (v != 0) {
      q = u / v;
      r = u % v;
      u = v;
      v = r;
      u0 = u2;
      v0 = v2;
      u2 =  u1 - q*u2;
      v2 = v1- q*v2;
      u1 = u0;
      v1 = v0;
   }

   if (aneg)
      u1 = -u1;

   if (bneg)
      v1 = -v1;

   d = u;
   s = u1;
   t = v1;
}
   
long InvModStatus(long& x, long a, long n)
{
   long d, s, t;

   XGCD(d, s, t, a, n);
   if (d != 1) {
      x = d;
      return 1;
   }
   else {
      if (s < 0)
         x = s + n;
      else
         x = s;

      return 0;
   }
}

long InvMod(long a, long n)
{
   long d, s, t;

   XGCD(d, s, t, a, n);
   if (d != 1) {
      InvModError("InvMod: inverse undefined");
   }
   if (s < 0)
      return s + n;
   else
      return s;
}


long PowerMod(long a, long ee, long n)
{
   long x, y;

   unsigned long e;

   if (ee < 0)
      e = - ((unsigned long) ee);
   else
      e = ee;

   x = 1;
   y = a;
   while (e) {
      if (e & 1) x = MulMod(x, y, n);
      y = MulMod(y, y, n);
      e = e >> 1;
   }

   if (ee < 0) x = InvMod(x, n);

   return x;
}

static
long MillerWitness_sp(long n, long x)
{
   long m, y, z;
   long j, k;

   if (x == 0) return 0;

   m = n - 1;
   k = 0;
   while((m & 1) == 0) {
      m = m >> 1;
      k++;
   }
   // n - 1 == 2^k * m, m odd

   z = PowerMod(x, m, n);
   if (z == 1) return 0;

   j = 0;
   do {
      y = z;
      z = MulMod(y, y, n);
      j++;
   } while (j != k && z != 1);

   if (z != 1 || y != n-1) return 1;
   return 0;
}

long ProbPrime(long n, long NumTrials)
{
   if (NumTrials < 0) NumTrials = 0;

   long m, x, y, z;
   long i, j, k;

   if (n <= 1) return 0;


   if (n == 2) return 1;
   if (n % 2 == 0) return 0;

   if (n == 3) return 1;
   if (n % 3 == 0) return 0;

   if (n == 5) return 1;
   if (n % 5 == 0) return 0;

   if (n == 7) return 1;
   if (n % 7 == 0) return 0;

   if (n == 11) return 1;
   if (n % 11 == 0) return 0;

   if (n == 13) return 1;
   if (n % 13 == 0) return 0;

   if (n >= NTL_SP_BOUND) {
      return ProbPrime(to_ZZ(n), NumTrials);
   }

   m = n - 1;
   k = 0;
   while((m & 1) == 0) {
      m = m >> 1;
      k++;
   }

   // n - 1 == 2^k * m, m odd

   for (i = 0; i < NumTrials+1; i++) {
      // for consistency with the multi-precision version,
      // we first see if 2 is a witness, so we really do 
      // NumTrials+1 tests

      if (i == 0) 
         x = 2;
      else {
	 do {
	    x = RandomBnd(n);
	 } while (x == 0);
         // x == 0 is not a useful candidate for a witness!
      }

      z = PowerMod(x, m, n);
      if (z == 1) continue;
   
      j = 0;
      do {
         y = z;
         z = MulMod(y, y, n);
         j++;
      } while (j != k && z != 1);

      if (z != 1 || y !=  n-1) return 0;
   }

   return 1;
}


long MillerWitness(const ZZ& n, const ZZ& x)
{
   if (n.SinglePrecision()) {
      return MillerWitness_sp(to_long(n), to_long(x));
   }

   ZZ m, y, z;

   long j, k;

   if (x == 0) return 0;

   add(m, n, -1);
   k = MakeOdd(m);
   // n - 1 == 2^k * m, m odd

   PowerMod(z, x, m, n);
   if (z == 1) return 0;

   j = 0;
   do {
      y = z;
      SqrMod(z, y, n);
      j++;
   } while (j != k && z != 1);

   if (z != 1) return 1;
   add(y, y, 1);
   if (y != n) return 1;
   return 0;
}


// ComputePrimeBound computes a reasonable bound for trial
// division in the Miller-Rabin test.
// It is computed a bit on the "low" side, since being a bit
// low doesn't hurt much, but being too high can hurt a lot.

// See the paper "Fast generation of prime numbers and secure
// public-key cryptographic parameters" by Ueli Maurer.
// In that paper, it is calculated that the optimal bound in
// roughly T_exp/T_div, where T_exp is the time for an exponentiation
// and T_div is is the time for a single precision division.
// Of course, estimating these times is a bit tricky, and
// the values we use are based on experimentation, assuming
// GMP is being used.  I've tested this on various bit lengths
// up to 16,000, and they seem to be pretty close to optimal. 

static
long ComputePrimeBound(long bn)
{
   long wn = (bn+NTL_ZZ_NBITS-1)/NTL_ZZ_NBITS;

   long fn;

   if (wn <= 36)
      fn = wn/4 + 1;
   else
      fn = long(1.67*sqrt(double(wn)));

   long prime_bnd;

   if (NumBits(bn) + NumBits(fn) > NTL_SP_NBITS)
      prime_bnd = NTL_SP_BOUND;
   else
      prime_bnd = bn*fn;

   return prime_bnd;


}


long ProbPrime(const ZZ& n, long NumTrials)
{
   if (NumTrials < 0) NumTrials = 0;

   if (n <= 1) return 0;

   if (n.SinglePrecision()) {
      return ProbPrime(to_long(n), NumTrials);
   }


   long prime_bnd = ComputePrimeBound(NumBits(n));


   PrimeSeq s;
   long p;

   p = s.next();
   while (p && p < prime_bnd) {
      if (rem(n, p) == 0)
         return 0;

      p = s.next();
   }

   ZZ W;
   W = 2;

   // first try W == 2....the exponentiation
   // algorithm runs slightly faster in this case

   if (MillerWitness(n, W))
      return 0;


   long i;

   for (i = 0; i < NumTrials; i++) {
      do {
         RandomBnd(W, n);
      } while (W == 0);
      // W == 0 is not a useful candidate for a witness!

      if (MillerWitness(n, W)) 
         return 0;
   }

   return 1;
}


static
void MultiThreadedRandomPrime(ZZ& n, long l, long NumTrials)
{

   long nt = AvailableThreads();



   const long LOCAL_ITER_BOUND = 8; 
   // since resetting the PRG comes at a certain cost,
   // we perform a few iterations with each reset to
   // amortize the reset cost. 

   unsigned long initial_counter = 0;
   ZZ seed;
   RandomBits(seed, 256);

   for (;;) {

      AtomicLowWaterMark low_water_mark(-1UL);
      AtomicCounter counter(initial_counter);

      Vec< UniquePtr<ZZ> > result(INIT_SIZE, nt);
      Vec<unsigned long> result_ctr(INIT_SIZE, nt, -1UL);

      NTL_EXEC_INDEX(nt, index)

         RandomStreamPush push;

	 SetSeed(seed);
	 RandomStream& stream = GetCurrentRandomStream();

	 ZZ cand;

	 while (low_water_mark == -1UL) {

	    unsigned long local_ctr = counter.inc();
            if (local_ctr >> (NTL_BITS_PER_NONCE-1)) {
               // counter overflow...rather academic
               break;
            }

	    stream.set_nonce(local_ctr);
	    
	    for (long iter = 0; iter < LOCAL_ITER_BOUND && 
				local_ctr <= low_water_mark; iter++) {

	       RandomLen(cand, l);
	       if (!IsOdd(cand)) add(cand, cand, 1);

	       if (ProbPrime(cand, 0)) { 
		  result[index].make(cand);
		  result_ctr[index] = local_ctr;
		  low_water_mark.UpdateMin(local_ctr);
		  break;
	       }
	    }
	 }

      NTL_EXEC_INDEX_END

      // find index of low_water_mark

      unsigned long low_water_mark1 = low_water_mark;
      long low_water_index = -1;

      for (long index = 0; index < nt; index++) {
	 if (result_ctr[index] == low_water_mark1) {
	    low_water_index = index;
	    break;
	 }
      }

      if (low_water_index == -1) {
         // counter overflow...rather academic
         initial_counter = 0;
         RandomBits(seed, 256);
         continue;
      }

      ZZ N;
      N = *result[low_water_index];

      Vec<ZZ> W(INIT_SIZE, NumTrials);

      for (long i = 0; i < NumTrials; i++) {
         do { 
            RandomBnd(W[i], N);
         } while (W[i] == 0);
      }

      AtomicBool tests_pass(true);

      NTL_EXEC_RANGE(NumTrials, first, last)

         for (long i = first; i < last && tests_pass; i++) {
            if (MillerWitness(N, W[i])) tests_pass = false;
         }

      NTL_EXEC_RANGE_END

      if (tests_pass) {
         n = N;
         return;
      }

      // very unlikey to get here
      initial_counter = low_water_mark1 + 1;
   }
}


void RandomPrime(ZZ& n, long l, long NumTrials)
{
   if (NumTrials < 0) NumTrials = 0;

   if (l >= 256) { 
      MultiThreadedRandomPrime(n, l, NumTrials); 
      return;
   }

   if (l <= 1)
      LogicError("RandomPrime: l out of range");

   if (l == 2) {
      if (RandomBnd(2))
         n = 3;
      else
         n = 2;

      return;
   }

   do {
      RandomLen(n, l);
      if (!IsOdd(n)) add(n, n, 1);
   } while (!ProbPrime(n, NumTrials));
}

void OldRandomPrime(ZZ& n, long l, long NumTrials)
{
   if (l <= 1)
      LogicError("RandomPrime: l out of range");

   if (l == 2) {
      if (RandomBnd(2))
         n = 3;
      else
         n = 2;

      return;
   }

   do {
      RandomLen(n, l);
      if (!IsOdd(n)) add(n, n, 1);
   } while (!ProbPrime(n, NumTrials));
}

void NextPrime(ZZ& n, const ZZ& m, long NumTrials)
{
   ZZ x;

   if (m <= 2) {
      n = 2;
      return;
   }

   x = m;

   while (!ProbPrime(x, NumTrials))
      add(x, x, 1);

   n = x;
}

long NextPrime(long m, long NumTrials)
{
   long x;

   if (m <= 2) 
      return 2;

   x = m;

   while (x < NTL_SP_BOUND && !ProbPrime(x, NumTrials))
      x++;

   if (x >= NTL_SP_BOUND)
      ResourceError("NextPrime: no more primes");

   return x;
}



long NextPowerOfTwo(long m)
{
   long k; 
   unsigned long n, um;

   if (m < 0) return 0;

   um = m;
   n = 1;
   k = 0;

   while (n < um) {
      n = n << 1;
      k++;
   }

   if (k >= NTL_BITS_PER_LONG-1)
      ResourceError("NextPowerOfTwo: overflow");

   return k;
}




long bit(long a, long k)
{
   unsigned long aa;
   if (a < 0)
      aa = - ((unsigned long) a);
   else
      aa = a;

   if (k < 0 || k >= NTL_BITS_PER_LONG) 
      return 0;
   else
      return long((aa >> k) & 1);
}



long divide(ZZ& q, const ZZ& a, const ZZ& b)
{
   NTL_ZZRegister(qq);
   NTL_ZZRegister(r);

   if (IsZero(b)) {
      if (IsZero(a)) {
         clear(q);
         return 1;
      }
      else
         return 0;
   }


   if (IsOne(b)) {
      q = a;
      return 1;
   }

   DivRem(qq, r, a, b);
   if (!IsZero(r)) return 0;
   q = qq;
   return 1;
}

long divide(const ZZ& a, const ZZ& b)
{
   NTL_ZZRegister(r);

   if (IsZero(b)) return IsZero(a);
   if (IsOne(b)) return 1;

   rem(r, a, b);
   return IsZero(r);
}

long divide(ZZ& q, const ZZ& a, long b)
{
   NTL_ZZRegister(qq);

   if (!b) {
      if (IsZero(a)) {
         clear(q);
         return 1;
      }
      else
         return 0;
   }

   if (b == 1) {
      q = a;
      return 1;
   }

   long r = DivRem(qq, a, b);
   if (r) return 0;
   q = qq;
   return 1;
}

long divide(const ZZ& a, long b)
{
   if (!b) return IsZero(a);
   if (b == 1) {
      return 1;
   }

   long r = rem(a,  b);
   return (r == 0);
}


void InvMod(ZZ& x, const ZZ& a, const ZZ& n)
{
   // NOTE: the underlying LIP routines write to the first argument,
   // even if inverse is undefined

   NTL_ZZRegister(xx);
   if (InvModStatus(xx, a, n)) 
      InvModError("InvMod: inverse undefined", a, n);
   x = xx;
}

void PowerMod(ZZ& x, const ZZ& a, const ZZ& e, const ZZ& n)
{
   // NOTE: this ensures that all modular inverses are computed
   // in the routine InvMod above, rather than the LIP-internal
   // modular inverse routine
   if (e < 0) {
      ZZ a_inv;
      ZZ e_neg;

      InvMod(a_inv, a, n);
      negate(e_neg, e);
      LowLevelPowerMod(x, a_inv, e_neg, n);
   }
   else
      LowLevelPowerMod(x, a, e, n); 
}
   
#ifdef NTL_EXCEPTIONS

void InvModError(const char *s, const ZZ& a, const ZZ& n)
{
   throw InvModErrorObject(s, a, n); 
}

#else

void InvModError(const char *s, const ZZ& a, const ZZ& n)
{
   TerminalError(s);
}


#endif

long RandomPrime_long(long l, long NumTrials)
{
   if (NumTrials < 0) NumTrials = 0;

   if (l <= 1 || l >= NTL_BITS_PER_LONG)
      ResourceError("RandomPrime: length out of range");

   long n;
   do {
      n = RandomLen_long(l);
   } while (!ProbPrime(n, NumTrials));

   return n;
}


static Lazy< Vec<char> > lowsieve_storage;
// This is a GLOBAL VARIABLE


PrimeSeq::PrimeSeq()
{
   movesieve = 0;
   pshift = -1;
   pindex = -1;
   exhausted = 0;
}


long PrimeSeq::next()
{
   if (exhausted) {
      return 0;
   }

   if (pshift < 0) {
      shift(0);
      return 2;
   }

   for (;;) {
      const char *p = movesieve;
      long i = pindex;

      while ((++i) < NTL_PRIME_BND) {
         if (p[i]) {
            pindex = i;
            return pshift + 2 * i + 3;
         }
      }

      long newshift = pshift + 2*NTL_PRIME_BND;

      if (newshift > 2 * NTL_PRIME_BND * (2 * NTL_PRIME_BND + 1)) {
         /* end of the road */
         exhausted = 1;
         return 0;
      }

      shift(newshift);
   }
}

void PrimeSeq::shift(long newshift)
{
   long i;
   long j;
   long jstep;
   long jstart;
   long ibound;
   char *p;

   if (!lowsieve_storage.built())
      start();

   const char *lowsieve = lowsieve_storage->elts();


   if (newshift < 0) {
      pshift = -1;
   }
   else if (newshift == 0) {
      pshift = 0;
      movesieve = lowsieve;
   } 
   else if (newshift != pshift) {
      if (movesieve_mem.length() == 0) {
         movesieve_mem.SetLength(NTL_PRIME_BND);
      }

      pshift = newshift;
      movesieve = p = movesieve_mem.elts();
      for (i = 0; i < NTL_PRIME_BND; i++)
         p[i] = 1;

      jstep = 3;
      ibound = pshift + 2 * NTL_PRIME_BND + 1;
      for (i = 0; jstep * jstep <= ibound; i++) {
         if (lowsieve[i]) {
            if (!((jstart = (pshift + 2) / jstep + 1) & 1))
               jstart++;
            if (jstart <= jstep)
               jstart = jstep;
            jstart = (jstart * jstep - pshift - 3) / 2;
            for (j = jstart; j < NTL_PRIME_BND; j += jstep)
               p[j] = 0;
         }
         jstep += 2;
      }
   }

   pindex = -1;
   exhausted = 0;
}


void PrimeSeq::start()
{
   long i;
   long j;
   long jstep;
   long jstart;
   long ibnd;
   char *p;

   do {
      Lazy< Vec<char> >::Builder builder(lowsieve_storage);
      if (!builder()) break;

      UniquePtr< Vec<char> > ptr;
      ptr.make();
      ptr->SetLength(NTL_PRIME_BND);

      p = ptr->elts();

      for (i = 0; i < NTL_PRIME_BND; i++)
         p[i] = 1;
         
      jstep = 1;
      jstart = -1;
      ibnd = (SqrRoot(2 * NTL_PRIME_BND + 1) - 3) / 2;
      for (i = 0; i <= ibnd; i++) {
         jstart += 2 * ((jstep += 2) - 1);
         if (p[i])
            for (j = jstart; j < NTL_PRIME_BND; j += jstep)
               p[j] = 0;
      }

      builder.move(ptr);
   } while (0);

}

void PrimeSeq::reset(long b)
{
   if (b > (2*NTL_PRIME_BND+1)*(2*NTL_PRIME_BND+1)) {
      exhausted = 1;
      return;
   }

   if (b <= 2) {
      shift(-1);
      return;
   }

   if ((b & 1) == 0) b++;

   shift(((b-3) / (2*NTL_PRIME_BND))* (2*NTL_PRIME_BND));
   pindex = (b - pshift - 3)/2 - 1;
}
 
long Jacobi(const ZZ& aa, const ZZ& nn)
{
   ZZ a, n;
   long t, k;
   long d;

   a = aa;
   n = nn;
   t = 1;

   while (a != 0) {
      k = MakeOdd(a);
      d = trunc_long(n, 3);
      if ((k & 1) && (d == 3 || d == 5)) t = -t;

      if (trunc_long(a, 2) == 3 && (d & 3) == 3) t = -t;
      swap(a, n);
      rem(a, a, n);
   }

   if (n == 1)
      return t;
   else
      return 0;
}


void SqrRootMod(ZZ& x, const ZZ& aa, const ZZ& nn)
{
   if (aa == 0 || aa == 1) {
      x = aa;
      return;
   }

   // at this point, we must have nn >= 5

   if (trunc_long(nn, 2) == 3) {  // special case, n = 3 (mod 4)
      ZZ n, a, e, z;

      n = nn;
      a  = aa;

      add(e, n, 1);
      RightShift(e, e, 2);

      PowerMod(z, a, e, n);
      x = z;

      return;
   }

   ZZ n, m;
   int h, nlen;

   n = nn;
   nlen = NumBits(n);

   sub(m, n, 1);
   h = MakeOdd(m);  // h >= 2


   if (nlen > 50 && h < SqrRoot(nlen)) {
      long i, j;
      ZZ a, b, a_inv, c, r, m1, d;

      a = aa;
      InvMod(a_inv, a, n);

      if (h == 2) 
         b = 2;
      else {
         do {
            RandomBnd(b, n);
         } while (Jacobi(b, n) != -1);
      }


      PowerMod(c, b, m, n);
      
      add(m1, m, 1);
      RightShift(m1, m1, 1);
      PowerMod(r, a, m1, n);

      for (i = h-2; i >= 0; i--) {
         SqrMod(d, r, n);
         MulMod(d, d, a_inv, n);
         for (j = 0; j < i; j++)
            SqrMod(d, d, n);
         if (!IsOne(d))
            MulMod(r, r, c, n);
         SqrMod(c, c, n);
      } 

      x = r;
      return;
   } 





   long i, k;
   ZZ ma, t, u, v, e;
   ZZ t1, t2, t3, t4;

   n = nn;
   NegateMod(ma, aa, n);

   // find t such that t^2 - 4*a is not a square

   MulMod(t1, ma, 4, n);
   do {
      RandomBnd(t, n);
      SqrMod(t2, t, n);
      AddMod(t2, t2, t1, n);
   } while (Jacobi(t2, n) != -1);

   // compute u*X + v = X^{(n+1)/2} mod f, where f = X^2 - t*X + a

   add(e, n, 1);
   RightShift(e, e, 1);

   u = 0;
   v = 1;

   k = NumBits(e);

   for (i = k - 1; i >= 0; i--) {
      add(t2, u, v);
      sqr(t3, t2);  // t3 = (u+v)^2
      sqr(t1, u);
      sqr(t2, v);
      sub(t3, t3, t1);
      sub(t3, t3, t2); // t1 = u^2, t2 = v^2, t3 = 2*u*v
      rem(t1, t1, n);
      mul(t4, t1, t);
      add(t4, t4, t3);
      rem(u, t4, n);

      mul(t4, t1, ma);
      add(t4, t4, t2);
      rem(v, t4, n);
      
      if (bit(e, i)) {
         MulMod(t1, u, t, n);
         AddMod(t1, t1, v, n);
         MulMod(v, u, ma, n);
         u = t1;
      }

   }

   x = v;
}



// Chinese Remaindering.
//
// This version in new to v3.7, and is significantly
// simpler and faster than the previous version.
//
// This function takes as input g, a, G, p,
// such that a > 0, 0 <= G < p, and gcd(a, p) = 1.
// It computes a' = a*p and g' such that 
//   * g' = g (mod a);
//   * g' = G (mod p);
//   * -a'/2 < g' <= a'/2.
// It then sets g := g' and a := a', and returns 1 iff g has changed.
//
// Under normal use, the input value g satisfies -a/2 < g <= a/2;
// however, this was not documented or enforced in earlier versions,
// so to maintain backward compatability, no restrictions are placed
// on g.  This routine runs faster, though, if -a/2 < g <= a/2,
// and the first thing the routine does is to make this condition
// hold.
//
// Also, under normal use, both a and p are odd;  however, the routine
// will still work even if this is not so.
//
// The routine is based on the following simple fact.
//
// Let -a/2 < g <= a/2, and let h satisfy
//   * g + a h = G (mod p);
//   * -p/2 < h <= p/2.
// Further, if p = 2*h and g > 0, set
//   g' := g - a h;
// otherwise, set
//   g' := g + a h.
// Then g' so defined satisfies the above requirements.
//
// It is trivial to see that g's satisfies the congruence conditions.
// The only thing is to check that the "balancing" condition
// -a'/2 < g' <= a'/2 also holds.


long CRT(ZZ& gg, ZZ& a, long G, long p)
{
   if (p >= NTL_SP_BOUND) {
      ZZ GG, pp;
      conv(GG, G);
      conv(pp, p);
      return CRT(gg, a, GG, pp);
   }

   long modified = 0;

   NTL_ZZRegister(g);

   if (!CRTInRange(gg, a)) {
      modified = 1;
      ZZ a1;
      rem(g, gg, a);
      RightShift(a1, a, 1);
      if (g > a1) sub(g, g, a);
   }
   else
      g = gg;


   long p1;
   p1 = p >> 1;

   long a_inv;
   a_inv = rem(a, p);
   a_inv = InvMod(a_inv, p);

   long h;
   h = rem(g, p);
   h = SubMod(G, h, p);
   h = MulMod(h, a_inv, p);
   if (h > p1)
      h = h - p;

   if (h != 0) {
      modified = 1;

      if (!(p & 1) && g > 0 && (h == p1))
         MulSubFrom(g, a, h);
      else
         MulAddTo(g, a, h);
   }

   mul(a, a, p);
   gg = g;

   return modified;
}

long CRT(ZZ& gg, ZZ& a, const ZZ& G, const ZZ& p)
{
   long modified = 0;

   ZZ g;

   if (!CRTInRange(gg, a)) {
      modified = 1;
      ZZ a1;
      rem(g, gg, a);
      RightShift(a1, a, 1);
      if (g > a1) sub(g, g, a);
   }
   else
      g = gg;


   ZZ p1;
   RightShift(p1, p, 1);

   ZZ a_inv;
   rem(a_inv, a, p);
   InvMod(a_inv, a_inv, p);

   ZZ h;
   rem(h, g, p);
   SubMod(h, G, h, p);
   MulMod(h, h, a_inv, p);
   if (h > p1)
      sub(h, h, p);

   if (h != 0) {
      modified = 1;
      ZZ ah;
      mul(ah, a, h);

      if (!IsOdd(p) && g > 0 &&  (h == p1))
         sub(g, g, ah);
      else
         add(g, g, ah);
   }

   mul(a, a, p);
   gg = g;

   return modified;
}



void sub(ZZ& x, long a, const ZZ& b)
{
   NTL_ZZRegister(A);
   conv(A, a);
   sub(x, A, b);
}


void power2(ZZ& x, long e)
{
   if (e < 0) ArithmeticError("power2: negative exponent");
   set(x);
   LeftShift(x, x, e);
}

   

void bit_and(ZZ& x, const ZZ& a, long b)
{
   NTL_ZZRegister(B);
   conv(B, b);
   bit_and(x, a, B);
}

void bit_or(ZZ& x, const ZZ& a, long b)
{
   NTL_ZZRegister(B);
   conv(B, b);
   bit_or(x, a, B);
}

void bit_xor(ZZ& x, const ZZ& a, long b)
{
   NTL_ZZRegister(B);
   conv(B, b);
   bit_xor(x, a, B);
}


long power_long(long a, long e)
{
   if (e < 0) ArithmeticError("power_long: negative exponent");

   if (e == 0) return 1;

   if (a == 1) return 1;
   if (a == -1) {
      if (e & 1)
         return -1;
      else
         return 1;
   }

   // no overflow check --- result is computed correctly
   // modulo word size

   unsigned long res = 1;
   unsigned long aa = a;
   long i;

   for (i = 0; i < e; i++)
      res *= aa;

   return to_long(res);
}



// ======================= new PRG stuff ======================




#if (NTL_BITS_PER_INT32 == 32)
#define INT32MASK(x) (x)
#else
#define INT32MASK(x) ((x) & _ntl_uint32(0xffffffff))
#endif



// SHA256 code adapted from an implementauin by Brad Conte.
// The following is from his original source files.
/*********************************************************************
* Filename:   sha256.c
* Author:     Brad Conte (brad AT bradconte.com)
* Copyright:
* Disclaimer: This code is presented "as is" without any guarantees.
* Details:    Implementation of the SHA-256 hashing algorithm.
              SHA-256 is one of the three algorithms in the SHA2
              specification. The others, SHA-384 and SHA-512, are not
              offered in this implementation.
              Algorithm specification can be found here:
               * http://csrc.nist.gov/publications/fips/fips180-2/fips180-2withchangenotice.pdf
              This implementation uses little endian byte order.
*********************************************************************/

// And the following is from the description at 
// https://github.com/B-Con/crypto-algorithms

/*********************************************************************

These are basic implementations of standard cryptography algorithms, written by
Brad Conte (brad@bradconte.com) from scratch and without any cross-licensing.
They exist to provide publically accessible, restriction-free implementations
of popular cryptographic algorithms, like AES and SHA-1. These are primarily
intended for educational and pragmatic purposes (such as comparing a
specification to actual implementation code, or for building an internal
application that computes test vectors for a product). The algorithms have been
tested against standard test vectors.

This code is released into the public domain free of any restrictions. The
author requests acknowledgement if the code is used, but does not require it.
This code is provided free of any liability and without any quality claims by
the author.

Note that these are not cryptographically secure implementations. They have no
resistence to side-channel attacks and should not be used in contexts that need
cryptographically secure implementations.

These algorithms are not optimized for speed or space. They are primarily
designed to be easy to read, although some basic optimization techniques have
been employed.

*********************************************************************/






#define SHA256_BLOCKSIZE (64)
#define SHA256_HASHSIZE  (32)

// DBL_INT_ADD treats two unsigned ints a and b as one 64-bit integer and adds c to it
static inline
void DBL_INT_ADD(_ntl_uint32& a, _ntl_uint32& b, _ntl_uint32 c)
{
   _ntl_uint32 aa = INT32MASK(a);
   if (aa > INT32MASK(_ntl_uint32(0xffffffff) - c)) b++;
   a = aa + c;
}

#define ROTLEFT(a,b) (((a) << (b)) | (INT32MASK(a) >> (32-(b))))
#define ROTRIGHT(a,b) ((INT32MASK(a) >> (b)) | ((a) << (32-(b))))

#define CH(x,y,z) (((x) & (y)) ^ (~(x) & (z)))
#define MAJ(x,y,z) (((x) & (y)) ^ ((x) & (z)) ^ ((y) & (z)))
#define EP0(x) (ROTRIGHT(x,2) ^ ROTRIGHT(x,13) ^ ROTRIGHT(x,22))
#define EP1(x) (ROTRIGHT(x,6) ^ ROTRIGHT(x,11) ^ ROTRIGHT(x,25))
#define SIG0(x) (ROTRIGHT(x,7) ^ ROTRIGHT(x,18) ^ (INT32MASK(x) >> 3))
#define SIG1(x) (ROTRIGHT(x,17) ^ ROTRIGHT(x,19) ^ (INT32MASK(x) >> 10))

struct SHA256_CTX {
   unsigned char data[64];
   _ntl_uint32 datalen;
   _ntl_uint32 bitlen[2];
   _ntl_uint32 state[8];
};

static const _ntl_uint32 sha256_const[64] = {
   0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
   0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
   0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
   0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
   0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
   0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
   0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
   0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
};


static
void sha256_transform(SHA256_CTX& ctx, unsigned char *data)
{  
   _ntl_uint32 a,b,c,d,e,f,g,h,i,j,t1,t2,m[64];
      
   for (i=0,j=0; i < 16; ++i, j += 4)
      m[i] = (data[j] << 24) | (data[j+1] << 16) | (data[j+2] << 8) | (data[j+3]);
   for ( ; i < 64; ++i)
      m[i] = SIG1(m[i-2]) + m[i-7] + SIG0(m[i-15]) + m[i-16];

   a = ctx.state[0];
   b = ctx.state[1];
   c = ctx.state[2];
   d = ctx.state[3];
   e = ctx.state[4];
   f = ctx.state[5];
   g = ctx.state[6];
   h = ctx.state[7];
   
   for (i = 0; i < 64; ++i) {
      t1 = h + EP1(e) + CH(e,f,g) + sha256_const[i] + m[i];
      t2 = EP0(a) + MAJ(a,b,c);
      h = g;
      g = f;
      f = e;
      e = d + t1;
      d = c;
      c = b;
      b = a;
      a = t1 + t2;
   }
   
   ctx.state[0] += a;
   ctx.state[1] += b;
   ctx.state[2] += c;
   ctx.state[3] += d;
   ctx.state[4] += e;
   ctx.state[5] += f;
   ctx.state[6] += g;
   ctx.state[7] += h;
}  

static
void sha256_init(SHA256_CTX& ctx)
{  
   ctx.datalen = 0; 
   ctx.bitlen[0] = 0; 
   ctx.bitlen[1] = 0; 
   ctx.state[0] = 0x6a09e667;
   ctx.state[1] = 0xbb67ae85;
   ctx.state[2] = 0x3c6ef372;
   ctx.state[3] = 0xa54ff53a;
   ctx.state[4] = 0x510e527f;
   ctx.state[5] = 0x9b05688c;
   ctx.state[6] = 0x1f83d9ab;
   ctx.state[7] = 0x5be0cd19;
}

static
void sha256_update(SHA256_CTX& ctx, const unsigned char *data, _ntl_uint32 len)
{  
   _ntl_uint32 i;
   
   for (i=0; i < len; ++i) { 
      ctx.data[ctx.datalen] = data[i]; 
      ctx.datalen++; 
      if (ctx.datalen == 64) { 
         sha256_transform(ctx,ctx.data);
         DBL_INT_ADD(ctx.bitlen[0],ctx.bitlen[1],512); 
         ctx.datalen = 0; 
      }  
   }  
}  

static
void sha256_final(SHA256_CTX& ctx, unsigned char *hash, 
                  long hlen=SHA256_HASHSIZE)
{  
   _ntl_uint32 i, j; 
   
   i = ctx.datalen; 
   
   // Pad whatever data is left in the buffer. 
   if (ctx.datalen < 56) { 
      ctx.data[i++] = 0x80; 
      while (i < 56) 
         ctx.data[i++] = 0x00; 
   }  
   else { 
      ctx.data[i++] = 0x80; 
      while (i < 64) 
         ctx.data[i++] = 0x00; 
      sha256_transform(ctx,ctx.data);
      memset(ctx.data,0,56); 
   }  
   
   // Append to the padding the total message's length in bits and transform. 
   DBL_INT_ADD(ctx.bitlen[0],ctx.bitlen[1],ctx.datalen * 8);

   ctx.data[63] = ctx.bitlen[0]; 
   ctx.data[62] = ctx.bitlen[0] >> 8; 
   ctx.data[61] = ctx.bitlen[0] >> 16; 
   ctx.data[60] = ctx.bitlen[0] >> 24; 
   ctx.data[59] = ctx.bitlen[1]; 
   ctx.data[58] = ctx.bitlen[1] >> 8; 
   ctx.data[57] = ctx.bitlen[1] >> 16;  
   ctx.data[56] = ctx.bitlen[1] >> 24; 
   sha256_transform(ctx,ctx.data);
   
   for (i = 0; i < 8; i++) {
      _ntl_uint32 w = ctx.state[i];
      for (j = 0; j < 4; j++) {
         if (hlen <= 0) break;
         hash[4*i + j] = w >> (24-j*8); 
         hlen--;
      }
   }

}  



static
void sha256(const unsigned char *data, long dlen, unsigned char *hash, 
            long hlen=SHA256_HASHSIZE)
{
   if (dlen < 0) dlen = 0;
   if (hlen < 0) hlen = 0;

   SHA256_CTX ctx;
   sha256_init(ctx);

   const long BLKSIZE = 4096;

   long i;
   for (i = 0; i <= dlen-BLKSIZE; i += BLKSIZE) 
      sha256_update(ctx, data + i, BLKSIZE);

   if (i < dlen)
      sha256_update(ctx, data + i, dlen - i);

   sha256_final(ctx, hash, hlen);
}


static
void hmac_sha256(const unsigned char *key, long klen, 
                 const unsigned char *data, long dlen,
                 unsigned char *hash, long hlen=SHA256_HASHSIZE)
{
   if (klen < 0) klen = 0;
   if (dlen < 0) dlen = 0;
   if (hlen < 0) hlen = 0;

   unsigned char K[SHA256_BLOCKSIZE];
   unsigned char tmp[SHA256_HASHSIZE];

   long i;

   if (klen <= SHA256_BLOCKSIZE) {
      for (i = 0; i < klen; i++)
         K[i] = key[i];
      for (i = klen; i < SHA256_BLOCKSIZE; i++) 
         K[i] = 0;
   }
   else {
      sha256(key, klen, K, SHA256_BLOCKSIZE); 
      for (i = SHA256_HASHSIZE; i < SHA256_BLOCKSIZE; i++)
         K[i] = 0;
   }

   for (i = 0; i < SHA256_BLOCKSIZE; i++)
      K[i] ^= 0x36;

   SHA256_CTX ctx;
   sha256_init(ctx);
   sha256_update(ctx, K, SHA256_BLOCKSIZE);
   sha256_update(ctx, data, dlen);
   sha256_final(ctx, tmp);

   for (i = 0; i < SHA256_BLOCKSIZE; i++)
      K[i] ^= (0x36 ^ 0x5C);

   sha256_init(ctx);
   sha256_update(ctx, K, SHA256_BLOCKSIZE);
   sha256_update(ctx, tmp, SHA256_HASHSIZE);
   sha256_final(ctx, hash, hlen);
}


// This key derivation uses HMAC with a zero key to derive
// an intermediate key K from the data, and then uses HMAC
// as a PRF in counter mode with key K to derive the final key

void DeriveKey(unsigned char *key, long klen,  
               const unsigned char *data, long dlen)
{
   if (dlen < 0) LogicError("DeriveKey: bad args");
   if (klen < 0) LogicError("DeriveKey: bad args");

   long i, j;


   unsigned char K[SHA256_HASHSIZE];
   hmac_sha256(0, 0, data, dlen, K); 

   // initialize 64-bit counter to zero
   unsigned char counter[8];
   for (j = 0; j < 8; j++) counter[j] = 0;

   for (i = 0; i <= klen-SHA256_HASHSIZE; i += SHA256_HASHSIZE) {
      hmac_sha256(K, SHA256_HASHSIZE, counter, 8, key+i); 

      // increment counter
      for (j = 0; j < 8; j++) {
         counter[j]++;
         if (counter[j] != 0) break; 
      }
   }

   if (i < klen) 
      hmac_sha256(K, SHA256_HASHSIZE, counter, 8, key+i, klen-i);
}




// ******************** ChaCha20 stuff ***********************

// ============= old stuff

#define LE(p) (((_ntl_uint32)((p)[0])) + ((_ntl_uint32)((p)[1]) << 8) + \
    ((_ntl_uint32)((p)[2]) << 16) + ((_ntl_uint32)((p)[3]) << 24))

#define FROMLE(p, x) (p)[0] = (x), (p)[1] = ((x) >> 8), \
   (p)[2] = ((x) >> 16), (p)[3] = ((x) >> 24)


#define QUARTERROUND(x, a, b, c, d) \
    x[a] += x[b], x[d] = ROTLEFT(x[d] ^ x[a], 16), \
    x[c] += x[d], x[b] = ROTLEFT(x[b] ^ x[c], 12), \
    x[a] += x[b], x[d] = ROTLEFT(x[d] ^ x[a], 8), \
    x[c] += x[d], x[b] = ROTLEFT(x[b] ^ x[c], 7)


static
void salsa20_core(_ntl_uint32* data)
{
   long i;

   for (i = 0; i < 10; i++) {
      QUARTERROUND(data, 0, 4, 8, 12);
      QUARTERROUND(data, 1, 5, 9, 13);
      QUARTERROUND(data, 2, 6, 10, 14);
      QUARTERROUND(data, 3, 7, 11, 15);
      QUARTERROUND(data, 0, 5, 10, 15);
      QUARTERROUND(data, 1, 6, 11, 12);
      QUARTERROUND(data, 2, 7, 8, 13);
      QUARTERROUND(data, 3, 4, 9, 14);
   }
}


// key K must be exactly 32 bytes
static
void salsa20_init(_ntl_uint32 *state, const unsigned char *K)  
{
   static const _ntl_uint32 chacha_const[4] = 
      { 0x61707865, 0x3320646e, 0x79622d32, 0x6b206574 };

   long i;

   for (i = 0; i < 4; i++)
      state[i] = chacha_const[i];

   for (i = 4; i < 12; i++)
      state[i] = LE(K + 4*(i-4));

   for (i = 12; i < 16; i++)
      state[i] = 0;
}



// state and data are of length 16
static
void salsa20_apply(_ntl_uint32 *state, _ntl_uint32 *data)
{
   long i;

   for (i = 0; i < 16; i++) data[i] = state[i];

   salsa20_core(data);

   for (i = 0; i < 16; i++) data[i] += state[i];

   for (i = 12; i < 14; i++) {
      state[i]++;
      state[i] = INT32MASK(state[i]);
      if (state[i] != 0) break;
   }
}



old_RandomStream::old_RandomStream(const unsigned char *key)
{
   salsa20_init(state, key);
   pos = 64;
}


void old_RandomStream::do_get(unsigned char *res, long n)
{
   if (n < 0) LogicError("RandomStream::get: bad args");

   long i, j;

   if (n <= 64-pos) {
      for (i = 0; i < n; i++) res[i] = buf[pos+i];
      pos += n;
      return;
   }

   // read remainder of buffer
   for (i = 0; i < 64-pos; i++) res[i] = buf[pos+i];
   n -= 64-pos;
   res += 64-pos;
   pos = 64;

   _ntl_uint32 wdata[16];

   // read 64-byte chunks
   for (i = 0; i <= n-64; i += 64) {
      salsa20_apply(state, wdata);
      for (j = 0; j < 16; j++)
         FROMLE(res + i + 4*j, wdata[j]);
   }

   if (i < n) { 
      salsa20_apply(state, wdata);

      for (j = 0; j < 16; j++)
         FROMLE(buf + 4*j, wdata[j]);

      pos = n-i;
      for (j = 0; j < pos; j++)
         res[i+j] = buf[j];
   }
}

#if defined(NTL_RANDOM_AES256CTR)

/* Size must be a multiple of AES block-size (16 bytes). */
#define BUFSIZE                 4096
//#define BUFSIZE                 8192

static void
inc32(unsigned char ctr[16])
{
   int i, c = 1;

   for (i = 0; i < 4; i++) {
      c += ctr[15 - i];
      ctr[15 - i] = (unsigned char)c;
      c >>= 8;
   }
}

#if defined(NTL_HAVE_AES_NI) && defined(NTL_HAVE_AVX2)

/*****************************************************************
This optimized AES-256 implementation is derived from public
domain code.

Authors:
Romain Dolbeau

Obtained from:
https://github.com/floodyberry/supercop/blob/master/crypto_stream/aes256ctr/dolbeau/aesenc-int/aesenc-int.c
*/

#ifdef __INTEL_COMPILER
#define ALIGN16 __declspec(align(16))
#define ALIGN32 __declspec(align(32))
#define ALIGN64 __declspec(align(64))
#else // assume GCC
#define ALIGN16  __attribute__((aligned(16)))
#define ALIGN32  __attribute__((aligned(32)))
#define ALIGN64  __attribute__((aligned(64)))
#ifndef _bswap64
#define _bswap64(a) __builtin_bswap64(a)
#endif
#ifndef  _bswap
#define _bswap(a) __builtin_bswap(a)
#endif
#endif

static inline void aesni_key256_expand(const unsigned char* key, __m128i rkeys[16]) {
   __m128i key0 = _mm_loadu_si128((const __m128i *)(key+0));
   __m128i key1 = _mm_loadu_si128((const __m128i *)(key+16));
   __m128i temp0, temp1, temp2, temp4;
   int idx = 0;

   rkeys[idx++] = key0;
   temp0 = key0;
   temp2 = key1;

   /* blockshift-based block by Cedric Bourrasset & Romain Dolbeau */
#define BLOCK1(IMM)                              \
   temp1 = _mm_aeskeygenassist_si128(temp2, IMM); \
   rkeys[idx++] = temp2;                          \
   temp4 = _mm_slli_si128(temp0,4);               \
   temp0 = _mm_xor_si128(temp0,temp4);            \
   temp4 = _mm_slli_si128(temp0,8);               \
   temp0 = _mm_xor_si128(temp0,temp4);            \
   temp1 = _mm_shuffle_epi32(temp1,0xff);         \
   temp0 = _mm_xor_si128(temp0,temp1)

#define BLOCK2(IMM)                              \
   temp1 = _mm_aeskeygenassist_si128(temp0, IMM); \
   rkeys[idx++] = temp0;                          \
   temp4 = _mm_slli_si128(temp2,4);               \
   temp2 = _mm_xor_si128(temp2,temp4);            \
   temp4 = _mm_slli_si128(temp2,8);               \
   temp2 = _mm_xor_si128(temp2,temp4);            \
   temp1 = _mm_shuffle_epi32(temp1,0xaa);         \
   temp2 = _mm_xor_si128(temp2,temp1)

   BLOCK1(0x01);
   BLOCK2(0x01);

   BLOCK1(0x02);
   BLOCK2(0x02);

   BLOCK1(0x04);
   BLOCK2(0x04);

   BLOCK1(0x08);
   BLOCK2(0x08);

   BLOCK1(0x10);
   BLOCK2(0x10);

   BLOCK1(0x20);
   BLOCK2(0x20);

   BLOCK1(0x40);
   rkeys[idx++] = temp0;
}

/** single, by-the-book AES encryption with AES-NI */
static inline void aesni_encrypt1(unsigned char *out, unsigned char *n, __m128i rkeys[16]) {
   __m128i nv = _mm_load_si128((const __m128i *)n);
   int i;
   __m128i temp = _mm_xor_si128(nv, rkeys[0]);
#if 0
// This pragma is not recognized by GCC < 8
#pragma unroll(13)
   for (i = 1 ; i < 14 ; i++) {
      temp = _mm_aesenc_si128(temp, rkeys[i]);
   }
#else
   temp = _mm_aesenc_si128(temp, rkeys[ 1]);
   temp = _mm_aesenc_si128(temp, rkeys[ 2]);
   temp = _mm_aesenc_si128(temp, rkeys[ 3]);
   temp = _mm_aesenc_si128(temp, rkeys[ 4]);
   temp = _mm_aesenc_si128(temp, rkeys[ 5]);
   temp = _mm_aesenc_si128(temp, rkeys[ 6]);
   temp = _mm_aesenc_si128(temp, rkeys[ 7]);
   temp = _mm_aesenc_si128(temp, rkeys[ 8]);
   temp = _mm_aesenc_si128(temp, rkeys[ 9]);
   temp = _mm_aesenc_si128(temp, rkeys[10]);
   temp = _mm_aesenc_si128(temp, rkeys[11]);
   temp = _mm_aesenc_si128(temp, rkeys[12]);
   temp = _mm_aesenc_si128(temp, rkeys[13]);
#endif
   temp = _mm_aesenclast_si128(temp, rkeys[14]);
   _mm_store_si128((__m128i*)(out), temp);
}

/** increment the 16-bytes nonce ;
    this really should be improved somehow...
    but it's not yet time-critical, because we
    use the vector variant anyway  */
static inline void incle(unsigned char n[16]) {
/*   unsigned long long out; */
/*   unsigned char carry; */
   unsigned long long *n_ = (unsigned long long*)n;
   n_[1]++;
   if (n_[1] == 0)
      n_[0] ++;
  /* perhaps this will be efficient on broadwell ? */
  /*   carry = _addcarry_u64(0, n_[1], 1ULL, &out); */
  /*   carry = _addcarry_u64(carry, n_[0], 0ULL, &out); */
}

/** multiple-blocks-at-once AES encryption with AES-NI ;
    on Haswell, aesenc as a latency of 7 and a througput of 1
    so the sequence of aesenc should be bubble-free, if you
    have at least 8 blocks. Let's build an arbitratry-sized
    function */
/* Step 1 : loading the nonce */
/* load & increment the n vector (non-vectorized, unused for now) */
#define NVx(a)                                                  \
  __m128i nv##a = _mm_shuffle_epi8(_mm_load_si128((const __m128i *)n), _mm_set_epi8(8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7)); incle(n)
/* load the incremented n vector (vectorized, probably buggy) */
#define NVxV_DEC(a)                                                     \
  __m128i nv##a;
#define NVxV_NOWRAP(a)                                                  \
  nv##a = _mm_shuffle_epi8(_mm_add_epi64(nv0i, _mm_set_epi64x(a,0)), _mm_set_epi8(8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7))
#define NVxV_WRAP(a)                                                    \
  __m128i ad##a = _mm_add_epi64(nv0i, _mm_set_epi64x(a,a>=wrapnumber?1:0)); \
  nv##a = _mm_shuffle_epi8(ad##a, _mm_set_epi8(8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7))

/* Step 2 : define value in round one (xor with subkey #0, aka key) */
#define TEMPx(a)                                        \
  __m128i temp##a = _mm_xor_si128(nv##a, rkeys[0])

/* Step 3: one round of AES */
#define AESENCx(a)                                      \
  temp##a =  _mm_aesenc_si128(temp##a, rkeys[i]);

/* Step 4: last round of AES */
#define AESENCLASTx(a)                                  \
  temp##a = _mm_aesenclast_si128(temp##a, rkeys[14]);

/* Step 5: store result */
#define STOREx(a)                                       \
  _mm_store_si128((__m128i*)(out+(a*16)), temp##a);

/* all the MAKE* macros are for automatic explicit unrolling */
#define MAKE4(X)                                \
  X(0);X(1);X(2);X(3)

#define MAKE6(X)                                \
  X(0);X(1);X(2);X(3);                          \
  X(4);X(5)

#define MAKE7(X)                                \
  X(0);X(1);X(2);X(3);                          \
  X(4);X(5);X(6)

#define MAKE8(X)                                \
  X(0);X(1);X(2);X(3);                          \
  X(4);X(5);X(6);X(7)

#define MAKE10(X)                               \
  X(0);X(1);X(2);X(3);                          \
  X(4);X(5);X(6);X(7);                          \
  X(8);X(9)

#define MAKE12(X)                               \
  X(0);X(1);X(2);X(3);                          \
  X(4);X(5);X(6);X(7);                          \
  X(8);X(9);X(10);X(11)

#define MAKE16(X)                               \
  X(0);X(1);X(2);X(3);                          \
  X(4);X(5);X(6);X(7);                          \
  X(8);X(9);X(10);X(11);                        \
  X(12);X(13);X(14);X(15)

/* create a function of unrolling N ; the MAKEN is the unrolling
   macro, defined above. The N in MAKEN must match N, obviously. */
#define FUNC(N, MAKEN)                          \
  static inline void aesni_encrypt##N(unsigned char *out, unsigned char *n, __m128i rkeys[16]) { \
    __m128i nv0i = _mm_load_si128((const __m128i *)n);                  \
    long long nl = *(long long*)&n[8];                                  \
    MAKEN(NVxV_DEC);                                                    \
    /* check for nonce wraparound */                                    \
    if ((nl < 0) && (nl + N) >= 0) {                                \
      int wrapnumber = (int)(N - (nl+N));                               \
      MAKEN(NVxV_WRAP);                                                 \
      _mm_storeu_si128((__m128i*)n, _mm_add_epi64(nv0i, _mm_set_epi64x(N,1))); \
    } else {                                                            \
      MAKEN(NVxV_NOWRAP);                                               \
      _mm_storeu_si128((__m128i*)n, _mm_add_epi64(nv0i, _mm_set_epi64x(N,0))); \
    }                                                                   \
    int i;                                                              \
    MAKEN(TEMPx);                                                       \
    for (i = 1 ; i < 14 ; i++) {                                        \
      MAKEN(AESENCx);                                                   \
    }                                                                   \
    MAKEN(AESENCLASTx);                                                 \
    MAKEN(STOREx);                                                      \
  }

/* and now building our unrolled function is trivial */
FUNC(4, MAKE4)
FUNC(6, MAKE6)
FUNC(7, MAKE7)
FUNC(8, MAKE8)
FUNC(10, MAKE10)
FUNC(12, MAKE12)
FUNC(16, MAKE16)

void crypto_stream(
unsigned char *out,
unsigned long long outlen,
unsigned char *n,
const unsigned char *k
)
{
   __m128i rkeys[16];
   ALIGN16 unsigned char n2[16];
   unsigned long long i, j;
   aesni_key256_expand(k, rkeys);
   /* n2 is in byte-reversed (i.e., native little endian)
      order to make increment/testing easier */
   (*(unsigned long long*)&n2[8]) = _bswap64((*(unsigned long long*)&n[8]));
   (*(unsigned long long*)&n2[0]) = _bswap64((*(unsigned long long*)&n[0]));

#define LOOP(iter)                                          \
   int lb = iter * 16;                                      \
   for (i = 0 ; i < outlen ; i+= lb) {                      \
      ALIGN16 unsigned char outni[lb];                      \
      aesni_encrypt##iter(outni, n2, rkeys);                \
      unsigned long long mj = lb;                           \
      if ((i+mj)>=outlen)                                   \
         mj = outlen-i;                                     \
      for (j = 0 ; j < mj ; j++)                            \
         out[i+j] = outni[j];                               \
   }

   LOOP(8);

   (*(unsigned long long*)&n[8]) = _bswap64((*(unsigned long long*)&n2[8]));
   (*(unsigned long long*)&n[0]) = _bswap64((*(unsigned long long*)&n2[0]));
}

static void
aes256ctr_stream(unsigned char out[BUFSIZE], unsigned char iv[16], const unsigned char key[32])
{
   crypto_stream(out, BUFSIZE, iv, key);
}

/*****************************************************************/

#elif defined(NTL_HAVE_KMA)

static void
aes256ctr_stream(unsigned char out[BUFSIZE], unsigned char iv[16], const unsigned char key[32])
{
   static const unsigned char zerobuf[BUFSIZE] = {0};
   unsigned long fc = CPACF_KMA_GCM_AES_256 | CPACF_KMA_HS | CPACF_KMA_LAAD;
   struct {
      unsigned char reserved[12];
      unsigned int cv;
      unsigned char _[48];
      unsigned char j0[16];
      unsigned char k[32];
   } param;

   memcpy(&param.cv, &iv[12], sizeof(param.cv));
   param.cv--;
   memcpy(&param.j0[0], &iv[0], sizeof(param.j0) - sizeof(param.cv));
   memcpy(&param.j0[12], &param.cv, sizeof(param.cv));
   memcpy(param.k, key, sizeof(param.k));

   cpacf_kma(fc, &param, out, NULL, 0, zerobuf, sizeof(zerobuf));

   param.cv++;
   memcpy(&iv[12], &param.cv, sizeof(param.cv));
}

#else

/*****************************************************************
This AES-256 reference implementation is derived from
public domain code.

Authors:
Vincent Rijmen <vincent.rijmen@esat.kuleuven.ac.be>
Antoon Bosselaers <antoon.bosselaers@esat.kuleuven.ac.be>
Paulo Barreto <paulo.barreto@terra.com.br>

Obtained from:
https://github.com/zakird/zdlibc/blob/master/rijndael-alg-fst.c
*/

typedef uint8_t u8;
typedef uint32_t u32;

static const u32 Te0[256] = {
   0xc66363a5U, 0xf87c7c84U, 0xee777799U, 0xf67b7b8dU,
   0xfff2f20dU, 0xd66b6bbdU, 0xde6f6fb1U, 0x91c5c554U,
   0x60303050U, 0x02010103U, 0xce6767a9U, 0x562b2b7dU,
   0xe7fefe19U, 0xb5d7d762U, 0x4dababe6U, 0xec76769aU,
   0x8fcaca45U, 0x1f82829dU, 0x89c9c940U, 0xfa7d7d87U,
   0xeffafa15U, 0xb25959ebU, 0x8e4747c9U, 0xfbf0f00bU,
   0x41adadecU, 0xb3d4d467U, 0x5fa2a2fdU, 0x45afafeaU,
   0x239c9cbfU, 0x53a4a4f7U, 0xe4727296U, 0x9bc0c05bU,
   0x75b7b7c2U, 0xe1fdfd1cU, 0x3d9393aeU, 0x4c26266aU,
   0x6c36365aU, 0x7e3f3f41U, 0xf5f7f702U, 0x83cccc4fU,
   0x6834345cU, 0x51a5a5f4U, 0xd1e5e534U, 0xf9f1f108U,
   0xe2717193U, 0xabd8d873U, 0x62313153U, 0x2a15153fU,
   0x0804040cU, 0x95c7c752U, 0x46232365U, 0x9dc3c35eU,
   0x30181828U, 0x379696a1U, 0x0a05050fU, 0x2f9a9ab5U,
   0x0e070709U, 0x24121236U, 0x1b80809bU, 0xdfe2e23dU,
   0xcdebeb26U, 0x4e272769U, 0x7fb2b2cdU, 0xea75759fU,
   0x1209091bU, 0x1d83839eU, 0x582c2c74U, 0x341a1a2eU,
   0x361b1b2dU, 0xdc6e6eb2U, 0xb45a5aeeU, 0x5ba0a0fbU,
   0xa45252f6U, 0x763b3b4dU, 0xb7d6d661U, 0x7db3b3ceU,
   0x5229297bU, 0xdde3e33eU, 0x5e2f2f71U, 0x13848497U,
   0xa65353f5U, 0xb9d1d168U, 0x00000000U, 0xc1eded2cU,
   0x40202060U, 0xe3fcfc1fU, 0x79b1b1c8U, 0xb65b5bedU,
   0xd46a6abeU, 0x8dcbcb46U, 0x67bebed9U, 0x7239394bU,
   0x944a4adeU, 0x984c4cd4U, 0xb05858e8U, 0x85cfcf4aU,
   0xbbd0d06bU, 0xc5efef2aU, 0x4faaaae5U, 0xedfbfb16U,
   0x864343c5U, 0x9a4d4dd7U, 0x66333355U, 0x11858594U,
   0x8a4545cfU, 0xe9f9f910U, 0x04020206U, 0xfe7f7f81U,
   0xa05050f0U, 0x783c3c44U, 0x259f9fbaU, 0x4ba8a8e3U,
   0xa25151f3U, 0x5da3a3feU, 0x804040c0U, 0x058f8f8aU,
   0x3f9292adU, 0x219d9dbcU, 0x70383848U, 0xf1f5f504U,
   0x63bcbcdfU, 0x77b6b6c1U, 0xafdada75U, 0x42212163U,
   0x20101030U, 0xe5ffff1aU, 0xfdf3f30eU, 0xbfd2d26dU,
   0x81cdcd4cU, 0x180c0c14U, 0x26131335U, 0xc3ecec2fU,
   0xbe5f5fe1U, 0x359797a2U, 0x884444ccU, 0x2e171739U,
   0x93c4c457U, 0x55a7a7f2U, 0xfc7e7e82U, 0x7a3d3d47U,
   0xc86464acU, 0xba5d5de7U, 0x3219192bU, 0xe6737395U,
   0xc06060a0U, 0x19818198U, 0x9e4f4fd1U, 0xa3dcdc7fU,
   0x44222266U, 0x542a2a7eU, 0x3b9090abU, 0x0b888883U,
   0x8c4646caU, 0xc7eeee29U, 0x6bb8b8d3U, 0x2814143cU,
   0xa7dede79U, 0xbc5e5ee2U, 0x160b0b1dU, 0xaddbdb76U,
   0xdbe0e03bU, 0x64323256U, 0x743a3a4eU, 0x140a0a1eU,
   0x924949dbU, 0x0c06060aU, 0x4824246cU, 0xb85c5ce4U,
   0x9fc2c25dU, 0xbdd3d36eU, 0x43acacefU, 0xc46262a6U,
   0x399191a8U, 0x319595a4U, 0xd3e4e437U, 0xf279798bU,
   0xd5e7e732U, 0x8bc8c843U, 0x6e373759U, 0xda6d6db7U,
   0x018d8d8cU, 0xb1d5d564U, 0x9c4e4ed2U, 0x49a9a9e0U,
   0xd86c6cb4U, 0xac5656faU, 0xf3f4f407U, 0xcfeaea25U,
   0xca6565afU, 0xf47a7a8eU, 0x47aeaee9U, 0x10080818U,
   0x6fbabad5U, 0xf0787888U, 0x4a25256fU, 0x5c2e2e72U,
   0x381c1c24U, 0x57a6a6f1U, 0x73b4b4c7U, 0x97c6c651U,
   0xcbe8e823U, 0xa1dddd7cU, 0xe874749cU, 0x3e1f1f21U,
   0x964b4bddU, 0x61bdbddcU, 0x0d8b8b86U, 0x0f8a8a85U,
   0xe0707090U, 0x7c3e3e42U, 0x71b5b5c4U, 0xcc6666aaU,
   0x904848d8U, 0x06030305U, 0xf7f6f601U, 0x1c0e0e12U,
   0xc26161a3U, 0x6a35355fU, 0xae5757f9U, 0x69b9b9d0U,
   0x17868691U, 0x99c1c158U, 0x3a1d1d27U, 0x279e9eb9U,
   0xd9e1e138U, 0xebf8f813U, 0x2b9898b3U, 0x22111133U,
   0xd26969bbU, 0xa9d9d970U, 0x078e8e89U, 0x339494a7U,
   0x2d9b9bb6U, 0x3c1e1e22U, 0x15878792U, 0xc9e9e920U,
   0x87cece49U, 0xaa5555ffU, 0x50282878U, 0xa5dfdf7aU,
   0x038c8c8fU, 0x59a1a1f8U, 0x09898980U, 0x1a0d0d17U,
   0x65bfbfdaU, 0xd7e6e631U, 0x844242c6U, 0xd06868b8U,
   0x824141c3U, 0x299999b0U, 0x5a2d2d77U, 0x1e0f0f11U,
   0x7bb0b0cbU, 0xa85454fcU, 0x6dbbbbd6U, 0x2c16163aU,
};
static const u32 Te1[256] = {
   0xa5c66363U, 0x84f87c7cU, 0x99ee7777U, 0x8df67b7bU,
   0x0dfff2f2U, 0xbdd66b6bU, 0xb1de6f6fU, 0x5491c5c5U,
   0x50603030U, 0x03020101U, 0xa9ce6767U, 0x7d562b2bU,
   0x19e7fefeU, 0x62b5d7d7U, 0xe64dababU, 0x9aec7676U,
   0x458fcacaU, 0x9d1f8282U, 0x4089c9c9U, 0x87fa7d7dU,
   0x15effafaU, 0xebb25959U, 0xc98e4747U, 0x0bfbf0f0U,
   0xec41adadU, 0x67b3d4d4U, 0xfd5fa2a2U, 0xea45afafU,
   0xbf239c9cU, 0xf753a4a4U, 0x96e47272U, 0x5b9bc0c0U,
   0xc275b7b7U, 0x1ce1fdfdU, 0xae3d9393U, 0x6a4c2626U,
   0x5a6c3636U, 0x417e3f3fU, 0x02f5f7f7U, 0x4f83ccccU,
   0x5c683434U, 0xf451a5a5U, 0x34d1e5e5U, 0x08f9f1f1U,
   0x93e27171U, 0x73abd8d8U, 0x53623131U, 0x3f2a1515U,
   0x0c080404U, 0x5295c7c7U, 0x65462323U, 0x5e9dc3c3U,
   0x28301818U, 0xa1379696U, 0x0f0a0505U, 0xb52f9a9aU,
   0x090e0707U, 0x36241212U, 0x9b1b8080U, 0x3ddfe2e2U,
   0x26cdebebU, 0x694e2727U, 0xcd7fb2b2U, 0x9fea7575U,
   0x1b120909U, 0x9e1d8383U, 0x74582c2cU, 0x2e341a1aU,
   0x2d361b1bU, 0xb2dc6e6eU, 0xeeb45a5aU, 0xfb5ba0a0U,
   0xf6a45252U, 0x4d763b3bU, 0x61b7d6d6U, 0xce7db3b3U,
   0x7b522929U, 0x3edde3e3U, 0x715e2f2fU, 0x97138484U,
   0xf5a65353U, 0x68b9d1d1U, 0x00000000U, 0x2cc1ededU,
   0x60402020U, 0x1fe3fcfcU, 0xc879b1b1U, 0xedb65b5bU,
   0xbed46a6aU, 0x468dcbcbU, 0xd967bebeU, 0x4b723939U,
   0xde944a4aU, 0xd4984c4cU, 0xe8b05858U, 0x4a85cfcfU,
   0x6bbbd0d0U, 0x2ac5efefU, 0xe54faaaaU, 0x16edfbfbU,
   0xc5864343U, 0xd79a4d4dU, 0x55663333U, 0x94118585U,
   0xcf8a4545U, 0x10e9f9f9U, 0x06040202U, 0x81fe7f7fU,
   0xf0a05050U, 0x44783c3cU, 0xba259f9fU, 0xe34ba8a8U,
   0xf3a25151U, 0xfe5da3a3U, 0xc0804040U, 0x8a058f8fU,
   0xad3f9292U, 0xbc219d9dU, 0x48703838U, 0x04f1f5f5U,
   0xdf63bcbcU, 0xc177b6b6U, 0x75afdadaU, 0x63422121U,
   0x30201010U, 0x1ae5ffffU, 0x0efdf3f3U, 0x6dbfd2d2U,
   0x4c81cdcdU, 0x14180c0cU, 0x35261313U, 0x2fc3ececU,
   0xe1be5f5fU, 0xa2359797U, 0xcc884444U, 0x392e1717U,
   0x5793c4c4U, 0xf255a7a7U, 0x82fc7e7eU, 0x477a3d3dU,
   0xacc86464U, 0xe7ba5d5dU, 0x2b321919U, 0x95e67373U,
   0xa0c06060U, 0x98198181U, 0xd19e4f4fU, 0x7fa3dcdcU,
   0x66442222U, 0x7e542a2aU, 0xab3b9090U, 0x830b8888U,
   0xca8c4646U, 0x29c7eeeeU, 0xd36bb8b8U, 0x3c281414U,
   0x79a7dedeU, 0xe2bc5e5eU, 0x1d160b0bU, 0x76addbdbU,
   0x3bdbe0e0U, 0x56643232U, 0x4e743a3aU, 0x1e140a0aU,
   0xdb924949U, 0x0a0c0606U, 0x6c482424U, 0xe4b85c5cU,
   0x5d9fc2c2U, 0x6ebdd3d3U, 0xef43acacU, 0xa6c46262U,
   0xa8399191U, 0xa4319595U, 0x37d3e4e4U, 0x8bf27979U,
   0x32d5e7e7U, 0x438bc8c8U, 0x596e3737U, 0xb7da6d6dU,
   0x8c018d8dU, 0x64b1d5d5U, 0xd29c4e4eU, 0xe049a9a9U,
   0xb4d86c6cU, 0xfaac5656U, 0x07f3f4f4U, 0x25cfeaeaU,
   0xafca6565U, 0x8ef47a7aU, 0xe947aeaeU, 0x18100808U,
   0xd56fbabaU, 0x88f07878U, 0x6f4a2525U, 0x725c2e2eU,
   0x24381c1cU, 0xf157a6a6U, 0xc773b4b4U, 0x5197c6c6U,
   0x23cbe8e8U, 0x7ca1ddddU, 0x9ce87474U, 0x213e1f1fU,
   0xdd964b4bU, 0xdc61bdbdU, 0x860d8b8bU, 0x850f8a8aU,
   0x90e07070U, 0x427c3e3eU, 0xc471b5b5U, 0xaacc6666U,
   0xd8904848U, 0x05060303U, 0x01f7f6f6U, 0x121c0e0eU,
   0xa3c26161U, 0x5f6a3535U, 0xf9ae5757U, 0xd069b9b9U,
   0x91178686U, 0x5899c1c1U, 0x273a1d1dU, 0xb9279e9eU,
   0x38d9e1e1U, 0x13ebf8f8U, 0xb32b9898U, 0x33221111U,
   0xbbd26969U, 0x70a9d9d9U, 0x89078e8eU, 0xa7339494U,
   0xb62d9b9bU, 0x223c1e1eU, 0x92158787U, 0x20c9e9e9U,
   0x4987ceceU, 0xffaa5555U, 0x78502828U, 0x7aa5dfdfU,
   0x8f038c8cU, 0xf859a1a1U, 0x80098989U, 0x171a0d0dU,
   0xda65bfbfU, 0x31d7e6e6U, 0xc6844242U, 0xb8d06868U,
   0xc3824141U, 0xb0299999U, 0x775a2d2dU, 0x111e0f0fU,
   0xcb7bb0b0U, 0xfca85454U, 0xd66dbbbbU, 0x3a2c1616U,
};
static const u32 Te2[256] = {
   0x63a5c663U, 0x7c84f87cU, 0x7799ee77U, 0x7b8df67bU,
   0xf20dfff2U, 0x6bbdd66bU, 0x6fb1de6fU, 0xc55491c5U,
   0x30506030U, 0x01030201U, 0x67a9ce67U, 0x2b7d562bU,
   0xfe19e7feU, 0xd762b5d7U, 0xabe64dabU, 0x769aec76U,
   0xca458fcaU, 0x829d1f82U, 0xc94089c9U, 0x7d87fa7dU,
   0xfa15effaU, 0x59ebb259U, 0x47c98e47U, 0xf00bfbf0U,
   0xadec41adU, 0xd467b3d4U, 0xa2fd5fa2U, 0xafea45afU,
   0x9cbf239cU, 0xa4f753a4U, 0x7296e472U, 0xc05b9bc0U,
   0xb7c275b7U, 0xfd1ce1fdU, 0x93ae3d93U, 0x266a4c26U,
   0x365a6c36U, 0x3f417e3fU, 0xf702f5f7U, 0xcc4f83ccU,
   0x345c6834U, 0xa5f451a5U, 0xe534d1e5U, 0xf108f9f1U,
   0x7193e271U, 0xd873abd8U, 0x31536231U, 0x153f2a15U,
   0x040c0804U, 0xc75295c7U, 0x23654623U, 0xc35e9dc3U,
   0x18283018U, 0x96a13796U, 0x050f0a05U, 0x9ab52f9aU,
   0x07090e07U, 0x12362412U, 0x809b1b80U, 0xe23ddfe2U,
   0xeb26cdebU, 0x27694e27U, 0xb2cd7fb2U, 0x759fea75U,
   0x091b1209U, 0x839e1d83U, 0x2c74582cU, 0x1a2e341aU,
   0x1b2d361bU, 0x6eb2dc6eU, 0x5aeeb45aU, 0xa0fb5ba0U,
   0x52f6a452U, 0x3b4d763bU, 0xd661b7d6U, 0xb3ce7db3U,
   0x297b5229U, 0xe33edde3U, 0x2f715e2fU, 0x84971384U,
   0x53f5a653U, 0xd168b9d1U, 0x00000000U, 0xed2cc1edU,
   0x20604020U, 0xfc1fe3fcU, 0xb1c879b1U, 0x5bedb65bU,
   0x6abed46aU, 0xcb468dcbU, 0xbed967beU, 0x394b7239U,
   0x4ade944aU, 0x4cd4984cU, 0x58e8b058U, 0xcf4a85cfU,
   0xd06bbbd0U, 0xef2ac5efU, 0xaae54faaU, 0xfb16edfbU,
   0x43c58643U, 0x4dd79a4dU, 0x33556633U, 0x85941185U,
   0x45cf8a45U, 0xf910e9f9U, 0x02060402U, 0x7f81fe7fU,
   0x50f0a050U, 0x3c44783cU, 0x9fba259fU, 0xa8e34ba8U,
   0x51f3a251U, 0xa3fe5da3U, 0x40c08040U, 0x8f8a058fU,
   0x92ad3f92U, 0x9dbc219dU, 0x38487038U, 0xf504f1f5U,
   0xbcdf63bcU, 0xb6c177b6U, 0xda75afdaU, 0x21634221U,
   0x10302010U, 0xff1ae5ffU, 0xf30efdf3U, 0xd26dbfd2U,
   0xcd4c81cdU, 0x0c14180cU, 0x13352613U, 0xec2fc3ecU,
   0x5fe1be5fU, 0x97a23597U, 0x44cc8844U, 0x17392e17U,
   0xc45793c4U, 0xa7f255a7U, 0x7e82fc7eU, 0x3d477a3dU,
   0x64acc864U, 0x5de7ba5dU, 0x192b3219U, 0x7395e673U,
   0x60a0c060U, 0x81981981U, 0x4fd19e4fU, 0xdc7fa3dcU,
   0x22664422U, 0x2a7e542aU, 0x90ab3b90U, 0x88830b88U,
   0x46ca8c46U, 0xee29c7eeU, 0xb8d36bb8U, 0x143c2814U,
   0xde79a7deU, 0x5ee2bc5eU, 0x0b1d160bU, 0xdb76addbU,
   0xe03bdbe0U, 0x32566432U, 0x3a4e743aU, 0x0a1e140aU,
   0x49db9249U, 0x060a0c06U, 0x246c4824U, 0x5ce4b85cU,
   0xc25d9fc2U, 0xd36ebdd3U, 0xacef43acU, 0x62a6c462U,
   0x91a83991U, 0x95a43195U, 0xe437d3e4U, 0x798bf279U,
   0xe732d5e7U, 0xc8438bc8U, 0x37596e37U, 0x6db7da6dU,
   0x8d8c018dU, 0xd564b1d5U, 0x4ed29c4eU, 0xa9e049a9U,
   0x6cb4d86cU, 0x56faac56U, 0xf407f3f4U, 0xea25cfeaU,
   0x65afca65U, 0x7a8ef47aU, 0xaee947aeU, 0x08181008U,
   0xbad56fbaU, 0x7888f078U, 0x256f4a25U, 0x2e725c2eU,
   0x1c24381cU, 0xa6f157a6U, 0xb4c773b4U, 0xc65197c6U,
   0xe823cbe8U, 0xdd7ca1ddU, 0x749ce874U, 0x1f213e1fU,
   0x4bdd964bU, 0xbddc61bdU, 0x8b860d8bU, 0x8a850f8aU,
   0x7090e070U, 0x3e427c3eU, 0xb5c471b5U, 0x66aacc66U,
   0x48d89048U, 0x03050603U, 0xf601f7f6U, 0x0e121c0eU,
   0x61a3c261U, 0x355f6a35U, 0x57f9ae57U, 0xb9d069b9U,
   0x86911786U, 0xc15899c1U, 0x1d273a1dU, 0x9eb9279eU,
   0xe138d9e1U, 0xf813ebf8U, 0x98b32b98U, 0x11332211U,
   0x69bbd269U, 0xd970a9d9U, 0x8e89078eU, 0x94a73394U,
   0x9bb62d9bU, 0x1e223c1eU, 0x87921587U, 0xe920c9e9U,
   0xce4987ceU, 0x55ffaa55U, 0x28785028U, 0xdf7aa5dfU,
   0x8c8f038cU, 0xa1f859a1U, 0x89800989U, 0x0d171a0dU,
   0xbfda65bfU, 0xe631d7e6U, 0x42c68442U, 0x68b8d068U,
   0x41c38241U, 0x99b02999U, 0x2d775a2dU, 0x0f111e0fU,
   0xb0cb7bb0U, 0x54fca854U, 0xbbd66dbbU, 0x163a2c16U,
};
static const u32 Te3[256] = {
   0x6363a5c6U, 0x7c7c84f8U, 0x777799eeU, 0x7b7b8df6U,
   0xf2f20dffU, 0x6b6bbdd6U, 0x6f6fb1deU, 0xc5c55491U,
   0x30305060U, 0x01010302U, 0x6767a9ceU, 0x2b2b7d56U,
   0xfefe19e7U, 0xd7d762b5U, 0xababe64dU, 0x76769aecU,
   0xcaca458fU, 0x82829d1fU, 0xc9c94089U, 0x7d7d87faU,
   0xfafa15efU, 0x5959ebb2U, 0x4747c98eU, 0xf0f00bfbU,
   0xadadec41U, 0xd4d467b3U, 0xa2a2fd5fU, 0xafafea45U,
   0x9c9cbf23U, 0xa4a4f753U, 0x727296e4U, 0xc0c05b9bU,
   0xb7b7c275U, 0xfdfd1ce1U, 0x9393ae3dU, 0x26266a4cU,
   0x36365a6cU, 0x3f3f417eU, 0xf7f702f5U, 0xcccc4f83U,
   0x34345c68U, 0xa5a5f451U, 0xe5e534d1U, 0xf1f108f9U,
   0x717193e2U, 0xd8d873abU, 0x31315362U, 0x15153f2aU,
   0x04040c08U, 0xc7c75295U, 0x23236546U, 0xc3c35e9dU,
   0x18182830U, 0x9696a137U, 0x05050f0aU, 0x9a9ab52fU,
   0x0707090eU, 0x12123624U, 0x80809b1bU, 0xe2e23ddfU,
   0xebeb26cdU, 0x2727694eU, 0xb2b2cd7fU, 0x75759feaU,
   0x09091b12U, 0x83839e1dU, 0x2c2c7458U, 0x1a1a2e34U,
   0x1b1b2d36U, 0x6e6eb2dcU, 0x5a5aeeb4U, 0xa0a0fb5bU,
   0x5252f6a4U, 0x3b3b4d76U, 0xd6d661b7U, 0xb3b3ce7dU,
   0x29297b52U, 0xe3e33eddU, 0x2f2f715eU, 0x84849713U,
   0x5353f5a6U, 0xd1d168b9U, 0x00000000U, 0xeded2cc1U,
   0x20206040U, 0xfcfc1fe3U, 0xb1b1c879U, 0x5b5bedb6U,
   0x6a6abed4U, 0xcbcb468dU, 0xbebed967U, 0x39394b72U,
   0x4a4ade94U, 0x4c4cd498U, 0x5858e8b0U, 0xcfcf4a85U,
   0xd0d06bbbU, 0xefef2ac5U, 0xaaaae54fU, 0xfbfb16edU,
   0x4343c586U, 0x4d4dd79aU, 0x33335566U, 0x85859411U,
   0x4545cf8aU, 0xf9f910e9U, 0x02020604U, 0x7f7f81feU,
   0x5050f0a0U, 0x3c3c4478U, 0x9f9fba25U, 0xa8a8e34bU,
   0x5151f3a2U, 0xa3a3fe5dU, 0x4040c080U, 0x8f8f8a05U,
   0x9292ad3fU, 0x9d9dbc21U, 0x38384870U, 0xf5f504f1U,
   0xbcbcdf63U, 0xb6b6c177U, 0xdada75afU, 0x21216342U,
   0x10103020U, 0xffff1ae5U, 0xf3f30efdU, 0xd2d26dbfU,
   0xcdcd4c81U, 0x0c0c1418U, 0x13133526U, 0xecec2fc3U,
   0x5f5fe1beU, 0x9797a235U, 0x4444cc88U, 0x1717392eU,
   0xc4c45793U, 0xa7a7f255U, 0x7e7e82fcU, 0x3d3d477aU,
   0x6464acc8U, 0x5d5de7baU, 0x19192b32U, 0x737395e6U,
   0x6060a0c0U, 0x81819819U, 0x4f4fd19eU, 0xdcdc7fa3U,
   0x22226644U, 0x2a2a7e54U, 0x9090ab3bU, 0x8888830bU,
   0x4646ca8cU, 0xeeee29c7U, 0xb8b8d36bU, 0x14143c28U,
   0xdede79a7U, 0x5e5ee2bcU, 0x0b0b1d16U, 0xdbdb76adU,
   0xe0e03bdbU, 0x32325664U, 0x3a3a4e74U, 0x0a0a1e14U,
   0x4949db92U, 0x06060a0cU, 0x24246c48U, 0x5c5ce4b8U,
   0xc2c25d9fU, 0xd3d36ebdU, 0xacacef43U, 0x6262a6c4U,
   0x9191a839U, 0x9595a431U, 0xe4e437d3U, 0x79798bf2U,
   0xe7e732d5U, 0xc8c8438bU, 0x3737596eU, 0x6d6db7daU,
   0x8d8d8c01U, 0xd5d564b1U, 0x4e4ed29cU, 0xa9a9e049U,
   0x6c6cb4d8U, 0x5656faacU, 0xf4f407f3U, 0xeaea25cfU,
   0x6565afcaU, 0x7a7a8ef4U, 0xaeaee947U, 0x08081810U,
   0xbabad56fU, 0x787888f0U, 0x25256f4aU, 0x2e2e725cU,
   0x1c1c2438U, 0xa6a6f157U, 0xb4b4c773U, 0xc6c65197U,
   0xe8e823cbU, 0xdddd7ca1U, 0x74749ce8U, 0x1f1f213eU,
   0x4b4bdd96U, 0xbdbddc61U, 0x8b8b860dU, 0x8a8a850fU,
   0x707090e0U, 0x3e3e427cU, 0xb5b5c471U, 0x6666aaccU,
   0x4848d890U, 0x03030506U, 0xf6f601f7U, 0x0e0e121cU,
   0x6161a3c2U, 0x35355f6aU, 0x5757f9aeU, 0xb9b9d069U,
   0x86869117U, 0xc1c15899U, 0x1d1d273aU, 0x9e9eb927U,
   0xe1e138d9U, 0xf8f813ebU, 0x9898b32bU, 0x11113322U,
   0x6969bbd2U, 0xd9d970a9U, 0x8e8e8907U, 0x9494a733U,
   0x9b9bb62dU, 0x1e1e223cU, 0x87879215U, 0xe9e920c9U,
   0xcece4987U, 0x5555ffaaU, 0x28287850U, 0xdfdf7aa5U,
   0x8c8c8f03U, 0xa1a1f859U, 0x89898009U, 0x0d0d171aU,
   0xbfbfda65U, 0xe6e631d7U, 0x4242c684U, 0x6868b8d0U,
   0x4141c382U, 0x9999b029U, 0x2d2d775aU, 0x0f0f111eU,
   0xb0b0cb7bU, 0x5454fca8U, 0xbbbbd66dU, 0x16163a2cU,
};
static const u32 Te4[256] = {
   0x63636363U, 0x7c7c7c7cU, 0x77777777U, 0x7b7b7b7bU,
   0xf2f2f2f2U, 0x6b6b6b6bU, 0x6f6f6f6fU, 0xc5c5c5c5U,
   0x30303030U, 0x01010101U, 0x67676767U, 0x2b2b2b2bU,
   0xfefefefeU, 0xd7d7d7d7U, 0xababababU, 0x76767676U,
   0xcacacacaU, 0x82828282U, 0xc9c9c9c9U, 0x7d7d7d7dU,
   0xfafafafaU, 0x59595959U, 0x47474747U, 0xf0f0f0f0U,
   0xadadadadU, 0xd4d4d4d4U, 0xa2a2a2a2U, 0xafafafafU,
   0x9c9c9c9cU, 0xa4a4a4a4U, 0x72727272U, 0xc0c0c0c0U,
   0xb7b7b7b7U, 0xfdfdfdfdU, 0x93939393U, 0x26262626U,
   0x36363636U, 0x3f3f3f3fU, 0xf7f7f7f7U, 0xccccccccU,
   0x34343434U, 0xa5a5a5a5U, 0xe5e5e5e5U, 0xf1f1f1f1U,
   0x71717171U, 0xd8d8d8d8U, 0x31313131U, 0x15151515U,
   0x04040404U, 0xc7c7c7c7U, 0x23232323U, 0xc3c3c3c3U,
   0x18181818U, 0x96969696U, 0x05050505U, 0x9a9a9a9aU,
   0x07070707U, 0x12121212U, 0x80808080U, 0xe2e2e2e2U,
   0xebebebebU, 0x27272727U, 0xb2b2b2b2U, 0x75757575U,
   0x09090909U, 0x83838383U, 0x2c2c2c2cU, 0x1a1a1a1aU,
   0x1b1b1b1bU, 0x6e6e6e6eU, 0x5a5a5a5aU, 0xa0a0a0a0U,
   0x52525252U, 0x3b3b3b3bU, 0xd6d6d6d6U, 0xb3b3b3b3U,
   0x29292929U, 0xe3e3e3e3U, 0x2f2f2f2fU, 0x84848484U,
   0x53535353U, 0xd1d1d1d1U, 0x00000000U, 0xededededU,
   0x20202020U, 0xfcfcfcfcU, 0xb1b1b1b1U, 0x5b5b5b5bU,
   0x6a6a6a6aU, 0xcbcbcbcbU, 0xbebebebeU, 0x39393939U,
   0x4a4a4a4aU, 0x4c4c4c4cU, 0x58585858U, 0xcfcfcfcfU,
   0xd0d0d0d0U, 0xefefefefU, 0xaaaaaaaaU, 0xfbfbfbfbU,
   0x43434343U, 0x4d4d4d4dU, 0x33333333U, 0x85858585U,
   0x45454545U, 0xf9f9f9f9U, 0x02020202U, 0x7f7f7f7fU,
   0x50505050U, 0x3c3c3c3cU, 0x9f9f9f9fU, 0xa8a8a8a8U,
   0x51515151U, 0xa3a3a3a3U, 0x40404040U, 0x8f8f8f8fU,
   0x92929292U, 0x9d9d9d9dU, 0x38383838U, 0xf5f5f5f5U,
   0xbcbcbcbcU, 0xb6b6b6b6U, 0xdadadadaU, 0x21212121U,
   0x10101010U, 0xffffffffU, 0xf3f3f3f3U, 0xd2d2d2d2U,
   0xcdcdcdcdU, 0x0c0c0c0cU, 0x13131313U, 0xececececU,
   0x5f5f5f5fU, 0x97979797U, 0x44444444U, 0x17171717U,
   0xc4c4c4c4U, 0xa7a7a7a7U, 0x7e7e7e7eU, 0x3d3d3d3dU,
   0x64646464U, 0x5d5d5d5dU, 0x19191919U, 0x73737373U,
   0x60606060U, 0x81818181U, 0x4f4f4f4fU, 0xdcdcdcdcU,
   0x22222222U, 0x2a2a2a2aU, 0x90909090U, 0x88888888U,
   0x46464646U, 0xeeeeeeeeU, 0xb8b8b8b8U, 0x14141414U,
   0xdedededeU, 0x5e5e5e5eU, 0x0b0b0b0bU, 0xdbdbdbdbU,
   0xe0e0e0e0U, 0x32323232U, 0x3a3a3a3aU, 0x0a0a0a0aU,
   0x49494949U, 0x06060606U, 0x24242424U, 0x5c5c5c5cU,
   0xc2c2c2c2U, 0xd3d3d3d3U, 0xacacacacU, 0x62626262U,
   0x91919191U, 0x95959595U, 0xe4e4e4e4U, 0x79797979U,
   0xe7e7e7e7U, 0xc8c8c8c8U, 0x37373737U, 0x6d6d6d6dU,
   0x8d8d8d8dU, 0xd5d5d5d5U, 0x4e4e4e4eU, 0xa9a9a9a9U,
   0x6c6c6c6cU, 0x56565656U, 0xf4f4f4f4U, 0xeaeaeaeaU,
   0x65656565U, 0x7a7a7a7aU, 0xaeaeaeaeU, 0x08080808U,
   0xbabababaU, 0x78787878U, 0x25252525U, 0x2e2e2e2eU,
   0x1c1c1c1cU, 0xa6a6a6a6U, 0xb4b4b4b4U, 0xc6c6c6c6U,
   0xe8e8e8e8U, 0xddddddddU, 0x74747474U, 0x1f1f1f1fU,
   0x4b4b4b4bU, 0xbdbdbdbdU, 0x8b8b8b8bU, 0x8a8a8a8aU,
   0x70707070U, 0x3e3e3e3eU, 0xb5b5b5b5U, 0x66666666U,
   0x48484848U, 0x03030303U, 0xf6f6f6f6U, 0x0e0e0e0eU,
   0x61616161U, 0x35353535U, 0x57575757U, 0xb9b9b9b9U,
   0x86868686U, 0xc1c1c1c1U, 0x1d1d1d1dU, 0x9e9e9e9eU,
   0xe1e1e1e1U, 0xf8f8f8f8U, 0x98989898U, 0x11111111U,
   0x69696969U, 0xd9d9d9d9U, 0x8e8e8e8eU, 0x94949494U,
   0x9b9b9b9bU, 0x1e1e1e1eU, 0x87878787U, 0xe9e9e9e9U,
   0xcecececeU, 0x55555555U, 0x28282828U, 0xdfdfdfdfU,
   0x8c8c8c8cU, 0xa1a1a1a1U, 0x89898989U, 0x0d0d0d0dU,
   0xbfbfbfbfU, 0xe6e6e6e6U, 0x42424242U, 0x68686868U,
   0x41414141U, 0x99999999U, 0x2d2d2d2dU, 0x0f0f0f0fU,
   0xb0b0b0b0U, 0x54545454U, 0xbbbbbbbbU, 0x16161616U,
};
static const u32 rcon[] = {
   0x01000000, 0x02000000, 0x04000000, 0x08000000,
   0x10000000, 0x20000000, 0x40000000, 0x80000000,
   0x1B000000, 0x36000000,
};

#define SWAP(x) (_lrotl(x, 8) & 0x00ff00ff | _lrotr(x, 8) & 0xff00ff00)

#ifdef _MSC_VER
#define GETU32(p) SWAP(*((u32 *)(p)))
#define PUTU32(ct, st) { *((u32 *)(ct)) = SWAP((st)); }
#else
#define GETU32(pt) (((u32)(pt)[0] << 24) ^ ((u32)(pt)[1] << 16) ^ ((u32)(pt)[2] <<  8) ^ ((u32)(pt)[3]))
#define PUTU32(ct, st) { (ct)[0] = (u8)((st) >> 24); (ct)[1] = (u8)((st) >> 16); (ct)[2] = (u8)((st) >>  8); (ct)[3] = (u8)(st); }
#endif

/**
 * Expand the cipher key into the encryption key schedule.
 */
void AES256KeySetupEnc(u32 rk[60], const u8 cipherKey[32]) {
   int i = 0;
   u32 temp;

   rk[0] = GETU32(cipherKey     );
   rk[1] = GETU32(cipherKey +  4);
   rk[2] = GETU32(cipherKey +  8);
   rk[3] = GETU32(cipherKey + 12);
   rk[4] = GETU32(cipherKey + 16);
   rk[5] = GETU32(cipherKey + 20);
   rk[6] = GETU32(cipherKey + 24);
   rk[7] = GETU32(cipherKey + 28);

   for (;;) {
      temp = rk[ 7];

      rk[ 8] = rk[ 0] ^
               (Te4[(temp >> 16) & 0xff] & 0xff000000) ^
               (Te4[(temp >>  8) & 0xff] & 0x00ff0000) ^
               (Te4[(temp      ) & 0xff] & 0x0000ff00) ^
               (Te4[(temp >> 24)       ] & 0x000000ff) ^
               rcon[i];

      rk[ 9] = rk[ 1] ^ rk[ 8];
      rk[10] = rk[ 2] ^ rk[ 9];
      rk[11] = rk[ 3] ^ rk[10];

      if (++i == 7)
         return;

      temp = rk[11];

      rk[12] = rk[ 4] ^
               (Te4[(temp >> 24)       ] & 0xff000000) ^
               (Te4[(temp >> 16) & 0xff] & 0x00ff0000) ^
               (Te4[(temp >>  8) & 0xff] & 0x0000ff00) ^
               (Te4[(temp      ) & 0xff] & 0x000000ff);

      rk[13] = rk[ 5] ^ rk[12];
      rk[14] = rk[ 6] ^ rk[13];
      rk[15] = rk[ 7] ^ rk[14];

      rk += 8;
   }
}

void AES256Encrypt(const u32 rk[60], const u8 pt[16], u8 ct[16]) {
   u32 s0, s1, s2, s3, t0, t1, t2, t3;
   int r, Nr = 14;

   /*
    * map byte array block to cipher state
    * and add initial round key:
    */
   s0 = GETU32(pt     ) ^ rk[0];
   s1 = GETU32(pt +  4) ^ rk[1];
   s2 = GETU32(pt +  8) ^ rk[2];
   s3 = GETU32(pt + 12) ^ rk[3];

   /*
    * Nr - 1 full rounds:
    */
   r = Nr >> 1;

   for (;;) {
      t0 = Te0[(s0 >> 24)       ] ^
           Te1[(s1 >> 16) & 0xff] ^
           Te2[(s2 >>  8) & 0xff] ^
           Te3[(s3      ) & 0xff] ^
           rk[4];
      t1 = Te0[(s1 >> 24)       ] ^
           Te1[(s2 >> 16) & 0xff] ^
           Te2[(s3 >>  8) & 0xff] ^
           Te3[(s0      ) & 0xff] ^
           rk[5];
      t2 = Te0[(s2 >> 24)       ] ^
           Te1[(s3 >> 16) & 0xff] ^
           Te2[(s0 >>  8) & 0xff] ^
           Te3[(s1      ) & 0xff] ^
           rk[6];
      t3 = Te0[(s3 >> 24)       ] ^
           Te1[(s0 >> 16) & 0xff] ^
           Te2[(s1 >>  8) & 0xff] ^
           Te3[(s2      ) & 0xff] ^
           rk[7];

      rk += 8;

      if (--r == 0)
         break;

      s0 = Te0[(t0 >> 24)       ] ^
           Te1[(t1 >> 16) & 0xff] ^
           Te2[(t2 >>  8) & 0xff] ^
           Te3[(t3      ) & 0xff] ^
           rk[0];
      s1 = Te0[(t1 >> 24)       ] ^
           Te1[(t2 >> 16) & 0xff] ^
           Te2[(t3 >>  8) & 0xff] ^
           Te3[(t0      ) & 0xff] ^
           rk[1];
      s2 = Te0[(t2 >> 24)       ] ^
           Te1[(t3 >> 16) & 0xff] ^
           Te2[(t0 >>  8) & 0xff] ^
           Te3[(t1      ) & 0xff] ^
           rk[2];
      s3 = Te0[(t3 >> 24)       ] ^
           Te1[(t0 >> 16) & 0xff] ^
           Te2[(t1 >>  8) & 0xff] ^
           Te3[(t2      ) & 0xff] ^
           rk[3];
   }
   /*
    * apply last round and
    * map cipher state to byte array block:
    */
   s0 = (Te4[(t0 >> 24)       ] & 0xff000000) ^
        (Te4[(t1 >> 16) & 0xff] & 0x00ff0000) ^
        (Te4[(t2 >>  8) & 0xff] & 0x0000ff00) ^
        (Te4[(t3      ) & 0xff] & 0x000000ff) ^
        rk[0];

   PUTU32(ct     , s0);

   s1 = (Te4[(t1 >> 24)       ] & 0xff000000) ^
        (Te4[(t2 >> 16) & 0xff] & 0x00ff0000) ^
        (Te4[(t3 >>  8) & 0xff] & 0x0000ff00) ^
        (Te4[(t0      ) & 0xff] & 0x000000ff) ^
        rk[1];

   PUTU32(ct +  4, s1);

   s2 = (Te4[(t2 >> 24)       ] & 0xff000000) ^
        (Te4[(t3 >> 16) & 0xff] & 0x00ff0000) ^
        (Te4[(t0 >>  8) & 0xff] & 0x0000ff00) ^
        (Te4[(t1      ) & 0xff] & 0x000000ff) ^
        rk[2];

   PUTU32(ct +  8, s2);

   s3 = (Te4[(t3 >> 24)       ] & 0xff000000) ^
        (Te4[(t0 >> 16) & 0xff] & 0x00ff0000) ^
        (Te4[(t1 >>  8) & 0xff] & 0x0000ff00) ^
        (Te4[(t2      ) & 0xff] & 0x000000ff) ^
        rk[3];

   PUTU32(ct + 12, s3);
}

/*****************************************************************/

static void
aes256ctr_stream(unsigned char out[BUFSIZE], unsigned char iv[16], const unsigned char key[32])
{
   u32 rk[60];
   int i;

   AES256KeySetupEnc(rk, key);

   for (i = 0; i < BUFSIZE; i += 16) {
      AES256Encrypt(rk, iv, out + i);
      inc32(iv);
   }
}

#endif

struct RandomStream_impl {
   unsigned char key[32];
   unsigned char iv[16];
   unsigned char buf[BUFSIZE];

   explicit
   RandomStream_impl(const unsigned char *k)
   {
      memcpy(key, k, sizeof(key));
      memset(iv, 0, sizeof(iv));
      iv[15] = 1; // nonce = 1
   }

   const unsigned char *
   get_buf() const
   {
      return buf + sizeof(key);
   }

   long
   get_buf_len() const
   {
      return sizeof(buf) - sizeof(key);
   }

   long get_bytes(unsigned char *res, long n, long pos)
   {
      size_t len;

      if (n < 0)
         LogicError("RandomStream::get: bad args");

      if (n > 0 && sizeof(buf) - sizeof(key) - pos > 0) {
         len = min((size_t)n, sizeof(buf) - sizeof(key) - pos);
         memcpy(res, buf + sizeof(key) + pos, len);

         n -= len;
         res += len;
         pos += len;
      }

      while (n > 0) {
         aes256ctr_stream(buf, iv, key);
         memcpy(key, buf, sizeof(key));

         len = min((size_t)n, sizeof(buf) - sizeof(key));
         memcpy(res, buf + sizeof(key), len);

	 n -= len;
         res += len;
         pos = len;
      }

      return pos;
   }

   void set_nonce(unsigned long nonce)
   {
      // low-order  8 bytes of iv set to zero
      // high-order 8 bytes of iv set to nonce
      memset(iv, 0, sizeof(iv));
      iv[ 8] = (unsigned char) nonce; nonce >>= 8;
      iv[ 9] = (unsigned char) nonce; nonce >>= 8;
      iv[10] = (unsigned char) nonce; nonce >>= 8;
      iv[11] = (unsigned char) nonce; nonce >>= 8;
      iv[12] = (unsigned char) nonce; nonce >>= 8;
      iv[13] = (unsigned char) nonce; nonce >>= 8;
      iv[14] = (unsigned char) nonce; nonce >>= 8;
      iv[15] = (unsigned char) nonce; nonce >>= 8;
   }
};

#else // defined(NTL_RANDOM_AES256CTR)

#if (defined(NTL_HAVE_AVX2) || defined(NTL_HAVE_SSSE3))


/*****************************************************************

This AVX2 implementation is derived from public domain code 
originally developed by Martin Goll Shay Gueron, and obtained from 
here:

https://github.com/floodyberry/supercop/tree/master/crypto_stream/chacha20/goll_gueron

On a Haswell machine, ths code is about 4.x faster than the vanilla 
C code.

The following is the README from that page

==================================================================

This code implements Daniel J. Bernstein's ChaCha stream cipher in C,
targeting architectures with AVX2 and future AVX512 vector extensions.

The implementation improves the slightly modified implementations of Ted Krovetz in the Chromium Project
(http://src.chromium.org/viewvc/chrome/trunk/deps/third_party/nss/nss/lib/freebl/chacha20/chacha20_vec.c and
http://src.chromium.org/viewvc/chrome/trunk/deps/third_party/openssl/openssl/crypto/chacha/chacha_vec.c)
by using the Advanced Vector Extensions AVX2 and, if available in future, AVX512 to widen the vectorization
to 256-bit, respectively 512-bit.

On Intel's Haswell architecture this implementation (using AVX2) is almost ~2x faster than the fastest 
implementation here, when encrypting (decrypting) 2 blocks and more. Also, this implementation is expected 
to double the speed again, when encrypting (decrypting) 4 blocks and more, running on a future architecture
with support for AVX512.

Further details and our measurement results are provided in:
Goll, M., and Gueron,S.: Vectorization of ChaCha Stream Cipher. Cryptology ePrint Archive, 
Report 2013/759, November, 2013, http://eprint.iacr.org/2013/759.pdf

Developers and authors:
*********************************************************
Martin Goll (1) and Shay Gueron (2, 3), 
(1) Ruhr-University Bochum, Germany
(2) University of Haifa, Israel
(3) Intel Corporation, Israel Development Center, Haifa, Israel
*********************************************************

Intellectual Property Notices
-----------------------------

There are no known present or future claims by a copyright holder that the
distribution of this software infringes the copyright. In particular, the author
of the software is not making such claims and does not intend to make such
claims.

There are no known present or future claims by a patent holder that the use of
this software infringes the patent. In particular, the author of the software is
not making such claims and does not intend to make such claims.

Our implementation is in public domain.

*****************************************************************/


// round selector, specified values:
//  8:  low security - high speed
// 12:  mid security -  mid speed
// 20: high security -  low speed
#ifndef CHACHA_RNDS
#define CHACHA_RNDS 20
#endif


#if (defined(NTL_HAVE_AVX2))

typedef __m256i ivec_t;

#define DELTA	_mm256_set_epi64x(0,2,0,2)
#define START   _mm256_set_epi64x(0,1,0,0)
#define NONCE(nonce) _mm256_set_epi64x(nonce, 1, nonce, 0)   

#define STOREU_VEC(m,r)	_mm256_storeu_si256((__m256i*)(m), r)
#define STORE_VEC(m,r)	_mm256_store_si256((__m256i*)(m), r)

#define LOAD_VEC(r,m) r = _mm256_load_si256((const __m256i *)(m))
#define LOADU_VEC(r,m) r = _mm256_loadu_si256((const __m256i *)(m))

#define LOADU_VEC_128(r, m) r = _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i*)(m)))



#define ADD_VEC_32(a,b)	_mm256_add_epi32(a, b)
#define ADD_VEC_64(a,b)	_mm256_add_epi64(a, b)
#define XOR_VEC(a,b)	_mm256_xor_si256(a, b)


#define ROR_VEC_V1(x)	_mm256_shuffle_epi32(x,_MM_SHUFFLE(0,3,2,1))
#define ROR_VEC_V2(x)	_mm256_shuffle_epi32(x,_MM_SHUFFLE(1,0,3,2))
#define ROR_VEC_V3(x)	_mm256_shuffle_epi32(x,_MM_SHUFFLE(2,1,0,3))
#define ROL_VEC_7(x)	XOR_VEC(_mm256_slli_epi32(x, 7), _mm256_srli_epi32(x,25))
#define ROL_VEC_12(x)	XOR_VEC(_mm256_slli_epi32(x,12), _mm256_srli_epi32(x,20))

#define ROL_VEC_8(x)	_mm256_shuffle_epi8(x,_mm256_set_epi8(14,13,12,15,10,9,8,11,6,5,4,7,2,1,0,3,14,13,12,15,10,9,8,11,6,5,4,7,2,1,0,3))


#define ROL_VEC_16(x)	_mm256_shuffle_epi8(x,_mm256_set_epi8(13,12,15,14,9,8,11,10,5,4,7,6,1,0,3,2,13,12,15,14,9,8,11,10,5,4,7,6,1,0,3,2))



#define WRITEU_VEC(op, d, v0, v1, v2, v3)						\
    STOREU_VEC(op + (d + 0*4), _mm256_permute2x128_si256(v0, v1, 0x20));	\
    STOREU_VEC(op + (d + 8*4), _mm256_permute2x128_si256(v2, v3, 0x20));	\
    STOREU_VEC(op + (d +16*4), _mm256_permute2x128_si256(v0, v1, 0x31));	\
    STOREU_VEC(op + (d +24*4), _mm256_permute2x128_si256(v2, v3, 0x31));

#define WRITE_VEC(op, d, v0, v1, v2, v3)						\
    STORE_VEC(op + (d + 0*4), _mm256_permute2x128_si256(v0, v1, 0x20));	\
    STORE_VEC(op + (d + 8*4), _mm256_permute2x128_si256(v2, v3, 0x20));	\
    STORE_VEC(op + (d +16*4), _mm256_permute2x128_si256(v0, v1, 0x31));	\
    STORE_VEC(op + (d +24*4), _mm256_permute2x128_si256(v2, v3, 0x31));

#define SZ_VEC (32)

#define RANSTREAM_NCHUNKS (2)
// leads to a BUFSZ of 512


#elif defined(NTL_HAVE_SSSE3)

typedef __m128i ivec_t;

#define DELTA	_mm_set_epi32(0,0,0,1)
#define START   _mm_setzero_si128()
#define NONCE(nonce) _mm_set_epi64x(nonce,0)

#define STOREU_VEC(m,r)	_mm_storeu_si128((__m128i*)(m), r)
#define STORE_VEC(m,r)	_mm_store_si128((__m128i*)(m), r)

#define LOAD_VEC(r,m) r = _mm_load_si128((const __m128i *)(m))
#define LOADU_VEC(r,m) r = _mm_loadu_si128((const __m128i *)(m))

#define LOADU_VEC_128(r, m) r = _mm_loadu_si128((const __m128i*)(m))

#define ADD_VEC_32(a,b)	_mm_add_epi32(a, b)
#define ADD_VEC_64(a,b)	_mm_add_epi64(a, b)
#define XOR_VEC(a,b)	_mm_xor_si128(a, b)


#define ROR_VEC_V1(x)	_mm_shuffle_epi32(x,_MM_SHUFFLE(0,3,2,1))
#define ROR_VEC_V2(x)	_mm_shuffle_epi32(x,_MM_SHUFFLE(1,0,3,2))
#define ROR_VEC_V3(x)	_mm_shuffle_epi32(x,_MM_SHUFFLE(2,1,0,3))
#define ROL_VEC_7(x)	XOR_VEC(_mm_slli_epi32(x, 7), _mm_srli_epi32(x,25))
#define ROL_VEC_12(x)	XOR_VEC(_mm_slli_epi32(x,12), _mm_srli_epi32(x,20))

#define ROL_VEC_8(x)	_mm_shuffle_epi8(x,_mm_set_epi8(14,13,12,15,10,9,8,11,6,5,4,7,2,1,0,3))


#define ROL_VEC_16(x)	_mm_shuffle_epi8(x,_mm_set_epi8(13,12,15,14,9,8,11,10,5,4,7,6,1,0,3,2))


#define WRITEU_VEC(op, d, v0, v1, v2, v3)	\
    STOREU_VEC(op + (d + 0*4), v0);	\
    STOREU_VEC(op + (d + 4*4), v1);	\
    STOREU_VEC(op + (d + 8*4), v2);	\
    STOREU_VEC(op + (d +12*4), v3);

#define WRITE_VEC(op, d, v0, v1, v2, v3)	\
    STORE_VEC(op + (d + 0*4), v0);	\
    STORE_VEC(op + (d + 4*4), v1);	\
    STORE_VEC(op + (d + 8*4), v2);	\
    STORE_VEC(op + (d +12*4), v3);

#define SZ_VEC (16)

#define RANSTREAM_NCHUNKS (4)
// leads to a BUFSZ of 512

#else

#error "unsupported architecture"

#endif


#define DQROUND_VECTORS_VEC(a,b,c,d)				\
    a = ADD_VEC_32(a,b); d = XOR_VEC(d,a); d = ROL_VEC_16(d);	\
    c = ADD_VEC_32(c,d); b = XOR_VEC(b,c); b = ROL_VEC_12(b);	\
    a = ADD_VEC_32(a,b); d = XOR_VEC(d,a); d = ROL_VEC_8(d);	\
    c = ADD_VEC_32(c,d); b = XOR_VEC(b,c); b = ROL_VEC_7(b);	\
    b = ROR_VEC_V1(b); c = ROR_VEC_V2(c); d = ROR_VEC_V3(d);	\
    a = ADD_VEC_32(a,b); d = XOR_VEC(d,a); d = ROL_VEC_16(d);	\
    c = ADD_VEC_32(c,d); b = XOR_VEC(b,c); b = ROL_VEC_12(b);	\
    a = ADD_VEC_32(a,b); d = XOR_VEC(d,a); d = ROL_VEC_8(d);	\
    c = ADD_VEC_32(c,d); b = XOR_VEC(b,c); b = ROL_VEC_7(b);	\
    b = ROR_VEC_V3(b); c = ROR_VEC_V2(c); d = ROR_VEC_V1(d);



#define RANSTREAM_STATESZ (4*SZ_VEC)

#define RANSTREAM_CHUNKSZ (2*RANSTREAM_STATESZ)
#define RANSTREAM_BUFSZ   (RANSTREAM_NCHUNKS*RANSTREAM_CHUNKSZ)


struct RandomStream_impl {

   AlignedArray<unsigned char> state_store;
   AlignedArray<unsigned char> buf_store;
   long chunk_count;

   void allocate_space() 
   {
      state_store.SetLength(RANSTREAM_STATESZ);
      buf_store.SetLength(RANSTREAM_BUFSZ);
   }


   explicit
   RandomStream_impl(const unsigned char *key) 
   {
      allocate_space();

      unsigned char *state = state_store.elts();

      unsigned int chacha_const[] = {
	      0x61707865,0x3320646E,0x79622D32,0x6B206574
      };


      ivec_t d0, d1, d2, d3;
      LOADU_VEC_128(d0, chacha_const);
      LOADU_VEC_128(d1, key);
      LOADU_VEC_128(d2, key+16);

      d3 = START;


      STORE_VEC(state + 0*SZ_VEC, d0); 
      STORE_VEC(state + 1*SZ_VEC, d1); 
      STORE_VEC(state + 2*SZ_VEC, d2); 
      STORE_VEC(state + 3*SZ_VEC, d3); 

      chunk_count = 0;
   }

   RandomStream_impl(const RandomStream_impl& other) 
   {
      allocate_space();
      *this = other;
   }

   RandomStream_impl& operator=(const RandomStream_impl& other) 
   {
      std::memcpy(state_store.elts(), other.state_store.elts(), RANSTREAM_STATESZ);
      std::memcpy(buf_store.elts(), other.buf_store.elts(), RANSTREAM_BUFSZ);
      chunk_count = other.chunk_count;
      return *this;
   }

   const unsigned char *
   get_buf() const
   {
      return buf_store.elts(); 
   }

   long
   get_buf_len() const
   {
      return RANSTREAM_BUFSZ;
   }

   // bytes are generated in chunks of RANSTREAM_BUFSZ bytes, except that
   // initially, we may generate a few chunks of RANSTREAM_CHUNKSZ
   // bytes.  This optimizes a bit for short bursts following a reset.

   long
   get_bytes(unsigned char *NTL_RESTRICT res, 
             long n, long pos)
   {
      if (n < 0) LogicError("RandomStream::get: bad args");
      if (n == 0) return pos;

      unsigned char *NTL_RESTRICT buf = buf_store.elts();

      if (n <= RANSTREAM_BUFSZ-pos) {
	 std::memcpy(&res[0], &buf[pos], n);
	 pos += n;
	 return pos;
      }

      unsigned char *NTL_RESTRICT state = state_store.elts();

      ivec_t d0, d1, d2, d3;
      LOAD_VEC(d0, state + 0*SZ_VEC);
      LOAD_VEC(d1, state + 1*SZ_VEC);
      LOAD_VEC(d2, state + 2*SZ_VEC);
      LOAD_VEC(d3, state + 3*SZ_VEC);


      // read remainder of buffer
      std::memcpy(&res[0], &buf[pos], RANSTREAM_BUFSZ-pos);
      n -= RANSTREAM_BUFSZ-pos;
      res += RANSTREAM_BUFSZ-pos;
      pos = RANSTREAM_BUFSZ;

      long i = 0;
      for (;  i <= n-RANSTREAM_BUFSZ; i += RANSTREAM_BUFSZ) {

         chunk_count |= RANSTREAM_NCHUNKS;  // disable small buffer strategy

	 for (long j = 0; j < RANSTREAM_NCHUNKS; j++) {
	    ivec_t v0=d0, v1=d1, v2=d2, v3=d3;
	    ivec_t v4=d0, v5=d1, v6=d2, v7=ADD_VEC_64(d3, DELTA);

	    for (long k = 0; k < CHACHA_RNDS/2; k++) {
		    DQROUND_VECTORS_VEC(v0,v1,v2,v3)
		    DQROUND_VECTORS_VEC(v4,v5,v6,v7)
	    }

	    WRITEU_VEC(res+i+j*(8*SZ_VEC), 0, ADD_VEC_32(v0,d0), ADD_VEC_32(v1,d1), ADD_VEC_32(v2,d2), ADD_VEC_32(v3,d3))
	    d3 = ADD_VEC_64(d3, DELTA);
	    WRITEU_VEC(res+i+j*(8*SZ_VEC), 4*SZ_VEC, ADD_VEC_32(v4,d0), ADD_VEC_32(v5,d1), ADD_VEC_32(v6,d2), ADD_VEC_32(v7,d3))
	    d3 = ADD_VEC_64(d3, DELTA);

	 }

      }

      if (i < n) {

         long nchunks;

         if (chunk_count < RANSTREAM_NCHUNKS) {
            nchunks = long(cast_unsigned((n-i)+RANSTREAM_CHUNKSZ-1)/RANSTREAM_CHUNKSZ);
            chunk_count += nchunks;
         }
         else
            nchunks = RANSTREAM_NCHUNKS;

         long pos_offset = RANSTREAM_BUFSZ - nchunks*RANSTREAM_CHUNKSZ;
         buf += pos_offset;

	 for (long j = 0; j < nchunks; j++) {
	    ivec_t v0=d0, v1=d1, v2=d2, v3=d3;
	    ivec_t v4=d0, v5=d1, v6=d2, v7=ADD_VEC_64(d3, DELTA);

	    for (long k = 0; k < CHACHA_RNDS/2; k++) {
               DQROUND_VECTORS_VEC(v0,v1,v2,v3)
               DQROUND_VECTORS_VEC(v4,v5,v6,v7)
	    }

	    WRITE_VEC(buf+j*(8*SZ_VEC), 0, ADD_VEC_32(v0,d0), ADD_VEC_32(v1,d1), ADD_VEC_32(v2,d2), ADD_VEC_32(v3,d3))
	    d3 = ADD_VEC_64(d3, DELTA);
	    WRITE_VEC(buf+j*(8*SZ_VEC), 4*SZ_VEC, ADD_VEC_32(v4,d0), ADD_VEC_32(v5,d1), ADD_VEC_32(v6,d2), ADD_VEC_32(v7,d3))
	    d3 = ADD_VEC_64(d3, DELTA);
	 }

	 pos = n-i+pos_offset;
	 std::memcpy(&res[i], &buf[0], n-i);
      }

      STORE_VEC(state + 3*SZ_VEC, d3); 

      return pos;
   }

   void set_nonce(unsigned long nonce)
   {
      unsigned char *state = state_store.elts();
      ivec_t d3;
      d3 = NONCE(nonce);
      STORE_VEC(state + 3*SZ_VEC, d3);
      chunk_count = 0;
   }

};


#else

struct RandomStream_impl {
   _ntl_uint32 state[16];
   unsigned char buf[64];

   explicit
   RandomStream_impl(const unsigned char *key)
   {
      salsa20_init(state, key);
   }

   const unsigned char *
   get_buf() const
   {
      return &buf[0];
   }

   long
   get_buf_len() const
   {
      return 64;
   }

   long get_bytes(unsigned char *res, long n, long pos) 
   {
      if (n < 0) LogicError("RandomStream::get: bad args");

      long i, j;

      if (n <= 64-pos) {
	 for (i = 0; i < n; i++) res[i] = buf[pos+i];
	 pos += n;
	 return pos;
      }

      // read remainder of buffer
      for (i = 0; i < 64-pos; i++) res[i] = buf[pos+i];
      n -= 64-pos;
      res += 64-pos;
      pos = 64;

      _ntl_uint32 wdata[16];

      // read 64-byte chunks
      for (i = 0; i <= n-64; i += 64) {
	 salsa20_apply(state, wdata);
	 for (j = 0; j < 16; j++)
	    FROMLE(res + i + 4*j, wdata[j]);
      }

      if (i < n) { 
	 salsa20_apply(state, wdata);

	 for (j = 0; j < 16; j++)
	    FROMLE(buf + 4*j, wdata[j]);

	 pos = n-i;
	 for (j = 0; j < pos; j++)
	    res[i+j] = buf[j];
      }

      return pos;
   }

   void set_nonce(unsigned long nonce)
   {
      _ntl_uint32 nonce0, nonce1;

      nonce0 = nonce;
      nonce0 = INT32MASK(nonce0);

      nonce1 = 0;

#if (NTL_BITS_PER_LONG > 32)
      nonce1 = nonce >> 32;
      nonce1 = INT32MASK(nonce1);
#endif

      state[12] = 0;
      state[13] = 0;
      state[14] = nonce0;
      state[15] = nonce1;
   }
};


#endif
#endif // defined(NTL_RANDOMSTREAM_AES256CTR)



// Boilerplate PIMPL code

RandomStream_impl *
RandomStream_impl_build(const unsigned char *key)
{
   UniquePtr<RandomStream_impl> p;
   p.make(key);
   return p.release();
}

RandomStream_impl *
RandomStream_impl_build(const RandomStream_impl& other)
{
   UniquePtr<RandomStream_impl> p;
   p.make(other);
   return p.release();
}

void
RandomStream_impl_copy(RandomStream_impl& x, const RandomStream_impl& y)
{
   x = y;
}

void
RandomStream_impl_delete(RandomStream_impl* p)
{
   delete p;
}

const unsigned char *
RandomStream_impl_get_buf(const RandomStream_impl& x)
{
   return x.get_buf();
}

long
RandomStream_impl_get_buf_len(const RandomStream_impl& x)
{
   return x.get_buf_len();
}

long
RandomStream_impl_get_bytes(RandomStream_impl& impl, 
   unsigned char *res, long n, long pos)
{
   return impl.get_bytes(res, n, pos);
}


void 
RandomStream_impl_set_nonce(RandomStream_impl& impl, unsigned long nonce)
{
   impl.set_nonce(nonce);
}





NTL_TLS_GLOBAL_DECL(UniquePtr<RandomStream>,  CurrentRandomStream);


void SetSeed(const RandomStream& s)
{
   NTL_TLS_GLOBAL_ACCESS(CurrentRandomStream);

   if (!CurrentRandomStream)
      CurrentRandomStream.make(s);
   else
      *CurrentRandomStream = s;
}


void SetSeed(const unsigned char *data, long dlen)
{
   if (dlen < 0) LogicError("SetSeed: bad args");

   Vec<unsigned char> key;
   key.SetLength(NTL_PRG_KEYLEN);
   DeriveKey(key.elts(), NTL_PRG_KEYLEN, data, dlen);
 
   SetSeed(RandomStream(key.elts()));
}

void SetSeed(const ZZ& seed)
{
   long nb = NumBytes(seed);

   Vec<unsigned char> buf;
   buf.SetLength(nb);

   BytesFromZZ(buf.elts(), seed, nb);

   SetSeed(buf.elts(), nb);
}


static
void InitRandomStream()
{
   const std::string& id = UniqueID();
   SetSeed((const unsigned char *) id.c_str(), id.length());
}

static inline
RandomStream& LocalGetCurrentRandomStream()
{
   NTL_TLS_GLOBAL_ACCESS(CurrentRandomStream);

   if (!CurrentRandomStream) InitRandomStream();
   return *CurrentRandomStream;
}

RandomStream& GetCurrentRandomStream()
{
   return LocalGetCurrentRandomStream();
}







static inline
unsigned long WordFromBytes(const unsigned char *buf, long n)
{
   unsigned long res = 0;
   long i;

   for (i = n-1; i >= 0; i--)
      res = (res << 8) | buf[i];

   return res;
}


unsigned long RandomWord()
{
   RandomStream& stream = LocalGetCurrentRandomStream();
   unsigned char buf[NTL_BITS_PER_LONG/8];

   stream.get(buf, NTL_BITS_PER_LONG/8);
   return WordFromBytes(buf, NTL_BITS_PER_LONG/8);
}


void VectorRandomWord(long k, unsigned long* x)
{
   RandomStream& stream = LocalGetCurrentRandomStream();
   unsigned char buf[NTL_BITS_PER_LONG/8];

   for (long i = 0; i < k; i++) {
      stream.get(buf, NTL_BITS_PER_LONG/8);
      x[i] = WordFromBytes(buf, NTL_BITS_PER_LONG/8);
   }
}

long RandomBits_long(long l)
{
   if (l <= 0) return 0;
   if (l >= NTL_BITS_PER_LONG) 
      ResourceError("RandomBits: length too big");

   RandomStream& stream = LocalGetCurrentRandomStream();
   unsigned char buf[NTL_BITS_PER_LONG/8];
   long nb = (l+7)/8;
   stream.get(buf, nb);

   return long(WordFromBytes(buf, nb) & ((1UL << l)-1UL)); 
}

unsigned long RandomBits_ulong(long l)
{
   if (l <= 0) return 0;
   if (l > NTL_BITS_PER_LONG) 
      ResourceError("RandomBits: length too big");

   RandomStream& stream = LocalGetCurrentRandomStream();
   unsigned char buf[NTL_BITS_PER_LONG/8];
   long nb = (l+7)/8;
   stream.get(buf, nb);
   unsigned long res = WordFromBytes(buf, nb);
   if (l < NTL_BITS_PER_LONG)
      res = res & ((1UL << l)-1UL);
   return res;
}

long RandomLen_long(long l)
{
   if (l <= 0) return 0;
   if (l == 1) return 1;
   if (l >= NTL_BITS_PER_LONG) 
      ResourceError("RandomLen: length too big");

   RandomStream& stream = LocalGetCurrentRandomStream();
   unsigned char buf[NTL_BITS_PER_LONG/8];
   long nb = ((l-1)+7)/8;
   stream.get(buf, nb);
   unsigned long res = WordFromBytes(buf, nb);
   unsigned long mask = (1UL << (l-1)) - 1UL;
   return long((res & mask) | (mask+1UL)); 
}


long RandomBnd(long bnd)
{
   if (bnd <= 1) return 0;

   RandomStream& stream = LocalGetCurrentRandomStream();
   unsigned char buf[NTL_BITS_PER_LONG/8];
   long l = NumBits(bnd-1);
   long nb = (l+7)/8;

   long tmp;
   do {
      stream.get(buf, nb);
      tmp = long(WordFromBytes(buf, nb) & ((1UL << l)-1UL));
   } while (tmp >= bnd);

   return tmp;
}



void RandomBits(ZZ& x, long l)
{
   if (l <= 0) {
      x = 0;
      return;
   }

   if (NTL_OVERFLOW(l, 1, 0))
      ResourceError("RandomBits: length too big");

   RandomStream& stream = LocalGetCurrentRandomStream();

   long nb = (l+7)/8;
   unsigned long mask = (1UL << (8 - nb*8 + l)) - 1UL;

   NTL_TLS_LOCAL(Vec<unsigned char>, buf_mem);
   Vec<unsigned char>::Watcher watch_buf_mem(buf_mem);

   buf_mem.SetLength(nb);
   unsigned char *buf = buf_mem.elts();

   x.SetSize((l + NTL_ZZ_NBITS - 1)/NTL_ZZ_NBITS);
   // pre-allocate to ensure strong ES

   stream.get(buf, nb);
   buf[nb-1] &= mask;
   
   ZZFromBytes(x, buf, nb);
}


void RandomLen(ZZ& x, long l)
{
   if (l <= 0) {
      x = 0;
      return;
   }

   if (l == 1) {
      x = 1;
      return;
   }

   if (NTL_OVERFLOW(l, 1, 0))
      ResourceError("RandomLen: length too big");

   RandomStream& stream = LocalGetCurrentRandomStream();

   long nb = (l+7)/8;
   unsigned long mask = (1UL << (8 - nb*8 + l)) - 1UL;

   NTL_TLS_LOCAL(Vec<unsigned char>, buf_mem);
   Vec<unsigned char>::Watcher watch_buf_mem(buf_mem);

   buf_mem.SetLength(nb);
   unsigned char *buf = buf_mem.elts();

   x.SetSize((l + NTL_ZZ_NBITS - 1)/NTL_ZZ_NBITS);
   // pre-allocate to ensure strong ES

   stream.get(buf, nb);
   buf[nb-1] &= mask;
   buf[nb-1] |= ((mask >> 1) + 1UL);
   
   ZZFromBytes(x, buf, nb);
}





/**********************************************************

The following implementation of RandomBnd is designed
for speed.  It certainly is not resilient against a
timing side-channel attack (but then again, none of these
PRG routines are designed to be).

The naive strategy generates random candidates of the right 
bit length until the candidate < bnd.
The idea in this implementation is to generate the high
order two bytes of the candidate first, and compare this
to the high order two bytes of tmp.  We can discard the
candidate if this is already too large.

***********************************************************/

void RandomBnd(ZZ& x, const ZZ& bnd)
{
   if (bnd <= 1) {
      x = 0;
      return;
   }

   RandomStream& stream = LocalGetCurrentRandomStream();

   long l = NumBits(bnd);
   long nb = (l+7)/8;

   if (nb <= 3) {
      long lbnd = conv<long>(bnd);
      unsigned char lbuf[3];
      long ltmp;
      
      x.SetSize((l + NTL_ZZ_NBITS - 1)/NTL_ZZ_NBITS);
      // pre-allocate to ensure strong ES
      do {
         stream.get(lbuf, nb);
         ltmp = long(WordFromBytes(lbuf, nb) & ((1UL << l)-1UL));
      } while (ltmp >= lbnd);

     conv(x, ltmp);
     return;
   }

   // deal with possible alias
   NTL_ZZRegister(tmp_store);
   const ZZ& bnd_ref = ((&x == &bnd) ? (tmp_store = bnd) : bnd); 


   NTL_ZZRegister(hbnd);
   RightShift(hbnd, bnd_ref, (nb-2)*8);
   long lhbnd = conv<long>(hbnd);

   unsigned long mask = (1UL << (16 - nb*8 + l)) - 1UL;

   NTL_TLS_LOCAL(Vec<unsigned char>, buf_mem);
   Vec<unsigned char>::Watcher watch_buf_mem(buf_mem);
   buf_mem.SetLength(nb);
   unsigned char *buf = buf_mem.elts();

   unsigned char hbuf[2];

   x.SetSize((l + NTL_ZZ_NBITS - 1)/NTL_ZZ_NBITS);
   // pre-allocate to ensure strong ES
   for (;;) {
      stream.get(hbuf, 2);
      long hpart = long(WordFromBytes(hbuf, 2) & mask);

      if (hpart > lhbnd) continue;

      stream.get(buf, nb-2);
      buf[nb-2] = ((unsigned long) hpart);
      buf[nb-1] = ((unsigned long) hpart) >> 8; 

      ZZFromBytes(x, buf, nb);
      if (hpart < lhbnd || x < bnd_ref) break;
   }
}




// More prime generation stuff...

static
double Log2(double x)
{
   static const double log2 = log(2.0); // GLOBAL (relies on C++11 thread-safe init)
   return log(x)/log2;
}

// Define p(k,t) to be the conditional probability that a random, odd, k-bit 
// number is composite, given that it passes t iterations of the 
// Miller-Rabin test.
// If this routine returns a non-zero value, then
//    p(k,t) <= 2^{-n}.
// This basically encodes the estimates of Damgard, Landrock, and Pomerance;
// it uses floating point arithmetic, but is coded in such a way
// that its results should be correct, assuming that the log function
// is computed with reasonable precision.
// 
// It is assumed that k >= 3 and t >= 1; if this does not hold,
// then 0 is returned.

static
long ErrBoundTest(long kk, long tt, long nn)

{
   const double fudge = (1.0 + 1024.0/NTL_FDOUBLE_PRECISION);
   const double log2_3 = Log2(3.0);
   const double log2_7 = Log2(7.0);
   const double log2_20 = Log2(20.0);

   double k = kk;
   double t = tt;
   double n = nn;

   if (k < 3 || t < 1) return 0;
   if (n < 1) return 1;

   // the following test is largely academic
   if (9*t > NTL_FDOUBLE_PRECISION) LogicError("ErrBoundTest: t too big");

   double log2_k = Log2(k);

   if ((n + log2_k)*fudge <= 2*t)
      return 1;

   if ((2*log2_k + 4.0 + n)*fudge <= 2*sqrt(k))
      return 2;

   if ((t == 2 && k >= 88) || (3 <= t && 9*t <= k && k >= 21)) {
      if ((1.5*log2_k + t + 4.0 + n)*fudge <= 0.5*Log2(t) + 2*(sqrt(t*k)))
         return 3;
   }

   if (k <= 9*t && 4*t <= k && k >= 21) {
      if ( ((log2_3 + log2_7 + log2_k + n)*fudge <= log2_20 + 5*t)  &&
           ((log2_3 + (15.0/4.0)*log2_k + n)*fudge <= log2_7 + k/2 + 2*t) &&
           ((2*log2_3 + 2 + log2_k + n)*fudge <= k/4 + 3*t) )
         return 4; 
   }

   if (4*t >= k && k >= 21) {
      if (((15.0/4.0)*log2_k + n)*fudge <= log2_7 + k/2 + 2*t)
         return 5;
   }

   return 0;
}


void GenPrime(ZZ& n, long k, long err)
{
   if (k <= 1) LogicError("GenPrime: bad length");

   if (k > (1L << 20)) ResourceError("GenPrime: length too large");

   if (err < 1) err = 1;
   if (err > 512) err = 512;

   if (k == 2) {
      if (RandomBnd(2))
         n = 3;
      else
         n = 2;

      return;
   }


   long t;

   t = 1;
   while (!ErrBoundTest(k, t, err))
      t++;

   RandomPrime(n, k, t);
}


long GenPrime_long(long k, long err)
{
   if (k <= 1) LogicError("GenPrime: bad length");

   if (k >= NTL_BITS_PER_LONG) ResourceError("GenPrime: length too large");

   if (err < 1) err = 1;
   if (err > 512) err = 512;

   if (k == 2) {
      if (RandomBnd(2))
         return 3;
      else
         return 2;
   }

   long t;

   t = 1;
   while (!ErrBoundTest(k, t, err))
      t++;

   return RandomPrime_long(k, t);
}

void MultiThreadedGenGermainPrime(ZZ& n, long k, long err)
{
   long nt = AvailableThreads();


   long prime_bnd = ComputePrimeBound(k);

   if (NumBits(prime_bnd) >= k/2)
      prime_bnd = (1L << (k/2-1));

   ZZ two;
   two = 2;

   const long LOCAL_ITER_BOUND = 8;
   // since resetting the PRG comes at a certain cost,
   // we perform a few iterations with each reset to
   // amortize the reset cost. 

   unsigned long initial_counter = 0;
   ZZ seed;
   RandomBits(seed, 256);


   ZZ overflow_counter;

   for (;;) {

      AtomicLowWaterMark low_water_mark(-1UL);
      AtomicCounter counter(initial_counter);

      Vec< UniquePtr<ZZ> > result(INIT_SIZE, nt);
      Vec<unsigned long> result_ctr(INIT_SIZE, nt, -1UL);

      NTL_EXEC_INDEX(nt, index)

         RandomStreamPush push;

	 SetSeed(seed);
	 RandomStream& stream = GetCurrentRandomStream();

	 ZZ cand, n1;
         PrimeSeq s;

	 while (low_water_mark == -1UL) {

	    unsigned long local_ctr = counter.inc();
            if (local_ctr >> (NTL_BITS_PER_NONCE-1)) {
               // counter overflow...rather academic
               break;
            }

	    stream.set_nonce(local_ctr);
	    
	    for (long iter = 0; iter < LOCAL_ITER_BOUND && 
				local_ctr <= low_water_mark; iter++) {


	       RandomLen(cand, k);
	       if (!IsOdd(cand)) add(cand, cand, 1);

	       s.reset(3);
	       long p;

	       long sieve_passed = 1;

	       p = s.next();
	       while (p && p < prime_bnd) {
		  long r = rem(cand, p);

		  if (r == 0) {
		     sieve_passed = 0;
		     break;
		  }

		  // test if 2*r + 1 = 0 (mod p)
		  if (r == p-r-1) {
		     sieve_passed = 0;
		     break;
		  }

		  p = s.next();
	       }

               if (!sieve_passed) continue;


               if (MillerWitness(cand, two)) continue;

	       // n1 = 2*cand+1
	       mul(n1, cand, 2);
	       add(n1, n1, 1);


               if (MillerWitness(n1, two)) continue;

	       result[index].make(cand);
	       result_ctr[index] = local_ctr;
	       low_water_mark.UpdateMin(local_ctr);
	       break;
            }
         }

      NTL_EXEC_INDEX_END

      // find index of low_water_mark

      unsigned long low_water_mark1 = low_water_mark;
      long low_water_index = -1;

      for (long index = 0; index < nt; index++) {
	 if (result_ctr[index] == low_water_mark1) {
	    low_water_index = index;
	    break;
	 }
      }

      if (low_water_index == -1) {
         // counter overflow...rather academic
         overflow_counter++;
         initial_counter = 0;
         RandomBits(seed, 256);
         continue;
      }

      ZZ N;
      N = *result[low_water_index];

      ZZ iter = ((overflow_counter << (NTL_BITS_PER_NONCE-1)) +
                 conv<ZZ>(low_water_mark1) + 1)*LOCAL_ITER_BOUND;

      // now do t M-R iterations...just to make sure
 
      // First compute the appropriate number of M-R iterations, t
      // The following computes t such that 
      //       p(k,t)*8/k <= 2^{-err}/(5*iter^{1.25})
      // which suffices to get an overall error probability of 2^{-err}.
      // Note that this method has the advantage of not requiring 
      // any assumptions on the density of Germain primes.

      long err1 = max(1, err + 7 + (5*NumBits(iter) + 3)/4 - NumBits(k));
      long t;
      t = 1;
      while (!ErrBoundTest(k, t, err1))
         t++;

      Vec<ZZ> W(INIT_SIZE, t);

      for (long i = 0; i < t; i++) {
         do { 
            RandomBnd(W[i], N);
         } while (W[i] == 0);
      }

      AtomicBool tests_pass(true);

      NTL_EXEC_RANGE(t, first, last)

         for (long i = first; i < last && tests_pass; i++) {
            if (MillerWitness(N, W[i])) tests_pass = false;
         }

      NTL_EXEC_RANGE_END

      if (tests_pass) {
         n = N;
         return;
      }

      // very unlikey to get here
      initial_counter = low_water_mark1 + 1;
   }
}

void GenGermainPrime(ZZ& n, long k, long err)
{
   if (k <= 1) LogicError("GenGermainPrime: bad length");

   if (k > (1L << 20)) ResourceError("GenGermainPrime: length too large");

   if (err < 1) err = 1;
   if (err > 512) err = 512;

   if (k == 2) {
      if (RandomBnd(2))
         n = 3;
      else
         n = 2;

      return;
   }

   if (k >= 192) {
      MultiThreadedGenGermainPrime(n, k, err);
      return;
   }


   long prime_bnd = ComputePrimeBound(k);

   if (NumBits(prime_bnd) >= k/2)
      prime_bnd = (1L << (k/2-1));


   ZZ two;
   two = 2;

   ZZ n1;

   
   PrimeSeq s;

   ZZ iter;
   iter = 0;


   for (;;) {
      iter++;

      RandomLen(n, k);
      if (!IsOdd(n)) add(n, n, 1);

      s.reset(3);
      long p;

      long sieve_passed = 1;

      p = s.next();
      while (p && p < prime_bnd) {
         long r = rem(n, p);

         if (r == 0) {
            sieve_passed = 0;
            break;
         }

         // test if 2*r + 1 = 0 (mod p)
         if (r == p-r-1) {
            sieve_passed = 0;
            break;
         }

         p = s.next();
      }

      if (!sieve_passed) continue;


      if (MillerWitness(n, two)) continue;

      // n1 = 2*n+1
      mul(n1, n, 2);
      add(n1, n1, 1);


      if (MillerWitness(n1, two)) continue;

      // now do t M-R iterations...just to make sure
 
      // First compute the appropriate number of M-R iterations, t
      // The following computes t such that 
      //       p(k,t)*8/k <= 2^{-err}/(5*iter^{1.25})
      // which suffices to get an overall error probability of 2^{-err}.
      // Note that this method has the advantage of not requiring 
      // any assumptions on the density of Germain primes.

      long err1 = max(1, err + 7 + (5*NumBits(iter) + 3)/4 - NumBits(k));
      long t;
      t = 1;
      while (!ErrBoundTest(k, t, err1))
         t++;

      ZZ W;
      long MR_passed = 1;

      long i;
      for (i = 1; i <= t; i++) {
         do {
            RandomBnd(W, n);
         } while (W == 0);
         // W == 0 is not a useful candidate witness!

         if (MillerWitness(n, W)) {
            MR_passed = 0;
            break;
         }
      }

      if (MR_passed) break;
   }
}

void OldGenGermainPrime(ZZ& n, long k, long err)
{
   if (k <= 1) LogicError("GenGermainPrime: bad length");

   if (k > (1L << 20)) ResourceError("GenGermainPrime: length too large");

   if (err < 1) err = 1;
   if (err > 512) err = 512;

   if (k == 2) {
      if (RandomBnd(2))
         n = 3;
      else
         n = 2;

      return;
   }


   long prime_bnd = ComputePrimeBound(k);

   if (NumBits(prime_bnd) >= k/2)
      prime_bnd = (1L << (k/2-1));


   ZZ two;
   two = 2;

   ZZ n1;

   
   PrimeSeq s;

   ZZ iter;
   iter = 0;


   for (;;) {
      iter++;

      RandomLen(n, k);
      if (!IsOdd(n)) add(n, n, 1);

      s.reset(3);
      long p;

      long sieve_passed = 1;

      p = s.next();
      while (p && p < prime_bnd) {
         long r = rem(n, p);

         if (r == 0) {
            sieve_passed = 0;
            break;
         }

         // test if 2*r + 1 = 0 (mod p)
         if (r == p-r-1) {
            sieve_passed = 0;
            break;
         }

         p = s.next();
      }

      if (!sieve_passed) continue;


      if (MillerWitness(n, two)) continue;

      // n1 = 2*n+1
      mul(n1, n, 2);
      add(n1, n1, 1);


      if (MillerWitness(n1, two)) continue;

      // now do t M-R iterations...just to make sure
 
      // First compute the appropriate number of M-R iterations, t
      // The following computes t such that 
      //       p(k,t)*8/k <= 2^{-err}/(5*iter^{1.25})
      // which suffices to get an overall error probability of 2^{-err}.
      // Note that this method has the advantage of not requiring 
      // any assumptions on the density of Germain primes.

      long err1 = max(1, err + 7 + (5*NumBits(iter) + 3)/4 - NumBits(k));
      long t;
      t = 1;
      while (!ErrBoundTest(k, t, err1))
         t++;

      ZZ W;
      long MR_passed = 1;

      long i;
      for (i = 1; i <= t; i++) {
         do {
            RandomBnd(W, n);
         } while (W == 0);
         // W == 0 is not a useful candidate witness!

         if (MillerWitness(n, W)) {
            MR_passed = 0;
            break;
         }
      }

      if (MR_passed) break;
   }
}

long GenGermainPrime_long(long k, long err)
{
   if (k >= NTL_BITS_PER_LONG-1)
      ResourceError("GenGermainPrime_long: length too long");

   ZZ n;
   GenGermainPrime(n, k, err);
   return to_long(n);
}


NTL_END_IMPL
