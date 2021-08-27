#ifndef LINUX_S390X_H
#define LINUX_S390X_H

#if defined(__s390x__) && defined(__linux__) \
    && (defined(__GNUC__) || defined(__clang__))

#define LINUX_S390X

#include <sys/auxv.h>

/* message-security-assist extension 8 */
#define MSA8                     146
/* Map a facility bit number or function code to its bit mask. */
#define MASK64(n)                  \
           (1ULL << (63 - (n) % 64))
/* Map a facility bit number or function code to its offset. */
#define OFF64(n)            (n / 64)

/* Function codes */
#define CPACF_KMA_QUERY            0
#define CPACF_KMA_GCM_AES_256     20

/* Function code flags */
#define CPACF_KMA_LAAD         0x200  /* Last-AAD */
#define CPACF_KMA_HS           0x400  /* Hash-subkey Supplied */

static inline unsigned long
stfle(unsigned long flist[], unsigned long nmemb)
{
   register unsigned long r0 __asm__("0") = (unsigned long)nmemb - 1;

   __asm__ volatile(
      ".insn s,%[opc]<<16,0(%[flist])"
      : "+d" (r0)
      : [flist] "a" (flist), [opc] "i" (0xb2b0)
      : "memory", "cc"
   );

   return r0 + 1;
}

/*  KMA (cipher message with authentication) */
static inline void
cpacf_kma(unsigned long fc, void *param, unsigned char *out, const unsigned char *aad,
    unsigned long aadlen, const unsigned char *in, unsigned long inlen)
{
   register unsigned long r0 __asm__("0") = (unsigned long)fc;
   register unsigned long r1 __asm__("1") = (unsigned long)param;
   register unsigned long r2 __asm__("2") = (unsigned long)in;
   register unsigned long r3 __asm__("3") = (unsigned long)inlen;
   register unsigned long r4 __asm__("4") = (unsigned long)aad;
   register unsigned long r5 __asm__("5") = (unsigned long)aadlen;
   register unsigned long r6 __asm__("6") = (unsigned long)out;

   __asm__ volatile(
      "0:     .insn   rrf,%[opc]<<16,%[out],%[in],%[aad],0\n"
      "       brc     1,0b\n"     /* partial completion */
      : [out] "+a" (r6),
        [in] "+a" (r2), [inlen] "+d" (r3),
        [aad] "+a" (r4), [aadlen] "+d" (r5)
      : [fc] "d" (r0), [param] "a" (r1), [opc] "i" (0xb929)
      : "cc", "memory"
   );
}

#endif
#endif
