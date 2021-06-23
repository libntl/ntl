#include <NTL/ctools.h>
#include <NTL/linux_s390x.h>

#include <cstdlib>
#include <iostream>


using namespace std;

#if !defined(LINUX_S390X)
#error "KMA not supported"
#endif

int main()
{
#if defined(AT_HWCAP) && defined(HWCAP_S390_STFLE)
   unsigned long hwcap, facility_list_nmemb;
   uint64_t status_word[2], facility_list[3];

   /* Check for STFLE. */
   hwcap = getauxval(AT_HWCAP);
   if (!(hwcap & HWCAP_S390_STFLE))
      return -1;

   /* Query facility list. */
   facility_list_nmemb = stfle(facility_list, 3);

   /* Check MSA8. */
   if (facility_list_nmemb >= OFF64(MSA8) + 1
       && (facility_list[OFF64(MSA8)] & MASK64(MSA8))) {
         cpacf_kma(CPACF_KMA_QUERY, &status_word, NULL, NULL, 0, NULL, 0);

         if (status_word[OFF64(CPACF_KMA_GCM_AES_256)]
             & MASK64(CPACF_KMA_GCM_AES_256)) {
            return 0;
         }
   }
#endif
   return -1;
}
