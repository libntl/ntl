#include <iostream>

/* output (compiler_name,language_standard,cpu_type)

   compiler_name:
   Right now, we just recognize "gcc", "clang", and "icc".
   Other compilers are named "unknown".

   language_standard:
   As of 2018, the available language standards are
   199711, 201103, 201402, 201703.

   cpu_type:
   Right now, we just recognize x86 and x86-64, and both are named "x86".
   Other CPUs are named "unknown".

*/

#ifndef __cplusplus
#define __cplusplus 1
#endif


int main()
{
   long language_standard = __cplusplus;
   // convert to one of 0, 1997, 2011, 2014, 2017
   if      (language_standard >= 201703) language_standard = 2017;
   else if (language_standard >= 201402) language_standard = 2014;
   else if (language_standard >= 201103) language_standard = 2011;
   else if (language_standard >= 199711) language_standard = 1997;
   else                                  language_standard = 0;

   const char *compiler_name = "unknown";
   const char *cpu_type = "unknown";
   const char *os_name = "unknown";

#if defined(__INTEL_COMPILER)
   compiler_name = "icc";
#elif defined(__clang__)
   compiler_name = "clang";
#elif defined (__GNUC__)
   compiler_name = "gcc";
#else
   compiler_name = "";
#endif

#if defined(__x86_64__) || defined(__x86_64) || defined(__i386__) || defined(__i386)
   cpu_type = "x86";
#elif defined(__s390x__)
   cpu_type = "s390x";
#endif

#if defined(__linux__)
   os_name = "linux";
#endif

   std::cout << "(" << compiler_name << "," << language_standard 
             << "," << cpu_type << "," << os_name << ")\n";

}
