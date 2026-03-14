#include <iostream>
#include <cstdlib>

/* output (compiler_name,language_standard,cpu_type)

   compiler_name:
   We recognize "gcc", "clang", and "icc".
   Other compilers are named "unknown".

   language_standard:
   one of 0, 1997, 2011, 2014, 2017, ...

   cpu_type:
   We recognize all x86 machines as "x86".
   We recognize all arm machines as "arm".
   By special request, we recognize IBM Z as "s390x"
   Other CPUs are named "unknown".

   os_name:
   We recognize "linux".
   Other operating systems are named "unknown".


*/

#ifndef __cplusplus
#define __cplusplus 0
#endif

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)


using namespace std;

int main()
{
   long language_standard = atoi(TOSTRING(__cplusplus));
   language_standard /= 100;
   if (language_standard < 1997) language_standard = 0;

   const char *compiler_name = "unknown";
   const char *cpu_type = "unknown";
   const char *os_name = "unknown";

#if defined(__INTEL_COMPILER)
   compiler_name = "icc";
#elif defined(__clang__)
   compiler_name = "clang";
#elif defined (__GNUC__)
   compiler_name = "gcc";
#endif

#if defined(__x86_64__) || defined(__x86_64) || defined(__i386__) || defined(__i386)
   cpu_type = "x86";
#elif defined(__ARM_ARCH) || defined(_M_ARM) || defined(_M_ARM64) || defined(_M_ARM64EC) || defined(__arm__) || defined(__thumb__) || defined(__TARGET_ARCH_ARM) || defined(_ARM) || defined(__aarch64__)  
   cpu_type="arm";
#elif defined(__s390x__)
   cpu_type = "s390x";
#endif

#if defined(__linux__)
   os_name = "linux";
#endif

   cout << "(" << compiler_name << "," << language_standard 
             << "," << cpu_type << "," << os_name << ")\n";

}
