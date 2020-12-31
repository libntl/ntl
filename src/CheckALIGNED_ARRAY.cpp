

#define NTL_HAVE_ALIGNED_ARRAY
// DIRT: we need to define this here so that ctools.h
// tries the aligned array code...we only check
// that it compiles without error

#include <NTL/ctools.h>
#include <cstdlib>


using namespace std;


int main()
{
   return 0;
}



