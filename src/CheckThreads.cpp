#include <NTL/config.h>
#include <NTL/new.h>
#include <atomic>
#include <thread>
#include <utility>
#include <cstdlib>

#include <iostream>

void TerminalError(const char *s)
{
   std::cerr << s << "\n";
   std::abort();
}

void MemoryError() { TerminalError("out of memory"); }
void ResourceError(const char *msg) { TerminalError(msg); }


#if (defined(NTL_THREADS) && defined(NTL_TLS_HACK)) 
#include <pthread.h>
#endif

#define NTL_THREAD_LOCAL thread_local

#ifdef __GNUC__
#define NTL_CHEAP_THREAD_LOCAL __thread
#else
#define NTL_CHEAP_THREAD_LOCAL thread_local
#endif


#define NTL_DETAILS_PTHREAD details_pthread


#if (defined(NTL_THREADS) && defined(NTL_TLS_HACK)) 

namespace details_pthread {


struct Node {
   Node *next;

   Node() : next(0) { }
   virtual ~Node() { }
};

template<class T>
struct DerivedNode : Node {
   T t;

   template<class... Args>
   DerivedNode(Args&&... args) : t(std::forward<Args>(args)...) { }
};

inline void 
delete_node(Node *p) noexcept { delete p;  }
// an exception here would likely lead to a complete mess...
// the noexcept specification should force an immediate termination

inline void 
delete_list(void *vp)
{
   Node *p = (Node *) vp;
   while (p) {
      Node *tmp = p;
      p = p->next;
      delete_node(tmp);
   }
}


using namespace std;
// I'm not sure if pthread stuff might be placed in namespace std

struct key_wrapper {
   pthread_key_t key;

   key_wrapper(void (*destructor)(void*))
   {
      if (pthread_key_create(&key, destructor))
         ResourceError("pthread_key_create failed");
   }
};


inline void
push_node(Node *p)
// this pushes a new node to the front to the list
// of objects that need to be deleted
{
   if (!p) MemoryError();

   static key_wrapper wkey(delete_list);
   // This relies on C++11 thread-safe static initialization.
   // It also relies on the guarantee that there is just one
   // global key (this requirement is only needed to 
   // limit the number of keys, not for correctness).

   p->next = (Node *) pthread_getspecific(wkey.key);

   if (pthread_setspecific(wkey.key, p)) {
      delete_node(p);
      ResourceError("pthread_setspecific failed");
   }
}

}


#define NTL_TLS_LOCAL_INIT(type, var, init)  \
   static NTL_CHEAP_THREAD_LOCAL NTL_DETAILS_PTHREAD::DerivedNode<type> *_ntl_hidden_variable_tls_local_ptr_ ## var = 0;  \
   NTL_DETAILS_PTHREAD::DerivedNode<type> *_ntl_hidden_variable_tls_local_ptr1_ ## var = _ntl_hidden_variable_tls_local_ptr_ ## var;  \
   if (!_ntl_hidden_variable_tls_local_ptr1_ ## var) {  \
      NTL_DETAILS_PTHREAD::DerivedNode<type> *_ntl_hidden_variable_tls_local_ptr2_ ## var = NTL_NEW_OP NTL_DETAILS_PTHREAD::DerivedNode<type> init;  \
      NTL_DETAILS_PTHREAD::push_node(_ntl_hidden_variable_tls_local_ptr2_ ## var); \
      _ntl_hidden_variable_tls_local_ptr1_ ## var = _ntl_hidden_variable_tls_local_ptr2_ ## var;  \
      _ntl_hidden_variable_tls_local_ptr_ ## var = _ntl_hidden_variable_tls_local_ptr1_ ## var;  \
   }  \
   type &var = _ntl_hidden_variable_tls_local_ptr1_ ## var->t  \



#else


// NOTE: this definition of NTL_TLS_LOCAL_INIT ensures that var names
// a local reference, regardless of the implementation
#define NTL_TLS_LOCAL_INIT(type,var,init) \
    static NTL_THREAD_LOCAL type _ntl_hidden_variable_tls_local ## var init; \
    type &var = _ntl_hidden_variable_tls_local ## var




#endif

#define NTL_EMPTY_ARG
#define NTL_TLS_LOCAL(type,var) NTL_TLS_LOCAL_INIT(type,var,NTL_EMPTY_ARG)

#define NTL_TLS_GLOBAL_DECL_INIT(type,var,init)  \
   typedef type _ntl_hidden_typedef_tls_access_ ## var;  \
   static inline  \
   type& _ntl_hidden_function_tls_access_ ## var() {  \
      NTL_TLS_LOCAL_INIT(type,var,init);  \
      return var;  \
   }  \


#define NTL_TLS_GLOBAL_DECL(type,var) NTL_TLS_GLOBAL_DECL_INIT(type,var,NTL_EMPTY_ARG)

#define NTL_TLS_GLOBAL_ACCESS(var) \
_ntl_hidden_typedef_tls_access_ ## var & var = _ntl_hidden_function_tls_access_ ## var()


//=======================

std::atomic_long count_con(0);
std::atomic_long count_des(0);
std::atomic_long count1(0);

struct X {
   long d;

   X() { d = count1++; count_con++; }
   ~X() { count_des++; }
};

NTL_TLS_GLOBAL_DECL(X,x)

void task(long *v)
{
   NTL_TLS_GLOBAL_ACCESS(x);
   *v = x.d;
}


#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

int main()
{
   long v1, v2, v3;
   std::thread t1(task, &v1);
   std::thread t2(task, &v2);
   std::thread t3(task, &v3);

   t1.join();
   t2.join();
   t3.join();

   //std::cout << count_con << "\n";
   //std::cout << count_des << "\n";
   //std::cout << v1 << " " << v2 << " " << v3 << "\n";

   long s1, s2, s3;
   s1 = MIN(MIN(v1,v2),v3);
   s3 = MAX(MAX(v1,v2),v3);
   s2 = v1+v2+v3-s1-s3;

   if (count_con != 3 || count_des != 3 || s1 != 0 || s2 != 1 || s3 != 2) {
      return -1;
   }
   return 0;
}
