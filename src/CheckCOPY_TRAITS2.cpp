#include <type_traits>

template <unsigned statement, typename out>
struct Relocate_aux_Failable
{
     typedef out type;
};

template <typename T>
struct Relocate_aux_has_copy {
     static const T *MakeT();

     template <typename U> // U and T are the same type
     static typename Relocate_aux_Failable<(unsigned(sizeof U(*MakeT()))), std::true_type>::type copy(int);

     template <typename U>
     static typename Relocate_aux_Failable<0U, std::false_type>::type copy(...);

     constexpr static bool value = decltype( copy<T>(0) )::value;
};

template<class T>
constexpr bool Relocate_aux_has_trivial_copy(T*)
{
   return  __has_trivial_copy(T) &&
           __has_trivial_destructor(T) &&
           Relocate_aux_has_copy<T>::value;
}

struct A { A(const A&) = delete; };
struct B { ~B() {} };
struct C { virtual void foo() { } };
struct D { private: D(const D&); };
struct E { private: E(const E&) = default; };
struct AA { A a; };
struct BB { B b; };
struct CC { C c; };
struct DD { D d; };
struct EE { E e; };

struct X { int x; };
struct XX { private: X x; };

int main()
{
   if (
         !Relocate_aux_has_trivial_copy((A*)0) &&
         !Relocate_aux_has_trivial_copy((B*)0) &&
         !Relocate_aux_has_trivial_copy((C*)0) &&
         !Relocate_aux_has_trivial_copy((D*)0) &&    
         !Relocate_aux_has_trivial_copy((E*)0) &&    
         !Relocate_aux_has_trivial_copy((AA*)0) &&
         !Relocate_aux_has_trivial_copy((BB*)0) &&
         !Relocate_aux_has_trivial_copy((CC*)0) &&
         !Relocate_aux_has_trivial_copy((DD*)0) &&
#ifndef __INTEL_COMPILER
         !Relocate_aux_has_trivial_copy((EE*)0) &&
// Ignore a bug in the Intel compiler.  It really does
// allow the use of this copy constructor (which should have been
// deleted), so the problem is not in the *detection* of this
// property.
#endif
          Relocate_aux_has_trivial_copy((X*)0) &&
          Relocate_aux_has_trivial_copy((XX*)0) &&
          Relocate_aux_has_trivial_copy((int*)0)
    )

      return 0;
   else
      return -1;
}

