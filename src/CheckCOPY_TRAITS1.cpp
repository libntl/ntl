#include <type_traits>

template<class T>
constexpr bool Relocate_aux_has_trivial_copy(T*)
{
   return  std::is_trivially_copyable<T>::value &&
           std::is_trivially_destructible<T>::value &&
           std::is_copy_constructible<T>::value;
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
         !Relocate_aux_has_trivial_copy((EE*)0) &&
          Relocate_aux_has_trivial_copy((X*)0) &&
          Relocate_aux_has_trivial_copy((XX*)0) &&
          Relocate_aux_has_trivial_copy((int*)0)
    )

      return 0;
   else
      return -1;
}

