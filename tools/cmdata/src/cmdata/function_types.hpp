#ifndef _CMDATA_FUNCTION_TYPES_HPP
#define _CMDATA_FUNCTION_TYPES_HPP

namespace cmdata::ftypes
{
  template<typename T>
  struct function_traits;
  
  // Specialization for function pointers
  template<typename Ret, typename... Args>
  struct function_traits<Ret(*)(Args...)> {
      using return_type = Ret;
      using args_tuple = std::tuple<Args...>;
      // Define a type for the function signature
      using signature = Ret(Args...);
  };
  
  // Define a template function `do_nothing` that returns a no-op lambda
  template<typename function_traits>
  auto do_nothing() {
      return [](auto&&... args) -> typename function_traits::return_type {};
  }
}

#endif // _CMDATA_FUNCTION_TYPES_HPP