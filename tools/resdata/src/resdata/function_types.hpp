#ifndef _RESDATA_FUNCTION_TYPES_HPP
#define _RESDATA_FUNCTION_TYPES_HPP

namespace resdata::ftypes
{
  template<typename T>
  struct function_traits;
  
  template<typename Ret, typename... Args>
  struct function_traits<Ret(*)(Args...)> {
      using return_type = Ret;
      using args_tuple = std::tuple<Args...>;
      // Define a type for the function signature
      using signature = Ret(Args...);
  };
  
  // does nothing while taking the same arguments as the function
  template<typename function_traits>
  auto do_nothing() {
      return [](auto&&... args) -> typename function_traits::return_type {};
  }
}

#endif // _RESDATA_FUNCTION_TYPES_HPP