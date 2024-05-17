#include <hip/hip_runtime.h>

#include <iostream>

constexpr int error_exit_code = -1;

// Checks if the provided error code is hipSuccess and if not,
// prints an error message to the standard error output and terminates the
// program with an error code.
#define HIP_CHECK(condition)                                              \
  {                                                                       \
    const hipError_t error = condition;                                   \
    if (error != hipSuccess) {                                            \
      std::cerr << "An error encountered: \"" << hipGetErrorString(error) \
                << "\" at " << __FILE__ << ':' << __LINE__ << std::endl;  \
      std::exit(error_exit_code);                                         \
    }                                                                     \
  }