#ifndef TIME_H_
#define TIME_H_

#ifdef USE_MPI
#include <mpi.h>
#define TIMERINIT(label)   \
  double start##label = 0; \
  double total##label = 0;

#define TIMERSTART(label) start##label = MPI_Wtime();

#define TIMERSTOP(label) total##label += MPI_Wtime() - start##label;

#define TIMERPRINT(label)                                                   \
  std::cout << "# elapsed time (" << #label << "): " << total##label << "s" \
            << std::endl;

#else
#include <chrono>
#define TIMERSTART(label)                                                \
  std::chrono::time_point<std::chrono::system_clock> a##label, b##label; \
  a##label = std::chrono::system_clock::now();

#define TIMERSTOP(label)                                                     \
  b##label = std::chrono::system_clock::now();                               \
  std::chrono::duration<double> delta##label = b##label - a##label;          \
  std::cout << "# elapsed time (" << #label << "): " << delta##label.count() \
            << "s" << std::endl;
#endif

#endif
