//#include <argparse/argparse.hpp>
//#include <fstream>

//#include "octree.h"
//#include "sequential_bh_avx.h"
//#include "time_utils.h"
#include "sequential_ap.h"

int main(int argc, char** argv) {
  SequentialAPSimulate2(2048, 0.01f, 10.0f, 1);
  // -xAVX2 and -xMIC-AVX512 flags force the compiler to generate AVX2
  // and AVX-512 SIMD instructions, respectively. AVX2 extensions accelerated
  // the previous version by a factor of 7.4× while AVX-512 instructions
  // achieved a speedup of 15.1×
  // http://sedici.unlp.edu.ar/bitstream/handle/10915/95855/Documento_completo.pdf?sequence=1

  return 0;
}