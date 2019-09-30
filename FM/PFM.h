// (C) 2018 University of NKU. Free for used
// Author: stoneboat@mail.nankai.edu.cn

/*
 * PFM.h
 *
 */

#ifndef PFM_H_
#define PFM_H_

#include <string>
#include <vector>

#include <stdint.h>
#include <math.h>

// TODO: What was the point of using boost::scoped_array?
#include <openssl/rand.h>
#include <openssl/err.h>

/*
*   @params N the number of the unique item, we test 20k, 1 million, 1 billion
*   @params w the length of the flajoletmartin sketch respectly with N
*   @params M the trails of the flajoletmartin approximation
*/

namespace FlajoletMartin{
  constexpr uint32_t  N[3] = {20000,1000000,1000000000};
  constexpr int       w[3] = {19,24,34};
  constexpr int       M[3] = {500,1000,2000};

  class randomSketchGenerator{
  private:
    bool initialized;
    unsigned char key[32];
  public:
    randomSketchGenerator() {
      RAND_poll();
      initialized = true;
    }

    void reseed() {
      int written = RAND_bytes(key, sizeof(key));
      RAND_seed(key, written);
    }

    bool random_bytes(const size_t &byte_count, unsigned char* out)
    {
      int rc = RAND_bytes(out, byte_count);
      unsigned long err = ERR_get_error();

      if(rc != 0 && rc != 1) {
        fprintf(stderr, "Generating Random Bytes failed, err = 0x%lx\n", err);
        return false;
      }

      return true;
    }




  };
}




/*
 * @Params of the Class PFM  (bit-structed plaintext FlajoletMartin sketch)
 *    @w          the length of the flajoletmartin sketch respectly with N
 *    @M          the trails of the flajoletmartin approximation
 *    @median     this capture the case when M or N is too large determined by threshold
 *
 */

// =============================================================================

class BitsFlajoletMartin {
private:
  int M;
  int w;
  int N;
  uint64_t mask; // mask an uint64_t type to w-bits wise

  // for multi-thread
  int n_threads;

public:
  static double R;
  static double PHI;
  static int Go;
  static pthread_mutex_t mutex_go;
  /*
  * @Func       
  *            return the index of the least significant 1 bit in n (= count trailing zeros)
  * @Params
  *            @n  hash of the item
  */
  static int ls1(uint64_t n);

  /*
  * @Func       
  *            return the index of the least significant 0 bit in n (= count trailing ones)
  * @Params
  *            @n  hash of the item
  */
  static int ls0(uint64_t n);

  /*
  * @Func       
  *            return the 2 power of the index of the least significant 1 bit in n (= count trailing zeros)
                      i.e. in bits mode
  * @Params
  *            @n  hash of the item
  */
  static uint64_t vls1(uint64_t n);

  /*
  * @Func       
  *            return the 2 power of the index of the least significant 0 bit in n (= count trailing ones)
                      i.e. in bits mode
  * @Params
  *            @n  hash of the item
  */
  static uint64_t vls0(uint64_t n);

  BitsFlajoletMartin(const int _M, const int _N, int _n_threads = 4);
  int get_M() { return M; }
  int get_N() { return N; }
  int get_w() { return w; }

  void reset_N(const int _N) {
    N = _N;
    w = ceil(log2 (N))+4;

    mask = 0;
    int k =w;
    while(k--){
      mask = (mask << 1);
      mask |= 1;  
    }
  }
  void set_M(const int _M) {M = _M;}
  void set_N(const int _N) {
    reset_N(_N);
  }

  /*
  * for benchmark used 
  */
  double randomFM_approximateCount();
  double single_thread_randomFM_approximateCount();
  static void* multi_thread_randomFM_approximateCount(void* arg);

  /*
  * @Func generate the fake FlajoletMartin sketch structure
  * @Input Params 
  *        @count_num    the number of the unique items
  *        @fm_sketch    the placeholder store the output random fake sketch
  */
  void generate_fake_FM(const int& count_num, uint64_t& fm_sketch);
  void generate_fake_FM(const int& count_num, const int n_fm, std::vector<uint64_t>& fm_sketch);
  void check_fake_FM();

};


#endif /* PFM_H_ */