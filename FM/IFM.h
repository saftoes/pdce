// (C) 2018 University of NKU. Free for used
// Author: stoneboat@mail.nankai.edu.cn

/*
 * OFM.h
 *
 */

#ifndef IFM_H_
#define IFM_H_

#include "PFM.h"

/*
 * @Params of the Class IFM  (integer-structed FlajoletMartin sketch)
 *    @fill_IFS          whether the IFS_group has been filled
 *
 */

// =============================================================================

class IntegerFlajoletMartin : public BitsFlajoletMartin {
private:
  int M;
  int w;
  int N;

  void init_IFS();

  bool fill_IFS;

public:
  std::vector<std::vector<int>> IFS_group;
  
  IntegerFlajoletMartin(const int _M, const int _N, int _n_threads = 4): BitsFlajoletMartin(_M, _N, _n_threads) { 
        M = get_M();
        N = get_N();
        w = get_w();
        init_IFS(); 
        fill_IFS = false;
      };

  void set_N(const int _N) {
    reset_N(_N);
    IFS_group.clear();
    fill_IFS = false;
    init_IFS();
  }

  void set_M(const int _M) {
    M = _M;
    IFS_group.clear();
    fill_IFS = false;
    init_IFS();
  }

  /*
  * @Func generate the fake Integer FlajoletMartin sketch structure
  * @Input Params 
  *        @count_num    the number of the unique items
  */
  void init_fake_IFM();
  void check_fake_IFM();

};


#endif /* OFM_H_ */