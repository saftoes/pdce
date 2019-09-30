// (C) 2018 University of NKU. Free for used
// Author: stoneboat@mail.nankai.edu.cn

/*
 * PFM.cpp
 *
 */

#include "PFM.h"

#include <farmhash.h>
#include <stdexcept>
#include <iostream>
#include <pthread.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>



// =============================================================================
double BitsFlajoletMartin::R = 0;
int BitsFlajoletMartin::Go = 0;
pthread_mutex_t BitsFlajoletMartin::mutex_go;
double BitsFlajoletMartin::PHI = 0.77351;

// Universal hashing: the point is to have a different hash function for each seed m
uint64_t hash(std::string s, int m) {
  return util::Hash64WithSeed(s.c_str(), s.length(), m);
}

uint64_t hash(char* s, int len, int m) {
  return util::Hash64WithSeed(s, len, m);
}


int BitsFlajoletMartin::ls1(uint64_t n) {
  // note __builtin_ctz() takes unsigned int, not long, so use __builtin_ctzl()
  return (n == 0 ? 8*sizeof(n) : __builtin_ctzl(n));
}

int BitsFlajoletMartin::ls0(uint64_t n) {
  // note __builtin_ctz() takes unsigned int, not long, so use __builtin_ctzl()
  return (n == 0 ? 8*sizeof(n) : __builtin_ctzl(~n));
}

uint64_t BitsFlajoletMartin::vls1(uint64_t n) {
  return (n & (-n));
}

uint64_t BitsFlajoletMartin::vls0(uint64_t n) {
  return ((~n) & (n + 1));
}

BitsFlajoletMartin::BitsFlajoletMartin(const int _M, const int _N, int _n_threads){
    M = _M;
    N = _N;
    n_threads = _n_threads;
    w = ceil(log2 (N))+4;

    mask = 0;
    int k =w;
    while(k--){
      mask = (mask << 1);
      mask |= 1;  
    }
}

/*
* for multi threads
*/
typedef struct randomFM_struct {
   int thread_id;
   int m_start;
   int m_interval;
   int N;
   uint64_t mask;
   double R;
   int w;
} thread_info;


// multi thread version, almostly the same as the single thread
void* BitsFlajoletMartin::multi_thread_randomFM_approximateCount(void* arg){
  thread_info* data = (thread_info*) arg;
  const int thread_M = data->m_interval;
  const int m_index = data->m_start;
  const int thread_N = data->N;
  const int byte_len = sizeof(int);
  int threshold = 10000000;
  threshold = threshold>thread_N ? thread_N : threshold;
  //char* string_tmp = new char[byte_len];
  char* string_tmp = new char[byte_len*threshold];
  char* char_p = string_tmp;
  const uint64_t thread_mask = data->mask;
  const int t_id = data->thread_id;
  double work_percentage = 0;

  std::vector<uint64_t> trails_fm_sketch(thread_M);

  for (int m = 0; m < thread_M; m++) {
    trails_fm_sketch[m] = 0;
  }

  // initial the random FM
  uint64_t tmp = 0;

  // time controling
  struct timeval t1;
  struct timeval t2;
  //double cpu_time_used;
  gettimeofday(&t1, NULL);



  for(int n=0; n< thread_N; n+=threshold){

    // generate the unique item
    for (int i = n; i < n+threshold; i++){
      memcpy((void*)char_p, (void*)&i,byte_len);
      char_p+=byte_len;
    }
    //reset the pointer
    char_p = string_tmp;

    // generate the random sketch
    for (int i=0; i<threshold; i++){
      for (int m = 0; m < thread_M; m++){
        tmp = hash(char_p,byte_len,m+m_index) & thread_mask;
        trails_fm_sketch[m] |= vls1(tmp);
        
      }
      char_p+=byte_len;
    }
    //reset the pointer
    char_p = string_tmp;

    // report the work percentage
    if(t_id == 0 && ((double(n)/thread_N - work_percentage) > 0.01)){
        work_percentage = double(n)/thread_N;
        gettimeofday(&t2, NULL);
        printf("the work has complete %f using %ld (s) \n",work_percentage,(t2.tv_sec-t1.tv_sec));
    }
    
  }
  

  // calculate from the estimator
  double thread_R = 0;
  for (int m = 0; m < thread_M; m++) {
    thread_R += ls0(trails_fm_sketch[m]); //trailing ones
  }
  

  /*
  *  Report the end of the thread
  */
  pthread_mutex_lock(&mutex_go); 
  R+= thread_R;
  //std::cout<<" thread id\t" << data->thread_id << " R "<<thread_R <<"estimated count\t" <<pow(2.0, thread_R/thread_M)/PHI<< std::endl;
  Go++;
  pthread_mutex_unlock(&mutex_go);

  pthread_exit(NULL);
}

double BitsFlajoletMartin::single_thread_randomFM_approximateCount(){
  std::vector<uint64_t> trails_fm_sketch(M);

  //pseudorandom generator
  const int byte_len = ceil(sizeof(uint64_t)/sizeof(unsigned char));
  //const int valid_byte = 64*16;
  int threshold = 10000000;
  threshold = threshold>N ? N : threshold;
  unsigned char* string_tmp = new unsigned char[threshold*byte_len];
  unsigned char* unchar_p = string_tmp;

  FlajoletMartin::randomSketchGenerator rsg;

  // initial the random FM
  uint64_t tmp = 0;
  for (int m=0; m<M; m++){
    rsg.random_bytes(threshold*byte_len,string_tmp);
    unchar_p = string_tmp;
    for( int i=0; i<N; i++){
      memcpy((void*)&tmp, (void*)unchar_p,byte_len);
      trails_fm_sketch[m] |= vls1(tmp);
      unchar_p+= byte_len;
    }
  }


  // calculate from the estimator
  double R = 0;
  for (int m = 0; m < M; m++) {
    R += ls0(trails_fm_sketch[m]); //trailing ones
  }
  R /= M;
  return pow(2.0, R)/PHI;
}
double BitsFlajoletMartin::randomFM_approximateCount(){
  if(n_threads == 1){
    return single_thread_randomFM_approximateCount();
  }else{
    R = 0;
    Go = 0;
    mutex_go = PTHREAD_MUTEX_INITIALIZER;
    pthread_t* t = new pthread_t[n_threads];
    //double* n_R = new double[n_threads];
    thread_info* randomFM_data = new thread_info[n_threads];

    int m_interval = ceil(M/n_threads);
    int m_start = 0;
    for (int i=0;i <n_threads;i++){
      randomFM_data[i].m_start = m_start;
      randomFM_data[i].m_interval = m_interval;
      randomFM_data[i].N = N;
      randomFM_data[i].mask = mask;
      randomFM_data[i].w = w;
      randomFM_data[i].thread_id = i;
      // update the task id
      m_start+=m_interval;
    }

    for (int i=0;i <n_threads;i++){
      // setup the task
      pthread_create(&t[i], NULL, multi_thread_randomFM_approximateCount, (void*) &randomFM_data[i]);
    }

    while(Go < n_threads){
      usleep(100);
    }

    R /= M;

    double ret = pow(2.0, R)/PHI;


    return ret;
  }
}

void BitsFlajoletMartin::generate_fake_FM(const int& count_num, uint64_t& fm_sketch){
  FlajoletMartin::randomSketchGenerator rsg;

  int R = round(log2 (count_num*PHI));
  uint64_t trailing_mask = 0;
  fm_sketch = 0;

  unsigned char buffer[9];
  rsg.random_bytes(9,buffer);

  int rnd_bit = buffer[8] & 1;
  R-=rnd_bit;
  int k =R;
  while(k--){
    trailing_mask = (trailing_mask << 1);
    trailing_mask |= 1;  
  }

  memcpy((void*)&fm_sketch, (void*)buffer,8);

  fm_sketch |= trailing_mask;
}

void BitsFlajoletMartin::generate_fake_FM(const int& count_num, const int n_fm, std::vector<uint64_t>& fm_sketch){
  FlajoletMartin::randomSketchGenerator rsg;
  const int byte_len = sizeof(uint64_t)/sizeof(char);

  int R = floor(log2 (count_num*PHI));
  uint64_t trailing_mask = 0;
  uint64_t tmp = 0;

  unsigned char* buffer = new unsigned char[byte_len*n_fm+1];
  rsg.random_bytes(byte_len*n_fm+1,buffer);

  int rnd_bit = buffer[byte_len*n_fm] & 1;
  R+=rnd_bit;
  int k =R;
  while(k--){
    trailing_mask = (trailing_mask << 1);
    trailing_mask |= 1;  
  }

  int blk_size = 256;
  int threshold = (floor(n_fm/blk_size)) * blk_size;
  unsigned char* unchar_p = buffer;
  for(int i=0; i<threshold; i+=blk_size){
    for(int j=0; j<blk_size; j++){
      memcpy((void*)&tmp, (void*)unchar_p,byte_len);
      unchar_p+=byte_len;
      tmp |= trailing_mask;
      tmp &= mask;
      fm_sketch[i+j] = tmp;
    } 
  }

  for(int i=threshold; i<n_fm; i++){
    memcpy((void*)&tmp, (void*)unchar_p,byte_len);
    unchar_p+=byte_len;
    tmp |= trailing_mask;
    tmp &= mask;
    fm_sketch[i] = tmp;
  }
}

void BitsFlajoletMartin::check_fake_FM(){
  std::vector<uint64_t> fm_sketch(M);
  double R = 0;
  double ret = 0;

  generate_fake_FM(N,M,fm_sketch);

  for (int m = 0; m < M; m++) {
    R += ls0(fm_sketch[m]); //trailing ones
  }
  R /= M;
  ret = pow(2.0, R)/PHI;
  std::cout<<" number of unique items is "<<N<<"\t estimate size is "<< ret<<" with "<< M <<"\t trails\n";
}
  



// =============================================================================
