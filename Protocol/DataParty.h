// (C) 2018 University of NKU. Free for used

/*
 * DataParty.h
 *
 */

#ifndef DATA_PARTY_H_
#define DATA_PARTY_H_


/*
 * @Params of the Class DP
 *		@n_computation_machines 			number of the computation parties
 *		@n_data_machines 					    number of the data parties
 *    @PortnumBase  						    Port number base to attempt to start connections from 
 *		@mynum								        the NO. mark of the current party(numbered in data party)
 *		@my_port							        the port that the server sokcet listening
 *		@hostname							        Host where CP.x is running to collecting data (default: localhost). Ignored if --ip-file-name is used.
 *		@ipFileName         				  Filename containing list of party ip addresses. Alternative to --hostname and running CP.x  data collection.
 *    @lgp                          length of the prime number
 *    @oM                           Number of FlajoletMartin trails 
 *    @uN                           number of unique items
 *
 */

#include <vector>
#include <string>
#include <pthread.h>
#include "Tools/sha1.h"
#include "Tools/time-func.h"
#include "Tools/octetStream.h"
#include "Math/gfp.h"


class DP
{
private:
	int n_computation_machines;
	int n_data_machines;
  int PortnumBase;
  int mynum;
  int my_port;
  int lgp; 
  std::vector<int> socket_num;
  std::string hostname;
  std::string ipFileName;
  std::string PREP_DATA_PREFIX;

  int oM;
  int uN;
  int ofs_w;
  int nthreads;

  // for safe socket communication
  bool connection_active;
  mutable blk_SHA_CTX ctx;
  mutable size_t sent;
  mutable Timer timer;

  // for multithread
  pthread_mutex_t mutex_go;
  int Go;

public:

	static const int DEFAULT_PORT = -1;
	static const int OFFSET_PORT = 2000;

  DP(int argc, const char** argv);
  void start();

  /*
  * @Phase: Data Collection Phase-- In this phase data party emulate sending mM trails of oblivious flajolet martin structures with N computation parties
  * @Return of the phase
  *     @bool           the phase successfully is running or not
  */
  bool data_collection();

  // this auxilary protocol help computation parties to share zero-array with length of M for the protocol of the Extract Z
  void share_zero_array();

  /*
  * @Protocol: data party run this protocol to distribute the share of a plaintext x (in gfp field) to all parties
  * @Params of the protocol share_value
  *     @x              the data party to share value (we use in typo: gfp)
  * @Return of the protocol share_value
  *     @bool           the protocol successfully is running or not
  */
  template <class T>
  bool share_value(T& x);

  template <class T>
  bool bucket_share_value(std::vector<std::vector<T>>& des_vec);

  template <class T>
  bool multi_thread_share_value(T& x);

  template <class T>
  bool multi_thread_bucket_share_value(std::vector<std::vector<T>>& des_vec);


  /*
  * for multi-thread purpose communication functionality
  */
  static void* thread_receive_player(void* arg);
  static void* thread_bucket_receive_player(void* arg);
  static void* thread_send_to(void* arg);
  static void* thread_bucket_send_to(void* arg);

  template <class T>
  //void run_thread_receive_player(int& player_no, T& data, octetStream& o,const int thread_id);
  void run_thread_receive_player(int& player_no, T& data, octetStream& o);

  template <class T>
  void run_thread_bucket_receive_player(int& player_no, std::vector<T>& data, octetStream& o);

  template <class T>
  //void run_thread_send_to(int& player_no, T& data, octetStream& o,const int thread_id);
  void run_thread_send_to(int& player_no, T& data, octetStream& o);

  void run_thread_bucket_send_to(int& player_no, octetStream& o);


  /*
  * basic communication functionality
  */
  // receive an octetstream from computation player i
  void receive_player(int i,octetStream& o,bool donthash=false) const;

  // Send an octetStream to computation player i
  void send_to(int i,const octetStream& o,bool donthash=false) const;
};



#endif /* DATA_PARTY_H_ */

