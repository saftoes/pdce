// (C) 2018 University of NKU. Free for used

/*
 * ComputationParty.h
 *
 */

#ifndef COMPUTATION_PARTY_H_
#define COMPUTATION_PARTY_H_


/*
 * @Func Computation Party
 * @Params of the Class CP 
 * 		@lgp 				length of the prime number
 *		@nparties 			number of the parties
 *		@ndparties 			number of the data parties
 *      @pnbase  			Port number base to attempt to start connections from 
 *		@mynum			    the NO. mark of the current party
 *		@my_port			the port that the server sokcet listening
 *		@hostname			Host where Server.x is running to coordinate startup (default: localhost). Ignored if --ip-file-name is used.
 *		@PREP_DATA_PREFIX   The position where the offline data is stored 
 *		@playerNames		Names server, help to identify the other machine in the whole connected network
 *		@player  			The wrapped socket for the p2p communication
 *		@dataF       		read the share of the triple and random from the offline file
 *		@keyp       		the mac key in typo of gfp
 *		@alphai       		this party's share of the mac key in typo of gfp
 *    	@oM                 Number of FlajoletMartin trails 
 *    	@uN                 number of unique items
 */

#include <string>
#include <vector>
#include <array>
#include <fstream>
#include <pthread.h>

#include "Math/gfp.h"
#include "Math/Setup.h"
#include "Networking/Player.h"
#include "Processor/Data_Files.h"
#include "Networking/ServerSocket.h"
#include "Tools/time-func.h"
#include "Tools/sha1.h"

namespace PSUCA{
	constexpr int max_w2c = 300000;
	typedef struct CP_thread_wait2check_mac {
		gfp value;
		gfp mac;
	}W2C_mac;
}

class CP
{
private:
	int lgp; 
	int sec;
	int tau;
	int nparties; 
	int ndparties;
	int pnbase; 
	int mynum; 
	int my_port; 
	gfp keyp;
	gfp alphai;
	Data_Files* dataF;
	std::string hostname;
	std::string ipFileName;
	std::string PREP_DATA_PREFIX;

	/*
	*	Network Parameters
	*/
	bool local_envs;
	Player* player;
	std::vector<Player*> thread_player;
	vector<vector<int>> CP_sockets;
	mutable vector<Timer> commu_timer;

	int oM;
  	int uN;
  	int nthreads;

	std::vector<std::vector<Share<gfp>>> share_FS;
	std::vector<std::vector<Share<gfp>>> share_BFS;
	Share<gfp> share_Z;

	bool playerone;

	/*
	* use for the polynomial function
	*/
	vector<gfp> polyzero_f;
	std::vector<vector<gfp>> poly_cdb_f;
	vector<gfp> clookup_d2;
	vector<gfp> clookup_d3;


	/*
	*  Offline file handle
	*/
    std::vector<std::ifstream> rand_pool;
    std::vector<std::ifstream> rand2_pool;
    std::vector<std::ifstream> triple_pool;
    // default reading randExp maximum is 18
    std::vector<std::ifstream> randExp_pool;
    // for the zerotest reading, randExp maximum is 65
    std::vector<std::ifstream> randExp64_pool;


    // for multithread
  	pthread_mutex_t mutex_go;
  	int Go;
  	vector<pthread_mutex_t> mutex_local_go;
  	vector<int> local_Go;

  	static constexpr int port_increase = 75;
  	// for communcation thread: 0 for send and 1 for receive
  	static constexpr int communication_multiplier = 2;

  	// this is the d of the randExp refered protocol 7, the number is hard coded, acturally it is decided by the offline data file
  	int power_l_RExp;
  	int l_RExp;


public:
	void benchmark();

	size_t report_size();
 	/*
 	* @Protocol
 	* @Functionality: get the share of 0 or 1; if x is zero then get the share of 0, otherwise 1 Implementation of Protocol 2
 	* @Params of the protocol ZeroTest
 	* 		@[x] 					share of the value to be tested 
 	* @author Helger Lipmaa and Tomas 2013
 	*/
 	template <class T>
 	void lipmaa_zeroTest(const Share<T>& t, Share<T>& t0, int thread_num=0);

 	template <class T>
 	void bucket_lipmaa_zeroTest(const vector<vector<Share<T>>>& src_vec, std::vector<std::vector<Share<T>>>& t0, const  int m_start, int m_internal, int thread_num = 0, bool multi_thread =false);

 	template <class T>
 	void bucket_lipmaa_zeroTest_bench(const vector<vector<Share<T>>>& src_vec, std::vector<std::vector<Share<T>>>& t0, const  int m_start, int m_internal, int thread_num = 0, bool multi_thread =false);

 	static void* thread_lipmaa_zeroTest(void* arg);

 	template <class T>
 	void multi_thread_bucket_lipmaa_zeroTest(std::vector<std::vector<Share<T>>>& src_vec, std::vector<std::vector<Share<T>>>& t0);

 	template <class T>
 	void bucket_zeroTest(const std::vector<std::vector<Share<T>>>& src_vec, std::vector<std::vector<Share<T>>>& t0, bool multi_thread =false);


 	/*
 	* @Protocol
 	* @Functionality: for the merge share step just for benchmark
 	*/
 	void mergeShare();

 	/*
 	* @Protocol
 	* @Functionality: get the share of the summed estimator with size of M  Implementation of Protocol 7
 	* @Params of the protocol basic version ExtractZ
 	* 		@[fm_sketch] 			oblivious bfs with w length
 	* 		@z 						the estimator of the fm sketch
 	*/
 	template <class T>
 	void bExtractZ(std::vector<Share<T>> fm_sketch, Share<T>& z);

 	template <class T>
 	void bucket_bExtractZ(const std::vector<std::vector<Share<T>>>& src_fm, Share<T>& Z, const int m_start, int m_internal, int thread_num = 0, bool multi_thread =false);

 	template <class T>
 	void multi_thread_bucket_bExtractZ(std::vector<std::vector<Share<T>>>& src_fm, Share<T>& Z);

 	static void* thread_bExtractZ(void* arg);


 	/*
	*@Protocol
	*@Functionality: extract the estimator from the share-type vector, the protocol 8
	*					the basic version is just for clearer representing the functionality
	*					to achieve the higher performance for networking and batch operations we finally run the implementation multi_thread_bucket_CExtractZ
	*/
 	template <class T>
 	void CExtractZ_basic(const std::vector<Share<T>>& src_fm, Share<T>& z);

 	template <class T>
 	void CExtractZ(const std::vector<Share<T>>& src_fm, Share<T>& z, int thread_num = 0);

 	template <class T>
 	void bucket_CExtractZ_basic(const std::vector<std::vector<Share<T>>>& fm_sketch, Share<T>& Z, const int m_start, int m_internal, int thread_num = 0, bool multi_thread =false);

 	template <class T>
 	void bucket_CExtractZ_bench(const std::vector<std::vector<Share<T>>>& fm_sketch, Share<T>& Z, const int m_start, int m_internal, int thread_num = 0, bool multi_thread =false);

 	template <class T>
 	void multi_thread_bucket_CExtractZ(std::vector<std::vector<Share<T>>>& src_fm, Share<T>& Z);

 	static void* thread_CExtractZ(void* arg);

 	/*
	*@Auxiliary Protocol
	*@Functionality: evaluate the share-type array on polynomial function which return 1 if x are all 0 bits otherwise 0
	*/
 	template <class T>
 	void bdb_CExtractZ(const std::vector<Share<T>>& x, Share<T>& y, const int end, int thread_num = 0);

 	/*
	*@Auxiliary Protocol
	*@Functionality: evaluate the share-type array on polynomial function which return the index of first 1 bit
	*/
 	template <class T>
 	void lookup_x_CExtractZ(const std::vector<Share<T>>& x, Share<T>& y, const int start, const int d = 3, int thread_num = 0);


 	/*
 	* @Functionality:	get the offline data starting from R^-1, R^1, R^2, ..., R^len_group where len_group is at most l_RExp
 	*/
 	template <class T>
 	void getR_group(std::vector<Share<T>>& R_group, const int len_group, int thread_num = 0, int l_R = 19);

 	/*
 	* @Functionality:	get the offline data random bits each group has lgp+1 elements
 	*/
 	template <class T>
 	void getRB_group(std::vector<Share<T>>& rb_group, int thread_num = 0);

 	/*
 	* @Functionality:	get the share of the z which equals to xy where x,y in type of Share<T>
 	*					following implementation 
 										1)	for single xy=z
 										2)	for array-type X Y = Z
 										3)  for array-type X Y = Z where X = <mask, mask, ... ,mask>
 	*/
 	template <class T>
	void dot_product(Share<T>& x, Share<T>& y, Share<T>& z, int thread_num=0);

	template <class T>
	void dot_product(std::vector<Share<T>>& x, std::vector<Share<T>>& y, std::vector<Share<T>>& z, int thread_num=0);

	template <class T>
	void _dot_product(std::vector<Share<T>>& x, std::vector<Share<T>>& y, std::vector<Share<T>>& z, \
					 int n_vectors, int start_x = 0, int start_y = 0, int start_z = 0, int thread_num=0);

	template <class T>
	void _dot_product(Share<T>& mask, std::vector<Share<T>>& y, std::vector<Share<T>>& z, \
					 int n_vectors, int start_y = 0, int start_z = 0, int thread_num=0);

	template <class T>
	void _dot_product(vector<Share<T>>& mask, vector<std::vector<Share<T>>>& y, vector<std::vector<Share<T>>>& z, \
					 int n_vectors, int start_y = 0, int start_z = 0, int thread_num=0);

	/*
	*@Protocol
	*@Functionality: reveal the estimator Z of M trials and return the approximate count
	*/
	template <class T>
	double fm_approximateCount(const Share<T>& Z, int M);

	/*
	*@Functionality: startpoint and constructor, deconstructor
	*/
    CP(int argc, const char** argv);
    void start();
    ~CP();

private:
	/*
 	* wait-to-check mac : wait until up to threshold number of mac then to check it
 	*/
 	std::vector<int> n_wait2Check;
 	std::vector<std::vector<PSUCA::W2C_mac>> W2C_mac_queue;

	/*
 	* @Functionality: share-based mpc auxilary function
 	* 				  @open_and_check: directly get the part of mac from other parties, then check the correctness
 	*				  @batch_mac_check: check the correctness of the mac from the local baked data which get from other parties in advance
 	*				  @delay_open_and_check: get the part of mac from other parties, then put off the check phase in batch_mac_check or up to a threshold number
 	*/
	template <class T>
	void open_and_check(const Share<T>& share_value, T& a, T& mac);

	void batch_mac_check(int thread_num = -1);

	template <class T>
	void delay_open_and_check(const Share<T>& share_value, T& a, int thread_num = 0);

	template <class T>
	void delay_open_and_check(const std::vector<Share<T>>& share_value, std::vector<T>& a, int thread_num = 0);

	template <class T>
	void _delay_open_and_check(const std::vector<Share<T>>& share_value, std::vector<T>& a, int n_elements, int start_share = 0, int start_a = 0, int thread_num = 0);

	template <class T>
	void open_and_check(const std::vector<Share<T>>& share_value, std::vector<T>& a, std::vector<T>& mac);

	template <class T>
	void thread_open_and_check(const std::vector<Share<T>>& share_value, std::vector<T>& a, std::vector<T>& mac, int thread_num=0);

	/*
 	* @Functionality: network communication part; receive and send from other parties
 	*			Alghough SPDZ provides the same functionality, we re-implement it and use in the WAN environments
 	*				  This is because, the broadcast functionality of SPDZ cannot make full use of the bandwitdh of network, which
 	*			leads to the overhead in WAN much larger, thus we send and receive each socket parallelly in multi-thread mode
 	*				  To note that this optimizition only works in the case:
 	*																1) the total communication overhead is small, but communication times is big
 	*																2) the send and receiver buffer of tcp are big(you can tune it, commonly 6-12MB is enough)
 	*/
	void Broadcast_S_and_R(vector<octetStream>& o, int thread_num=0);
	void Broadcast_S_single(octetStream& o, int father_num, int player_no, bool multi_thread = false);
	void Broadcast_R_single(octetStream& o, int father_num, int player_no, bool multi_thread = false);

	static void* thread_Broadcast_single(void* arg);

	/*
 	* @Functionality: initialization auxiliary protocol for 
 	* 				  		1) opening offline data reader
 	*				  		2) close offline data reading stream
 	*				  		3) load the polynomial functions; if they didn't exist, generate first
 	*						4) initialize the mac check buffer 
 	*/
	void init_offline_data();
	void close_offline_data_stream();
	void init_Polynomial_coef(std::string folder_name, int tau);
	void init_wait_to_check_buffer();

	/*
	* @Functionality: get polynomial function
	*						1)	the auxiliary function for advanced version of extract-Z: CExtract
	*						2)	the auxiliary function for zerotest protocol
	*/
	const vector<gfp>& get_CExtractZ_PolyNomailFunc(int d);
	const vector<gfp>& get_zerotest_PolyNomailFunc();
};


/*
 * @Func The wrapped communication class with data party
 * @Params of the Class DP_connector
 * 		@default_port 		default receive port of the DP_connector
 *		@nplayers 			number of the data parties
 *      @portnum_base  		Port number base to attempt to start connections from 
 *		@mynum				the NO. mark of the current party(numbered in computation party)
 *		@my_port			the server listening port
 *		@alphai       		this party's share of the mac key in typo of gfp
 *		@keyp       		the mac key in typo of gfp
 *		@player  			The wrapped socket for the p2p communication
 *    	@oM                 Number of FlajoletMartin trails 
 *    	@uN                 number of unique items
 */


class DP_connector
{
private:
	int default_port(int playerno) { return portnum_base + playerno; }
	int nplayers;
  	int portnum_base;
  	int mynum;
  	int my_port;
  	gfp keyp;
  	gfp alphai;
  	Data_Files* dataF;
  	Player* player;

  	int oM;
  	int uN;
  	int nthreads;

  	void setup_server();

  	// for safe socket communication
  	mutable blk_SHA_CTX ctx;
  	mutable size_t sent;
  	mutable Timer timer;
public:
	size_t report_size();
	mutable ServerSocket* server;

	DP_connector(int _mynum,int pnb,int my_port, int _num_dp_players, int oM, int uN, int _nthreads = 1);
	DP_connector() : nplayers(-1), portnum_base(-1), mynum(-1), server(0){ 
		oM = 0;
		uN = 0;
		nthreads =1;
		dataF=0; 
		player = 0;
		sent = 0;
	}

	void init(int _mynum,int pnb,int my_port, int _num_dp_players);
	void key_init(gfp& _alphai, gfp& _gfp keyp);
	void player_init(Player* _player);
	void n_threads_init(int _nthreads = 1) {nthreads = _nthreads;}
	void params_ofs_init(int _oM, int _uN) {oM = _oM; uN = _uN;}


  	~DP_connector();

  	static const int DEFAULT_PORT = -1;
  	static const int OFFSET_PORT = 2000;
  	int num_players() const { return nplayers; }
  	int my_num() const { return mynum; }
  	int get_portnum_base() const { return portnum_base; } 

  	void start(Data_Files* _df, std::vector<std::vector<Share<gfp>>>& share_FS);

  	/*
	* @Phase: Data Collection Phase-- In this phase data party emulate sending mM trails of oblivious flajolet martin structures with N computation parties
	* @Return of the phase
	*     @bool           the phase successfully is running or not
	*/
	bool data_collection(int socket_num, std::vector<std::vector<Share<gfp>>>& share_FS);
	// auxilary protocol for Protocol 4 : get the share of the 0 elements
	void share_zero_array(int socket_num, std::vector<Share<gfp>>& share_Z);

  	/*
 	* @Protocol: data party run this protocol to distribute the share of a plaintext x (in gfp field) to all parties
 	* @Params of the protocol share_value
 	* 		@x 						the data party to share value shares (we use in typo: gfp)
 	* 		@socket_num 			the marking of the current data partie, make sure run the get_connection_socket sucessfully first
 	* @Return of the protocol share_value
 	* 		@bool 					the protocol successfully is running or not
 	*/
 	template <class T>
 	bool share_value(Share<T>& x, int socket_num);

 	//bucket transfer the whole Oblivious FlajoletMartin shared-value data
 	template <class T>
 	bool share_value(std::vector<std::vector<Share<T>>>& des_vec, int socket_num);

 	// Send an octetStream to data player 
 	void send_to(int data_player,const octetStream& o,bool donthash=false) const;

 	// receive an octetstream from data player 
  	void receive_player(int data_player,octetStream& o,bool donthash=false) const;

};


#endif /* COMPUTATION_PARTY_H_ */
