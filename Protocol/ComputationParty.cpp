// (C) 2018 University of NKU. Free for used

/*
 * ComputationParty.cpp
 *
 */
#include "Protocol/ComputationParty.h"
#include "Polynomial.h"

#include "Math/Share.h"
#include "Math/bigint.h"
#include "Auth/fake-stuff.h"
#include "Tools/ezOptionParser.h"
#include "Exceptions/Exceptions.h"

#include "Tools/octetStream.h"
#include "Tools/int.h"


#include <vector>
#include <array>
#include <numeric>
#include <sstream>
#include <math.h>  

#include <sys/time.h>

template<class T>
struct CP_thread_bucket_mod2 {
   CP* obj;
   std::vector<std::vector<Share<T>>>*  src_vec;
   std::vector<std::vector<Share<T>>>*  t0;
   int m_start;
   int m_internal;
   int t_id;
} ;


template<class T>
struct CP_thread_bucket_bExtractZ {
   CP* obj;
   std::vector<std::vector<Share<T>>>*  src_fm;
   Share<T>* z;
   int m_start;
   int m_internal;
   int t_id;
   int n_bucket;
} ;

typedef struct CP_thread_sender {
   CP* obj;
   octetStream* o;
   int t_id;
   int player_no;
    bool send_func;
} thread_sender_info;

template<class T>
using thread_lipmaaZero_info = CP_thread_bucket_mod2<T>;

template<class T>
using thread_bExtractZ_info = CP_thread_bucket_bExtractZ<T>;


CP::CP(int argc, const char** argv)
{
    ez::ezOptionParser opt;

    opt.syntax = "./CP.x [OPTIONS]\n";
    opt.example = "./CP.x -lgp 64 -np 2 -p 0 -x 4 -ndp 1 \n";
    opt.add(
        "5", // Default.
        0, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "number of computation parties (default: 5)", // Help description.
        "-np", // Flag token.
        "--number_parties" // Flag token.
    );
    opt.add(
        "1", // Default.
        0, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "number of data parties (default: 1)", // Help description.
        "-ndp", // Flag token.
        "--number_data_parties" // Flag token.
    );
    opt.add(
        "100", // Default.
        0, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "number of the party (starting from 0)", // Help description.
        "-p", // Flag token.
        "--my_num" // Flag token.
    );
    opt.add(
        "128", // Default.
        0, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "Bit length of GF(p) field (default: 128)", // Help description.
        "-lgp", // Flag token.
        "--lgp" // Flag token.
    );
    opt.add(
        "40", // Default.
        0, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "security parameter (default: 40)", // Help description.
        "-sec", // Flag token.
        "--security-parameter" // Flag token.
    );
    opt.add(
            "", // Default.
            0, // Required?
            1, // Number of args expected.
            0, // Delimiter if expecting multiple args.
            "Directory containing the data (default: " PREP_DIR "<nparties>-<lgp>-<lg2>", // Help description.
            "-d", // Flag token.
            "--dir" // Flag token.
    );
    opt.add(
          "", // Default.
          0, // Required?
          1, // Number of args expected.
          0, // Delimiter if expecting multiple args.
          "Port to listen on (default: port number base + player number)", // Help description.
          "-mp", // Flag token.
          "--my-port" // Flag token.
    );
    opt.add(
          "5000", // Default.
          0, // Required?
          1, // Number of args expected.
          0, // Delimiter if expecting multiple args.
          "Port number base to attempt to start connections from (default: 5000)", // Help description.
          "-pn", // Flag token.
          "--portnumbase" // Flag token.
    );
    opt.add(
          "localhost", // Default.
          0, // Required?
          1, // Number of args expected.
          0, // Delimiter if expecting multiple args.
          "Host where Server.x is running to coordinate startup (default: localhost).", // Help description.
          "-h", // Flag token.
          "--hostname" // Flag token.
    );
    opt.add(
            "1", // Default.
            0, // Required?
            1, // Number of args expected.
            0, // Delimiter if expecting multiple args.
            "number of the threads(default 1)", // Help description.
            "-nt", // Flag token.
            "--number_threads" // Flag token.
    );
    opt.add(
            "500", // Default.
            0, // Required?
            1, // Number of args expected.
            0, // Delimiter if expecting multiple args.
            "number of FlajoletMartin trails", // Help description.
            "-oM", // Flag token.
            "--number_OFS_trails" // Flag token.
    );
    opt.add(
            "20000", // Default.
            0, // Required?
            1, // Number of args expected.
            0, // Delimiter if expecting multiple args.
            "number of unique items", // Help description.
            "-uN", // Flag token.
            "--number_unique_itmes" // Flag token.
    );
    opt.add(
            "18", // Default.
            0, // Required?
            1, // Number of args expected.
            0, // Delimiter if expecting multiple args.
            "the d of the RandExp generated by offline file", // Help description.
            "-L", // Flag token.
            "--length_RExp" // Flag token.
    );
    opt.add(
        "1", // Default.
        0, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "whether in local network environment, 0 means not and 1 means true", // Help description.
        "-lan", // Flag token.
        "--in lan envs" // Flag token.
    );




    opt.parse(argc, argv);

    string usage;
    string hostname;
    int in_local_envs;
    opt.get("--my_num")->getInt(mynum);
    opt.get("--lgp")->getInt(lgp);
    opt.get("--number_parties")->getInt(nparties);
    opt.get("--number_data_parties")->getInt(ndparties);
    opt.get("--portnumbase")->getInt(pnbase);
    opt.get("--hostname")->getString(hostname);
    opt.get("--number_OFS_trails")->getInt(oM);
    opt.get("--number_unique_itmes")->getInt(uN);
    opt.get("--number_threads")->getInt(nthreads);
    opt.get("--length_RExp")->getInt(l_RExp);
    opt.get("--security-parameter")->getInt(sec);
    opt.get("--in lan envs")->getInt(in_local_envs);

    local_envs = bool(in_local_envs);
    
    if(mynum){
        playerone = false;
    }else{
        playerone = true;
    }
    l_RExp+=1; // for R^{-1}
    tau = lgp - sec;     // sec = 40

    if(mynum>10)
        throw runtime_error("the player number have not been set!");
    if (opt.isSet("--dir"))
    {
        opt.get("--dir")->getString(PREP_DATA_PREFIX);
        PREP_DATA_PREFIX += "/";
    }
    else
        PREP_DATA_PREFIX = get_prep_dir(nparties, lgp, 40);


    ez::OptionGroup* mp_opt = opt.get("--my-port");
    if (mp_opt->isSet)
      mp_opt->getInt(my_port);
    else
      my_port = Names::DEFAULT_PORT;


    CommsecKeysPackage *keys = NULL;
    vector<Names> playerNames(nthreads*communication_multiplier);
    // // ******************************** Network Part ***********************************
    for(int i=0; i<nthreads*communication_multiplier; i++){
        playerNames[i].init(mynum, pnbase+port_increase*i, my_port, hostname.c_str());
        playerNames[i].set_keys(keys);
    }
    

    // ******************************** File Handler  ************************************
    read_setup(PREP_DATA_PREFIX);
    dataF = new Data_Files(mynum, nparties, PREP_DATA_PREFIX);
    init_offline_data();
    init_Polynomial_coef(PREP_DATA_PREFIX, tau);
    init_wait_to_check_buffer();

    // ******************************** Mac Key ************************************
    /* Find number players and MAC keys etc*/
    char filename[1024];
    gfp pp; 
    keyp.assign_zero();
    alphai.assign_zero();
    int N=1;
    ifstream inpf;
    for (int i= 0; i < nparties; i++)
    {
        sprintf(filename, (PREP_DATA_PREFIX + "Player-MAC-Keys-P%d").c_str(), i);
        inpf.open(filename);
        if (inpf.fail()) { throw file_error(filename); }
        inpf >> N;
        pp.input(inpf,true);
        //cout << " Key " << i << "\t p: " << pp  << endl;
        keyp.add(pp);
        inpf.close();
        if(i==mynum)
            alphai.add(pp);
    }
    //cout << "--------------\n";
    //cout << "Final Keys :\t p: " << keyp << endl;

    // p2p whole connected
    //player = new Player(playerNames[0], 0);
    thread_player.resize(nthreads*communication_multiplier);
    commu_timer.resize(nparties*nthreads*communication_multiplier);
    for(int i=0; i<nthreads*communication_multiplier; i++){
        thread_player[i] = new Player(playerNames[i], 0);
    }
    player = thread_player[0];
    CP_sockets.resize(nthreads*communication_multiplier);
    for(int i=0; i<nthreads*communication_multiplier; i++){
        CP_sockets[i].resize(nparties);
        for(int k=0; k<nparties; k++){
            CP_sockets[i][k] = thread_player[i]->socket(k);
        }
    }


    //******************************** Data Collection Phase *****************************
    DP_connector dp_handle;
    dp_handle.init(mynum, pnbase, my_port,ndparties);
    dp_handle.key_init(alphai,keyp);
    dp_handle.player_init(player);
    dp_handle.n_threads_init(nthreads);
    dp_handle.params_ofs_init(oM, uN);
    dp_handle.start(dataF,share_FS);

    int w = share_FS[0].size();
    share_BFS.resize(oM);
    for(int m=0; m<oM; m++){
        share_BFS[m].resize(w);
    }

    // multithread 
    mutex_go = PTHREAD_MUTEX_INITIALIZER;
    local_Go.resize(nthreads);
    mutex_local_go.resize(nthreads);
    for(int i=0; i<nthreads; i++){
        mutex_local_go[i] = PTHREAD_MUTEX_INITIALIZER;
    }
}

CP::~CP(){
    delete dataF;
    //delete player;
    player = 0;
    for(int i=0; i<nthreads; i++){
        delete thread_player[i];
    }  
}
void CP::init_offline_data(){
    std::stringstream ss;

    rand_pool.resize(nthreads);
    rand2_pool.resize(nthreads);
    triple_pool.resize(nthreads);
    randExp_pool.resize(nthreads);
    randExp64_pool.resize(nthreads);

    for(int i=0; i<nthreads; i++){
        //read share of rand value
        ss.str("");
        ss << PREP_DATA_PREFIX << "/Rands-p-P" << mynum;
        if (i){
            ss << "-" << i;
        }
        rand_pool[i].open(ss.str().c_str(), ios::in | ios::binary);
        //read share of rand2 value
        ss.str("");
        ss << PREP_DATA_PREFIX << "/Rand2s-p-P" << mynum;
        if (i){
            ss << "-" << i;
        }
        rand2_pool[i].open(ss.str().c_str(), ios::in | ios::binary);
        //read share of triple value
        ss.str("");
        ss << PREP_DATA_PREFIX << "/Triples-p-P" << mynum;
        if (i){
            ss << "-" << i;
        }
        triple_pool[i].open(ss.str().c_str(), ios::in | ios::binary);
        //read share of randExp value 
        ss.str("");
        // ss << PREP_DATA_PREFIX << "/RandExps-p-18-P" << mynum;
        ss << PREP_DATA_PREFIX << "/RandExps-18-p-P" << mynum;
        if (i){
            ss << "-" << i;
        }
        randExp_pool[i].open(ss.str().c_str(), ios::in | ios::binary);
        //read share of randExp64 value 
        ss.str("");
        // ss << PREP_DATA_PREFIX << "/RandExps-p-"<<tau<<"-P" << mynum;
        ss << PREP_DATA_PREFIX << "/RandExps-p-P" << mynum;
        if (i){
            ss << "-" << i;
        }
        randExp64_pool[i].open(ss.str().c_str(), ios::in | ios::binary);
    }

}

void CP::close_offline_data_stream(){
    for(int i=0; i<nthreads; i++){
        rand_pool[i].close();
        rand2_pool[i].close();
        triple_pool[i].close();
        randExp_pool[i].close();
        randExp64_pool[i].close();
    }
    
}

void CP::init_Polynomial_coef(std::string folder_name, int tau){
    using namespace PSUCA::Polynomial;

    folder_name = folder_name + "../Polynomial-" + to_string(lgp) + "-" + to_string(mynum) + "/";
    bigint pr_p = gfp::pr();
    std::vector<std::pair<int,int>> points;
    int ofs_w = ceil(log2 (uN))+4;
    poly_cdb_f.clear();
    poly_cdb_f.resize(ofs_w);

    // If not prepared the polynomial, generate it first
    bool generated = true;
    try{
        // polyzero
        {
            Interpolatation_Polynomial poly_f(pr_p, folder_name, tau, poly_zero);
            poly_f.output(polyzero_f); 
        }
        
        // cdb
        int w = ceil(double(ofs_w)/2);
        while(w > 2){
            Interpolatation_Polynomial poly_f(pr_p, folder_name, w, poly_cdb);
            poly_f.output(poly_cdb_f[w]); 
            w = ceil(double(w)/2);
            

        }

        // lookup
        {
            Interpolatation_Polynomial poly_f(pr_p, folder_name, 3, poly_lookup);
            poly_f.output(clookup_d3); 
        }
        
        {
            Interpolatation_Polynomial poly_f(pr_p, folder_name, 2, poly_lookup);
            poly_f.output(clookup_d2); 
        }
    }catch(file_missing e){
        generated = false; 
    }

    if(!generated){
        {
            Interpolatation_Polynomial::write_points(folder_name, tau, poly_zero);
            Interpolatation_Polynomial::load_points(folder_name, tau, poly_zero, points);
            Interpolatation_Polynomial poly_f(pr_p, points.size(), points);
            poly_f.write_coeff(folder_name, poly_zero);
        }
        

        // cdb
        int w = ceil(log2 (uN))+4;
        w = ceil(double(w)/2);
        while(w > 2){
            Interpolatation_Polynomial::write_points(folder_name, w, poly_cdb);
            Interpolatation_Polynomial::load_points(folder_name, w, poly_cdb, points);
            Interpolatation_Polynomial poly_f(pr_p, points.size(), points);
            poly_f.write_coeff(folder_name, poly_cdb);
            //cout<<w<<" ";
            w = ceil(double(w)/2);
        }

        // lookup
        {
            Interpolatation_Polynomial::write_points(folder_name, 3, poly_lookup);
            Interpolatation_Polynomial::load_points(folder_name, 3, poly_lookup, points);
            Interpolatation_Polynomial poly_f(pr_p, points.size(), points);
            poly_f.write_coeff(folder_name, poly_lookup);
        }
        
        {
            Interpolatation_Polynomial::write_points(folder_name, 2, poly_lookup);
            Interpolatation_Polynomial::load_points(folder_name, 2, poly_lookup, points);
            Interpolatation_Polynomial poly_f(pr_p, points.size(), points);
            poly_f.write_coeff(folder_name, poly_lookup);
        }

        // reload to the polynomial
        // polyzero
        {
            Interpolatation_Polynomial poly_f(pr_p, folder_name, tau, poly_zero);
            poly_f.output(polyzero_f); 
        }
        
        // cdb
        w = ceil(double(ofs_w)/2);
        while(w > 2){
            Interpolatation_Polynomial poly_f(pr_p, folder_name, w, poly_cdb);
            poly_f.output(poly_cdb_f[w]); 
            w = ceil(double(w)/2);
        }

        // lookup
        {
            Interpolatation_Polynomial poly_f(pr_p, folder_name, 3, poly_lookup);
            poly_f.output(clookup_d3); 
        }
        
        {
            Interpolatation_Polynomial poly_f(pr_p, folder_name, 2, poly_lookup);
            poly_f.output(clookup_d2); 
        }
    }
}

template <class T>
void CP::lipmaa_zeroTest(const Share<T>& t, Share<T>& t0, int thread_num){
    std::vector<Share<T>> rb_group(lgp+1);
    std::vector<Share<T>> R_group(tau+1);
    std::vector<Share<T>> c(tau+1);
    T m, mac, tmp;
    bigint big_m, big_tmp, big_one(1);
    T T_negate_two(-2);

    Share<T> sh_m, sh_H;


    getRB_group(rb_group,thread_num);
    sh_m = rb_group[0]+t;
    open_and_check(sh_m,m,mac);

    to_bigint(big_m,m);
    for(int i=1; i<tau+1; i++){
        big_tmp =  big_m & 1;
        big_m = (big_m >> 1);
        sh_H += rb_group[i];
        if(big_tmp == big_one){
            sh_H += rb_group[i]*T_negate_two;
            if(mynum){
                sh_H.add(sh_H,1,false,alphai);
            }else{
                sh_H.add(sh_H,1,true,alphai);
            }
        }
    }
    
    getR_group(R_group,tau+1,thread_num,tau+1);
    dot_product(sh_H, R_group[0], sh_m,thread_num);  
    open_and_check(sh_m,m,mac);

    t0.assign_zero();
    tmp.assign(m);
    // get the x, x^2,..., x^d
    for(int i=1; i<tau+1; i++){
        c[i] = R_group[i]*tmp;
        tmp = tmp*m;
    }

    const vector<gfp>& polyzero_f = get_zerotest_PolyNomailFunc();

    for(int i=1; i<tau+1; i++){
        t0 += c[i]*polyzero_f[i];
    } 
}

template <class T>
void CP::bucket_lipmaa_zeroTest(const vector<vector<Share<T>>>& src_vec, std::vector<std::vector<Share<T>>>& t0, const int m_start, int m_internal, int thread_num, bool multi_thread){
    int ofs_M = src_vec.size();
    int ofs_w = src_vec[0].size();
    if(!ofs_M){
        return;
    }

    const vector<gfp>& polyzero_f = get_zerotest_PolyNomailFunc();
    m_internal = ofs_M < (m_start+m_internal) ? ofs_M - m_start : m_internal;
    int m_end = m_start+m_internal;

    int tmp_int = m_internal*ofs_w;
    int offset;
    vector<vector<Share<T>>> rb_group(tmp_int);
    vector<vector<Share<T>>> R_group(tmp_int);
    vector<Share<T>> sh_m(tmp_int), sh_H(tmp_int), R0(tmp_int);
    std::vector<T> m(tmp_int); //mac(tmp_int);

    tmp_int = lgp+1;
    for(int i=0; i<m_internal*ofs_w; i++){
        rb_group[i].resize(tmp_int);
        getRB_group(rb_group[i],thread_num);
        R_group[i].resize(tau+1);
        getR_group(R_group[i],tau+1,thread_num,tau+1);
        R0[i] = R_group[i][0];
    }

    Share<T> sh_tmp;

    T tmp;
    bigint big_m, big_tmp, big_one(1);
    T T_negate_two(-2);
       
    offset = 0;
    for(int k=m_start; k<m_end; k++){
        for(int i=0; i<ofs_w; i++){
            sh_m[offset] = rb_group[offset][0]+src_vec[k][i];
            offset++;
        }
    }

    _delay_open_and_check(sh_m, m, m_internal*ofs_w, 0, 0, thread_num);


    for(int j=0; j<m_internal*ofs_w; j++){
        to_bigint(big_m,m[j]);
        for(int i=1; i<tau+1; i++){
            big_tmp =  big_m & 1;
            big_m = (big_m >> 1);
            sh_H[j] += rb_group[j][i];
            if(big_tmp == big_one){
                sh_H[j] += rb_group[j][i]*T_negate_two;
                if(mynum){
                    sh_H[j].add(sh_H[j],1,false,alphai);
                }else{
                    sh_H[j].add(sh_H[j],1,true,alphai);
                }
            }
        }
    }

    
    dot_product(sh_H, R0, sh_m,thread_num);  
    _delay_open_and_check(sh_m, m, m_internal*ofs_w, 0, 0, thread_num);

    offset = 0;
    for(int k=m_start; k<m_end; k++){   
        for(int j=0; j<ofs_w; j++){
            T& m_tmp = m[offset];
            tmp.assign(m_tmp);
            t0[k][j].assign_zero();
            for(int i=1; i<tau+1; i++){
                sh_tmp = R_group[offset][i]*tmp;
                t0[k][j] += sh_tmp*polyzero_f[i];
                tmp = tmp*m_tmp;
            }
            offset++;
        }
    }


    if(multi_thread){
        pthread_mutex_lock(&mutex_go); 
        Go++;
        pthread_mutex_unlock(&mutex_go);
    }
}

template <class T>
void CP::bucket_lipmaa_zeroTest_bench(const vector<vector<Share<T>>>& src_vec, std::vector<std::vector<Share<T>>>& t0, const int m_start, int m_internal, int thread_num, bool multi_thread){
    int ofs_M = src_vec.size();
    int ofs_w = src_vec[0].size();
    if(!ofs_M){
        return;
    }

    const vector<gfp>& polyzero_f = get_zerotest_PolyNomailFunc();

    m_internal = ofs_M < (m_start+m_internal) ? ofs_M - m_start : m_internal;
    int m_end = m_start+m_internal;

    int tmp_int = m_internal*ofs_w;
    int offset;
    vector<Share<T>> sh_m(tmp_int), sh_H(tmp_int), R0(tmp_int);
    std::vector<T> m(tmp_int); //mac(tmp_int);

    vector<Share<T>> rb_group(lgp+1);
    vector<Share<T>> R_group(tau+1);

    getRB_group(rb_group,thread_num);
    getR_group(R_group,tau+1,thread_num,tau+1);
    for(int i=0; i<m_internal*ofs_w; i++){
        R0[i] = R_group[0];
    }

    Share<T> sh_tmp;

    T tmp;
    bigint big_m, big_tmp, big_one(1);
    T T_negate_two(-2);
       
    offset = 0;
    for(int k=m_start; k<m_end; k++){
        for(int i=0; i<ofs_w; i++){
            sh_m[offset++] = rb_group[0]+src_vec[k][i];
        }
    }


    _delay_open_and_check(sh_m, m, m_internal*ofs_w, 0, 0, thread_num);

    for(int j=0; j<m_internal*ofs_w; j++){
        to_bigint(big_m,m[j]);
        for(int i=1; i<tau+1; i++){
            big_tmp =  big_m & 1;
            big_m = (big_m >> 1);
            sh_H[j] += rb_group[i];
            if(big_tmp == big_one){
                sh_H[j] += rb_group[i]*T_negate_two;
                if(mynum){
                    sh_H[j].add(sh_H[j],1,false,alphai);
                }else{
                    sh_H[j].add(sh_H[j],1,true,alphai);
                }
            }
        }
    }

    dot_product(sh_H, R0, sh_m,thread_num);  
    _delay_open_and_check(sh_m, m, m_internal*ofs_w, 0, 0, thread_num);

    offset = 0;
    for(int k=m_start; k<m_end; k++){   
        for(int j=0; j<ofs_w; j++){
            T& m_tmp = m[offset];
            tmp.assign(m_tmp);
            t0[k][j].assign_zero();
            for(int i=1; i<tau+1; i++){
                sh_tmp = R_group[i]*tmp;
                t0[k][j] += sh_tmp*polyzero_f[i];
                tmp = tmp*m_tmp;
            }
            offset++;
        }
    }

    if(multi_thread){
        pthread_mutex_lock(&mutex_go); 
        Go++;
        pthread_mutex_unlock(&mutex_go);
    }
}

template <class T>
void CP::multi_thread_bucket_lipmaa_zeroTest(std::vector<std::vector<Share<T>>>& src_vec, std::vector<std::vector<Share<T>>>& t0){
    if(thread_player.size() == 0){
        throw runtime_error("The p2p socket has been expired");
    }

    int ofs_M = src_vec.size();
    if(!ofs_M){
        return;
    }

    std::vector<thread_lipmaaZero_info<T>> thread_lipmaaZero_data(nthreads);
    Go = 0;

    int m_start = 0;
    int m_internal = ceil(((double)ofs_M)/nthreads);
    //int m_internal = ceil(ofs_M/nthreads);

    pthread_t* t = new pthread_t[nthreads];
    for(int i=0; i<nthreads; i++){
        thread_lipmaaZero_data[i].obj = this;
        thread_lipmaaZero_data[i].src_vec = &src_vec;
        thread_lipmaaZero_data[i].t0 = &t0;
        thread_lipmaaZero_data[i].t_id = i;
        thread_lipmaaZero_data[i].m_start = m_start;
        m_start+=m_internal;
        thread_lipmaaZero_data[i].m_internal = m_internal;
        pthread_create(&t[i], NULL,thread_lipmaa_zeroTest, (void*) &thread_lipmaaZero_data[i]);
    }

    while(Go < nthreads){
        usleep(10);
    }

    delete t;
}

template <class T>
void CP::getR_group(std::vector<Share<T>>& R_group, const int len_group, int thread_num, int l_R){
    if(len_group > l_R){
        throw runtime_error("no valid RandExp files found, please generate it first or tune the input parameter of this program use -L");
    }

    if((int)R_group.size() < len_group){
        R_group.resize(len_group);
    }

    Share<T> tmp;
    // this function is ugly, but I don't have much time on it
    if(l_R==19){
        for(int i=0; i<l_RExp;i++){
            try{
                if(i<len_group){
                    R_group[i].input(randExp_pool[thread_num], false);
                }else{
                    tmp.input(randExp_pool[thread_num], false);
                }
                
            }catch( end_of_file e ){
                randExp_pool[thread_num].clear();
                randExp_pool[thread_num].seekg(0);
                if(i<len_group){
                    R_group[i].input(randExp_pool[thread_num], false);
                }else{
                    tmp.input(randExp_pool[thread_num], false);
                }
            }   
        }
    }

    if(l_R==tau+1){
        for(int i=0; i<tau+1;i++){
            try{
                if(i<len_group){
                    R_group[i].input(randExp64_pool[thread_num], false);
                }else{
                    tmp.input(randExp64_pool[thread_num], false);
                }
                
            }catch( end_of_file e ){
                randExp64_pool[thread_num].clear();
                randExp64_pool[thread_num].seekg(0);
                if(i<len_group){
                    R_group[i].input(randExp64_pool[thread_num], false);
                }else{
                    tmp.input(randExp64_pool[thread_num], false);
                }
            }   
        }
    }
}

template <class T>
void CP::getRB_group(std::vector<Share<T>>& rb_group, int thread_num){
    if((int)rb_group.size() < (lgp+1)){
        rb_group.resize(lgp+1);
    }

    rb_group[0].assign_zero();
    T tmp_int;
    tmp_int.assign(1);
    for(int j=1; j<lgp; j++){
        try{
            rb_group[j].input(rand2_pool[thread_num],false);
            rb_group[0] += rb_group[j]*tmp_int;
            tmp_int *= 2;
        }catch(end_of_file e){
            rand2_pool[thread_num].clear();
            rand2_pool[thread_num].seekg(0);
            rb_group[j].input(rand2_pool[thread_num],false);
            rb_group[0] += rb_group[j]*tmp_int;
            tmp_int *= 2;
        }
    }

}


const vector<gfp>& CP::get_zerotest_PolyNomailFunc(){
    if(polyzero_f.empty()){
        throw runtime_error("The zeroTest polynomial function has not been initialized");
    }else{
        return polyzero_f;
    }
}
const vector<gfp>& CP::get_CExtractZ_PolyNomailFunc(int d){
    switch(d){
        case 3:
        case 5:
        case 6:
        case 9:
        case 10:
        case 12:
        case 17:
            if(poly_cdb_f[d].empty()){
                throw runtime_error("The polynomial function has not been initialized");
            }else{
                return poly_cdb_f[d];
            }
            break;
        case 8:
            if(clookup_d3.empty()){
                throw runtime_error("The polynomial function has not been initialized");
            }else{
                return clookup_d3;
            }
            break;
        case 4:
            if(clookup_d2.empty()){
                throw runtime_error("The polynomial function has not been initialized");
            }else{
                return clookup_d2;
            }
            break;
        default:
            throw runtime_error("The polynomial function has not been initialized");
    }
    throw runtime_error("The CEXtract polynomial function has not been initialized");
    return polyzero_f; // just to pass around that treating warning as erro
}

template <class T>
void CP::bdb_CExtractZ(const std::vector<Share<T>>& b, Share<T>& y, const int d, int thread_num){
    vector<T> poly_f = get_CExtractZ_PolyNomailFunc(d);
    std::vector<Share<T>> R_group(d+1);
    std::vector<Share<T>> c(d+1);
    getR_group(R_group,d+1,thread_num);
    

    Share<T> sh_x,sh_m;

    // get the maping from 2^d number to [1,d] 0 for 1; otherwise 0
    for(int i=0; i<d; i++){
        sh_x+=b[i];
    }

    // multiply x*R^-1 then reveal it
    T m,mac,tmp;

    dot_product(sh_x, R_group[0], sh_m,thread_num);  
    open_and_check(sh_m,m,mac);


    y.assign_zero();
    tmp.assign(m);
    // get the x, x^2,..., x^d
    if(mynum){
        c[0].add(c[0],1,false,alphai);
    }else{
        c[0].add(c[0],1,true,alphai);
    }
    
    for(int i=1; i<d+1; i++){
        c[i] = R_group[i]*tmp;
        tmp = tmp*m;
    }

    for(int i=0; i<d+1; i++){
        y += c[i]*poly_f[i];
    }   
}

template <class T>
void CP::lookup_x_CExtractZ(const std::vector<Share<T>>& b, Share<T>& y, const int start, const int d , int thread_num ){
    int len_group = pow(2,d);
    std::vector<Share<T>> R_group(len_group);
    std::vector<Share<T>> c(len_group);
    getR_group(R_group,len_group,thread_num);
    vector<T> poly_f = get_CExtractZ_PolyNomailFunc(len_group);

    int tmp_int = 1;
    Share<T> sh_x, sh_m;
    T m,mac,tmp;

    int end =  start+d;
    sh_x.assign_zero();

    for(int i=end-1; i>=start; i--){
        sh_x += b[i]*tmp_int;
        tmp_int *= 2;
    } 

    // multiply x*R^-1 then reveal it
    dot_product(sh_x, R_group[0], sh_m,thread_num);  
    open_and_check(sh_m,m,mac);

    y.assign_zero();
    tmp.assign(m);
    // get the x, x^2,..., x^d
    if(mynum){
        c[0].add(c[0],1,false,alphai);
    }else{
        c[0].add(c[0],1,true,alphai);
    }
    for(int i=1; i<len_group; i++){
        c[i] = R_group[i]*tmp;
        tmp = tmp*m;
    }

    for(int i=0; i<len_group; i++){
        y += c[i]*poly_f[i];
    } 
}



template <class T>
void CP::dot_product(Share<T>& x, Share<T>& y, Share<T>& z, int thread_num){
    if(thread_player[thread_num] == 0){
        throw runtime_error("The p2p socket has been expired");
    }
    Share<T> a,b,c; // offline triple
    std::vector<Share<T>> transfer_v(2); // first is epsilon, second is rou;

    std::vector<T> transfer_vT(2); // first is epsilon, second is rou;
    std::vector<T> open_mac(2);

    try{
        a.input(triple_pool[thread_num], false);
        b.input(triple_pool[thread_num], false);
        c.input(triple_pool[thread_num], false);
    }catch( end_of_file e ){
        triple_pool[thread_num].clear();
        triple_pool[thread_num].seekg(0);

        a.input(triple_pool[thread_num], false);
        b.input(triple_pool[thread_num], false);
        c.input(triple_pool[thread_num], false);
    }
    

    transfer_v[0].sub(x,a);
    transfer_v[1].sub(y,b);

    thread_open_and_check(transfer_v,transfer_vT,open_mac,thread_num);

    z.assign(c);
    z += b*transfer_vT[0];
    z += a*transfer_vT[1];
    z.add(z,transfer_vT[0]*transfer_vT[1],playerone,alphai);
}

template <class T>
void CP::dot_product(std::vector<Share<T>>& x, std::vector<Share<T>>& y, std::vector<Share<T>>& z, int thread_num){
    if(thread_player[thread_num] == 0){
        throw runtime_error("The p2p socket has been expired");
    }

    if(x.size() != y.size() || x.size() > z.size()){
        throw runtime_error("unlegal dot product , two vector have different size");
    }

    int n_elements = x.size();

    
    std::vector<Share<T>> a(n_elements), b(n_elements), c(n_elements); // offline triple
    std::vector<Share<T>> transfer_v(2*n_elements); // even for epsilon, odd for rou;
    //std::vector<T> open_mac(2*n_elements);
    std::vector<T> transfer_vT(2*n_elements);

    for(int i=0; i<n_elements; i++){
        try{
            a[i].input(triple_pool[thread_num], false);
            b[i].input(triple_pool[thread_num], false);
            c[i].input(triple_pool[thread_num], false);
        }catch( end_of_file e ){
            triple_pool[thread_num].clear();
            triple_pool[thread_num].seekg(0);

            a[i].input(triple_pool[thread_num], false);
            b[i].input(triple_pool[thread_num], false);
            c[i].input(triple_pool[thread_num], false);
        }

        transfer_v[i*2].sub(x[i],a[i]);
        transfer_v[i*2+1].sub(y[i],b[i]);
    }

    //thread_open_and_check(transfer_v,transfer_vT,open_mac,thread_num);
    _delay_open_and_check(transfer_v,transfer_vT,2*n_elements,0,0,thread_num);

    for(int i=0; i<n_elements; i++){
        z[i].assign(c[i]);
        z[i] += b[i]*transfer_vT[i*2];
        z[i] += a[i]*transfer_vT[i*2+1];
        z[i].add(z[i],transfer_vT[i*2]*transfer_vT[i*2+1],playerone,alphai);
    }
}

/*
*   This funciton is not memory safe at all, just for benchmark use. Since I don't want to copy the same value to different vector
*/
template <class T>
void CP::_dot_product(std::vector<Share<T>>& x, std::vector<Share<T>>& y, std::vector<Share<T>>& z, \
                     int n_elements, int start_x, int start_y, int start_z, int thread_num){
    
    std::vector<Share<T>> a(n_elements), b(n_elements), c(n_elements); // offline triple
    std::vector<Share<T>> transfer_v(2*n_elements); // even for epsilon, odd for rou;
    // std::vector<T> open_mac(2*n_elements);
    std::vector<T> transfer_vT(2*n_elements);

    for(int i=0; i<n_elements; i++){
        try{
            a[i].input(triple_pool[thread_num], false);
            b[i].input(triple_pool[thread_num], false);
            c[i].input(triple_pool[thread_num], false);
        }catch( end_of_file e ){
            triple_pool[thread_num].clear();
            triple_pool[thread_num].seekg(0);

            a[i].input(triple_pool[thread_num], false);
            b[i].input(triple_pool[thread_num], false);
            c[i].input(triple_pool[thread_num], false);
        }

        transfer_v[i*2].sub(x[start_x+i],a[i]);
        transfer_v[i*2+1].sub(y[start_y+i],b[i]);
    }

    //thread_open_and_check(transfer_v,transfer_vT,open_mac,thread_num);
    delay_open_and_check(transfer_v,transfer_vT,thread_num);

    int offset = start_z; 
    for(int i=0; i<n_elements; i++){
        z[offset].assign(c[i]);
        z[offset] += b[i]*transfer_vT[i*2];
        z[offset] += a[i]*transfer_vT[i*2+1];
        z[offset].add(z[offset],transfer_vT[i*2]*transfer_vT[i*2+1],playerone,alphai);
        offset++;
    }
}

template <class T>
void CP::_dot_product(Share<T>& mask, std::vector<Share<T>>& y, std::vector<Share<T>>& z, \
                     int n_elements, int start_y , int start_z , int thread_num){
    std::vector<Share<T>> a(n_elements), b(n_elements), c(n_elements); // offline triple
    std::vector<Share<T>> transfer_v(2*n_elements); // even for epsilon, odd for rou;
    // std::vector<T> open_mac(2*n_elements);
    std::vector<T> transfer_vT(2*n_elements);

    int offset = start_y; 
    for(int i=0; i<n_elements; i++){
        try{
            a[i].input(triple_pool[thread_num], false);
            b[i].input(triple_pool[thread_num], false);
            c[i].input(triple_pool[thread_num], false);
        }catch( end_of_file e ){
            triple_pool[thread_num].clear();
            triple_pool[thread_num].seekg(0);

            a[i].input(triple_pool[thread_num], false);
            b[i].input(triple_pool[thread_num], false);
            c[i].input(triple_pool[thread_num], false);
        }

        transfer_v[i*2].sub(mask,a[i]);
        transfer_v[i*2+1].sub(y[offset++],b[i]);
    }

    //thread_open_and_check(transfer_v,transfer_vT,open_mac,thread_num);
    delay_open_and_check(transfer_v,transfer_vT,thread_num);

    offset = start_z; 
    for(int i=0; i<n_elements; i++){
        z[offset].assign(c[i]);
        z[offset] += b[i]*transfer_vT[i*2];
        z[offset] += a[i]*transfer_vT[i*2+1];
        z[offset].add(z[offset],transfer_vT[i*2]*transfer_vT[i*2+1],playerone,alphai);
        offset++;
    }
}

template <class T>
void CP::_dot_product(vector<Share<T>>& mask, vector<std::vector<Share<T>>>& y, vector<std::vector<Share<T>>>& z, \
                     int n_elements, int start_y , int start_z , int thread_num){
    int len = mask.size();
    int num = len*n_elements;
    std::vector<Share<T>> a(num), b(num), c(num); // offline triple
    std::vector<Share<T>> transfer_v(2*num); // even for epsilon, odd for rou;
    std::vector<T> transfer_vT(2*num);

    int out_offset = 0, offset;
    for(int k=0; k<len; k++){
        offset = start_y; 
        for(int i=0; i<n_elements; i++){
            try{
                a[out_offset].input(triple_pool[thread_num], false);
                b[out_offset].input(triple_pool[thread_num], false);
                c[out_offset].input(triple_pool[thread_num], false);
            }catch( end_of_file e ){
                triple_pool[thread_num].clear();
                triple_pool[thread_num].seekg(0);

                a[out_offset].input(triple_pool[thread_num], false);
                b[out_offset].input(triple_pool[thread_num], false);
                c[out_offset].input(triple_pool[thread_num], false);
            }

            transfer_v[2*out_offset].sub(mask[k],a[out_offset]);
            transfer_v[2*out_offset+1].sub(y[k][offset++],b[out_offset]);
            out_offset++;
        }
    }

    //thread_open_and_check(transfer_v,transfer_vT,open_mac,thread_num);
    delay_open_and_check(transfer_v,transfer_vT,thread_num);
   
    out_offset = 0;
    for(int k=0; k<len; k++){
        offset = start_z; 
        for(int i=0; i<n_elements; i++){
            z[k][offset].assign(c[out_offset]);
            z[k][offset] += b[out_offset]*transfer_vT[out_offset*2];
            z[k][offset] += a[out_offset]*transfer_vT[out_offset*2+1];
            z[k][offset].add(z[k][offset],transfer_vT[out_offset*2]*transfer_vT[out_offset*2+1],playerone,alphai);
            offset++;
            out_offset++;
        }
    }
    
}


template <class T>
void CP::bExtractZ(vector<Share<T>> fm_sketch, Share<T>& z){
    z.assign(fm_sketch[0]);
    int ofs_w = fm_sketch.size();

    for(int i=1; i<ofs_w; i++){
        dot_product(fm_sketch[i-1],fm_sketch[i],fm_sketch[i]);
        z+= fm_sketch[i];
    }
}

template <class T>
void CP::bucket_bExtractZ(const std::vector<std::vector<Share<T>>>& src_fm, Share<T>& Z, const int m_start, int m_internal, int thread_num, bool multi_thread){
    int ofs_M = src_fm.size();
    int ofs_w = src_fm[0].size();
    if(!ofs_M){
        return;
    }
    m_internal = ofs_M < (m_start+m_internal) ? ofs_M - m_start : m_internal;
    Z.assign_zero();

    std::vector<std::vector<Share<T>>> fm_sketch(ofs_w);
    for(int i=0; i<ofs_w; i++){
        fm_sketch[i].resize(m_internal);
    }
    int m_end = m_start+m_internal;
    int m_index;
    for(int m=m_start; m<m_end; m++){
        m_index = m-m_start;
        for(int i=0; i<ofs_w; i++){
            fm_sketch[i][m_index] = src_fm[m][i];
        }
        Z += fm_sketch[0][m_index];
    }

    for(int i=1; i<ofs_w; i++){
        dot_product(fm_sketch[i-1],fm_sketch[i],fm_sketch[i],thread_num);
        for(int j=0; j<m_internal; j++){
            Z += fm_sketch[i][j];
        }
    }

    if(multi_thread){
        pthread_mutex_lock(&mutex_go); 
        Go++;
        pthread_mutex_unlock(&mutex_go);
    }
}


template <class T>
void CP::CExtractZ_basic(const std::vector<Share<T>>& src_fm, Share<T>& z){
    /*
    *   PreProcessing data and dataStructure
    */
    z.assign_zero();
    int ofs_w = src_fm.size();
    int end = ofs_w;
    int t = ceil(ofs_w/2.0);
    std::vector<Share<T>> fm_sketch(2*t);

    // bit-wise negate
    Share<T> sh_tmp(1,mynum,alphai);
    Share<T> B0;
    T score;
    for(int i=0; i<ofs_w; i++){
        fm_sketch[i].sub(sh_tmp,src_fm[i]);
    }
    for(int i=ofs_w; i<2*t; i++){
        fm_sketch[i].assign(sh_tmp);
    }

    // temp variable
    std::vector<T> a_array(2*t),mac_array(2*t);
    int offset;

    /*
    *   Start Protocol
    */
    while(end>3){
        bdb_CExtractZ(fm_sketch,B0,t);
        score.assign(t);
        z+=score*B0;

        _dot_product(B0,fm_sketch,fm_sketch,t,t,t);

        offset = t;
        for(int i=0; i<t; i++){
            fm_sketch[i]+=fm_sketch[offset++];
        }
        end = t;
        t = ceil(end/2.0);
        for(int i=end; i<2*t; i++){
            fm_sketch[i].assign(sh_tmp);
        }
    }

    lookup_x_CExtractZ(fm_sketch,sh_tmp,0,end);
    z+=sh_tmp;
}


template <class T>
void CP::bucket_CExtractZ_bench(const std::vector<std::vector<Share<T>>>& src_fm, Share<T>& z, const int m_start, int m_internal, int thread_num, bool multi_thread){
    /*
    *   PreProcessing data and dataStructure
    */
    int ofs_M = src_fm.size();
    int ofs_w = src_fm[0].size();
    if(!ofs_M){
        return;
    }
    m_internal = ofs_M < (m_start+m_internal) ? ofs_M - m_start : m_internal;
    z.assign_zero();
    int end = ofs_w;
    int t = ceil(ofs_w/2.0);
    int tmp_int = 2*t;
    std::vector<std::vector<Share<T>>> fm_sketch(m_internal);
    for(int m=0; m<m_internal; m++){
        fm_sketch[m].resize(tmp_int);
    }

    // temp variable
    //std::vector<T> a_array(2*t),mac_array(2*t);
    int offset;
    tmp_int = ceil(double(m_internal)/communication_multiplier) * communication_multiplier;
    vector<Share<T>> sh_x(tmp_int),sh_m(tmp_int),R0(tmp_int);
    Share<T> sh_tmp;
    T mac,tmp;
    std::vector<T> m(tmp_int);

    // bit-wise negate
    Share<T> sh_one(1,mynum,alphai);
    vector<Share<T>> B0(m_internal);
    T score;
    tmp_int = m_start;
    for(int k=0; k<m_internal; k++){
        for(int i=0; i<ofs_w; i++){
            fm_sketch[k][i].sub(sh_one,src_fm[tmp_int][i]);
        }
        for(int i=ofs_w; i<2*t; i++){
            fm_sketch[k][i].assign(sh_one);
        }
        tmp_int++;
    }
    
    vector<vector<Share<T>>> R_group(m_internal);
    for(int k=0; k<m_internal; k++){
        R_group[k].resize(t+1);
        getR_group(R_group[k],t+1,thread_num);
        R0[k] = R_group[k][0];
    }

    /*
    * variable used for transfering
    */
    tmp_int = t*m_internal;
    std::vector<Share<T>> triple_a(tmp_int), triple_b(tmp_int), triple_c(tmp_int); // offline triple
    std::vector<Share<T>> transfer_v(2*tmp_int); // even for epsilon, odd for rou;
    std::vector<T> transfer_vT(2*tmp_int);

    /*
    *   Start Protocol
    */
    while(end>3){
        const vector<T>& poly_f = get_CExtractZ_PolyNomailFunc(t);  
        for(int k=0; k<m_internal; k++){
            sh_x[k].assign_zero();
            for(int i=0; i<t; i++){
                sh_x[k]+=fm_sketch[k][i];
            }
        }
        
        /*
        *   dot_product(sh_x, R0, sh_m,thread_num);  
        */
        /******************************************START to Implement dot_product*************************************/
        for(int i=0; i<m_internal; i++){
            try{
                triple_a[i].input(triple_pool[thread_num], false);
                triple_b[i].input(triple_pool[thread_num], false);
                triple_c[i].input(triple_pool[thread_num], false);
            }catch( end_of_file e ){
                triple_pool[thread_num].clear();
                triple_pool[thread_num].seekg(0);

                triple_a[i].input(triple_pool[thread_num], false);
                triple_b[i].input(triple_pool[thread_num], false);
                triple_c[i].input(triple_pool[thread_num], false);
            }

            transfer_v[i*2].sub(sh_x[i],triple_a[i]);
            transfer_v[i*2+1].sub(R0[i],triple_b[i]);
        }

        
        
        _delay_open_and_check(transfer_v,transfer_vT,2*m_internal,0,0,thread_num);
        

        for(int i=0; i<m_internal; i++){
            sh_m[i].assign(triple_c[i]);
            sh_m[i] += triple_b[i]*transfer_vT[i*2];
            sh_m[i] += triple_a[i]*transfer_vT[i*2+1];
            sh_m[i].add(sh_m[i],transfer_vT[i*2]*transfer_vT[i*2+1],playerone,alphai);
        }
        /******************************************END to Implement dot_product*************************************/
        _delay_open_and_check(sh_m,m,m_internal,0,0,thread_num);
        

        for(int k=0; k<m_internal; k++){
            B0[k].assign(sh_one);
            B0[k] += sh_x[k]*poly_f[1];
        }

        
        // get the x, x^2,..., x^d
        for(int k=0; k<m_internal; k++){
            T& tmp_m = m[k];
            tmp.assign(tmp_m);
            tmp = tmp*tmp_m;
            for(int i=2; i<t+1; i++){
                sh_tmp = R_group[k][i]*tmp;
                B0[k] += sh_tmp*poly_f[i];
                tmp = tmp*tmp_m;
            } 
        }
        
        score.assign(t);
        for(int k=0; k<m_internal; k++){
            z+=score*B0[k];
        }
        
        /*
        *   _dot_product(B0,fm_sketch,fm_sketch,t,t,t);
        */
        /******************************************START to Implement _dot_product*************************************/
        tmp_int = 0;
        for(int k=0; k<m_internal; k++){
            offset = t; 
            Share<T>& tmp_B0 = B0[k];
            for(int i=0; i<t; i++){
                try{
                    triple_a[tmp_int].input(triple_pool[thread_num], false);
                    triple_b[tmp_int].input(triple_pool[thread_num], false);
                    triple_c[tmp_int].input(triple_pool[thread_num], false);
                }catch( end_of_file e ){
                    triple_pool[thread_num].clear();
                    triple_pool[thread_num].seekg(0);

                    triple_a[tmp_int].input(triple_pool[thread_num], false);
                    triple_b[tmp_int].input(triple_pool[thread_num], false);
                    triple_c[tmp_int].input(triple_pool[thread_num], false);
                }


                transfer_v[2*tmp_int].sub(tmp_B0,triple_a[tmp_int]);
                transfer_v[2*tmp_int+1].sub(fm_sketch[k][offset++],triple_b[tmp_int]);
                tmp_int++;
            }
        }

        _delay_open_and_check(transfer_v,transfer_vT,m_internal*2*t,0,0,thread_num);
       
        tmp_int = 0;
        for(int k=0; k<m_internal; k++){
            offset = t; 
            for(int i=0; i<t; i++){
                Share<T>& tmp_z = fm_sketch[k][offset];
                tmp_z.assign(triple_c[tmp_int]);
                tmp_z += triple_b[tmp_int]*transfer_vT[tmp_int*2];
                tmp_z += triple_a[tmp_int]*transfer_vT[tmp_int*2+1];
                tmp_z.add(tmp_z,transfer_vT[tmp_int*2]*transfer_vT[tmp_int*2+1],playerone,alphai);
                offset++;
                tmp_int++;
            }
        }
        /******************************************END to Implement _dot_product*************************************/

        for(int k=0; k<m_internal; k++){
            offset = t;
            for(int i=0; i<t; i++){
                fm_sketch[k][i]+=fm_sketch[k][offset++];
            }
        }
        
        end = t;
        t = ceil(end/2.0);

        for(int k=0; k<m_internal; k++){
            for(int i=end; i<2*t; i++){
                fm_sketch[k][i].assign(sh_one);
            }
        }    
    }

    //lookup_x_CExtractZ(fm_sketch,sh_tmp,0,end);
    int len_group = pow(2,end);
    for(int k=0; k<m_internal; k++){
        sh_x[k].assign_zero();
    }

    tmp_int = 1;
    for(int i=end-1; i>=0; i--){
        for(int k=0; k<m_internal; k++){
            sh_x[k] += fm_sketch[k][i]*tmp_int;
        }
        tmp_int *= 2;
    }


    const vector<T>& poly_f = get_CExtractZ_PolyNomailFunc(len_group);
    /*
    *   dot_product(sh_x, R0, sh_m,thread_num);  
    */
    /******************************************START to Implement dot_product*************************************/
    for(int i=0; i<m_internal; i++){
        try{
            triple_a[i].input(triple_pool[thread_num], false);
            triple_b[i].input(triple_pool[thread_num], false);
            triple_c[i].input(triple_pool[thread_num], false);
        }catch( end_of_file e ){
            triple_pool[thread_num].clear();
            triple_pool[thread_num].seekg(0);

            triple_a[i].input(triple_pool[thread_num], false);
            triple_b[i].input(triple_pool[thread_num], false);
            triple_c[i].input(triple_pool[thread_num], false);
        }

        transfer_v[i*2].sub(sh_x[i],triple_a[i]);
        transfer_v[i*2+1].sub(R0[i],triple_b[i]);
    }

    _delay_open_and_check(transfer_v,transfer_vT,m_internal*2,0,0,thread_num);


    for(int i=0; i<m_internal; i++){
        sh_m[i].assign(triple_c[i]);
        sh_m[i] += triple_b[i]*transfer_vT[i*2];
        sh_m[i] += triple_a[i]*transfer_vT[i*2+1];
        sh_m[i].add(sh_m[i],transfer_vT[i*2]*transfer_vT[i*2+1],playerone,alphai);
    }
    /******************************************END to Implement dot_product*************************************/
    _delay_open_and_check(sh_m,m,m_internal,0,0,thread_num);

    Share<T> sh_end(end,mynum,alphai);

    for(int k=0; k<m_internal; k++){
        z += sh_end;
        z += sh_x[k]*poly_f[1];

        tmp.assign(m[k]);
        tmp = tmp*m[k];

        for(int i=2; i<len_group; i++){
            sh_tmp = R_group[k][i]*tmp;
            z += sh_tmp*poly_f[i];
            tmp = tmp*m[k];
        }
    }
    
 
    if(multi_thread){
        pthread_mutex_lock(&mutex_go); 
        Go++;
        pthread_mutex_unlock(&mutex_go);
    }
}


void* CP::thread_CExtractZ(void* arg){
    thread_bExtractZ_info<gfp>* data = static_cast<thread_bExtractZ_info<gfp>*>(arg);
    CP* obj = data->obj;
    obj->bucket_CExtractZ_bench(*(data->src_fm),*(data->z),data->m_start,data->m_internal,data->t_id,true);

    pthread_exit(NULL);
}

void CP::benchmark(){
    cout << "************** Start data Aggregation phase **************" << endl;

    struct timeval t1;
    struct timeval t2;
    double cpu_time_used;

    size_t send_bytes;
    double KB_bytes;

    gettimeofday(&t1, NULL);
    mergeShare();
    gettimeofday(&t2, NULL);
    cpu_time_used = (double)(t2.tv_sec-t1.tv_sec)*1000+(double)(t2.tv_usec-t1.tv_usec)/1000;
    //printf("mergeShare\n%d\t%d\t%f\n",uN,oM,cpu_time_used);
    printf("mergeShare:\ttime:\t%f\t(ms)\n\n",cpu_time_used);

    gettimeofday(&t1, NULL);
    multi_thread_bucket_lipmaa_zeroTest(share_FS,share_BFS);
    gettimeofday(&t2, NULL);
    cpu_time_used = (double)(t2.tv_sec-t1.tv_sec)*1000+(double)(t2.tv_usec-t1.tv_usec)/1000;
    //printf("zeroTest(cpu)\n%d\t%d\t%f\n",uN,oM,cpu_time_used);
    printf("zeroTest:\ttime:\t%f\t(ms)\t",cpu_time_used);
    send_bytes = report_size();
    KB_bytes = send_bytes/1000;
    printf("\tcommu:\t%f\t(KB)\n\n",KB_bytes);
    //printf("zeroTest(comm:)\n%d\t%d\t%zu\n",uN,oM,send_bytes);

    /*
    *   Test the protocol  Again modified-version of Advanced ExtractZ protocol 5
    */
    gettimeofday(&t1, NULL);
    multi_thread_bucket_CExtractZ(share_BFS,share_Z);
    gettimeofday(&t2, NULL);
    cpu_time_used = (double)(t2.tv_sec-t1.tv_sec)*1000+(double)(t2.tv_usec-t1.tv_usec)/1000;
    printf("CExtractZ:\ttime:\t%f\t(ms)\t",cpu_time_used);
    //printf("CExtractZ\n%d\t%d\t%f\n",uN,oM,cpu_time_used);
    send_bytes = report_size();
    KB_bytes = send_bytes/1000;
    printf("\tcommu:\t%f\t(KB)\n\n",KB_bytes);
    //printf("CExtractZ(comm:)\n%d\t%d\t%zu\n",uN,oM,send_bytes);
    batch_mac_check();
    fm_approximateCount(share_Z,oM);
    double unique_countC = fm_approximateCount(share_Z,oM);
    cout<<"approximate:\t"<<unique_countC<<endl;
    

    /*
    *   Test the protocol basic Extract_Z protocol 4
    */
    gettimeofday(&t1, NULL);
    multi_thread_bucket_bExtractZ(share_BFS,share_Z);
    gettimeofday(&t2, NULL);
    cpu_time_used = (double)(t2.tv_sec-t1.tv_sec)*1000+(double)(t2.tv_usec-t1.tv_usec)/1000;
    //printf("bExtractZ\n%d\t%d\t%f\n",uN,oM,cpu_time_used);
    printf("bExtractZ:\ttime:\t%f\t(ms)\t",cpu_time_used);
    send_bytes = report_size();
    KB_bytes = send_bytes/1000;
    printf("\tcommu:\t%f\t(KB)\n\n",KB_bytes);
    //printf("bExtractZ(comm:)\n%d\t%d\t%zu\n",uN,oM,send_bytes);
    batch_mac_check();
    fm_approximateCount(share_Z,oM);
    double unique_countb = fm_approximateCount(share_Z,oM);
    cout<<"approximate:\t"<<unique_countb<<endl;


    cout << "************** End data Aggregation phase **************" << endl;
}

size_t CP::report_size(){
    size_t sent = 0;
    for(int i=0; i<nthreads*communication_multiplier; i++){
        sent += thread_player[i]->sent;
        thread_player[i]->sent = 0;
    }
    return sent;
}

void* CP::thread_Broadcast_single(void* arg){
    thread_sender_info* data = static_cast<thread_sender_info*>(arg);
    CP* obj = data->obj;

    if(data->send_func){
        obj->Broadcast_S_single(*(data->o), data->t_id, data->player_no, true);
    }else{
        obj->Broadcast_R_single(*(data->o), data->t_id, data->player_no, true);
    }
    pthread_exit(NULL);
}


void CP::Broadcast_S_single(octetStream& o, int father_num, int player_no, bool multi_thread){
    TimeScope ts(commu_timer[2*father_num*nparties + player_no]);

    if(player_no > mynum){
        o.Send(CP_sockets[father_num][player_no]);
    }

    if(player_no < mynum){
        o.Send(CP_sockets[nthreads + father_num][player_no]);
    }

    //thread_player[father_num]->sent += o.get_length() * (nparties - 1); it's father thread do this functionality

    if(multi_thread){
        pthread_mutex_lock(&mutex_local_go[father_num]); 
        local_Go[father_num]++;
        pthread_mutex_unlock(&mutex_local_go[father_num]);
    }else{
        thread_player[father_num]->sent += o.get_length();
    }
}

void CP::Broadcast_R_single(octetStream& o, int father_num, int player_no, bool multi_thread){
    TimeScope ts(commu_timer[2*father_num*nparties + nparties + player_no]);
    o.reset_write_head();

    if(player_no < mynum){
        o.Receive(CP_sockets[father_num][player_no]);
    }

    if(player_no > mynum){
        o.Receive(CP_sockets[nthreads + father_num][player_no]);
    }

    if(multi_thread){
        pthread_mutex_lock(&mutex_local_go[father_num]); 
        local_Go[father_num]++;
        pthread_mutex_unlock(&mutex_local_go[father_num]);
    }  
}

void CP::Broadcast_S_and_R(vector<octetStream>& o, int thread_num){
    /*
    *   Implementation 3
    */
    local_Go[thread_num] = 0;

    pthread_t* t = new pthread_t[nparties*2];
    std::vector<thread_sender_info> thread_receive_data(nparties);
    std::vector<thread_sender_info> thread_send_data(nparties);
    for(int i=0; i<nparties; i++){
        if(i == mynum){
            continue;
        }

        thread_receive_data[i].obj = this;
        thread_receive_data[i].t_id = thread_num;
        thread_receive_data[i].o = &(o[i]);
        thread_receive_data[i].player_no = i;
        thread_receive_data[i].send_func = false;
        pthread_create(&t[i], NULL,thread_Broadcast_single, (void*) &thread_receive_data[i]);
    }

    for(int i=0; i<nparties; i++){
        if(i == mynum){
            continue;
        }

        thread_send_data[i].obj = this;
        thread_send_data[i].t_id = thread_num;
        thread_send_data[i].o = &(o[mynum]);
        thread_send_data[i].player_no = i;
        thread_send_data[i].send_func = true;
        pthread_create(&t[i+nparties], NULL,thread_Broadcast_single, (void*) &thread_send_data[i]);
    }



    while(local_Go[thread_num] < (nparties-1)*2){
        usleep(5);
    }

    delete t;
}

void CP::mergeShare(){
    int ofs_M = share_FS.size();
    int ofs_w = share_FS[0].size();
    Share<gfp> sh_tmp;
    sh_tmp.assign_zero();

    for(int i=0; i<ofs_M; i++){
        for(int j=0; j<ofs_w; j++){
            share_FS[i][j]+= sh_tmp;
        }
    }
}

template <class T>
void CP::multi_thread_bucket_CExtractZ(std::vector<std::vector<Share<T>>>& src_fm, Share<T>& Z){
    int ofs_M = src_fm.size();
    if(!ofs_M){
        return;
    }

    std::vector<Share<T>> z(nthreads);
    std::vector<thread_bExtractZ_info<T>> thread_CExtractZ_data(nthreads);

    Go = 0;

    int m_start = 0;
    //int m_internal = ceil(ofs_M/nthreads);
    int m_internal = ceil(((double)ofs_M)/nthreads);

    pthread_t* t = new pthread_t[nthreads];
    for(int i=0; i<nthreads; i++){
        thread_CExtractZ_data[i].obj = this;
        thread_CExtractZ_data[i].src_fm = &src_fm;
        thread_CExtractZ_data[i].z = &z[i];
        thread_CExtractZ_data[i].t_id = i;
        thread_CExtractZ_data[i].m_start = m_start;
        m_start+=m_internal;
        thread_CExtractZ_data[i].m_internal = m_internal;
        pthread_create(&t[i], NULL,thread_CExtractZ, (void*) &thread_CExtractZ_data[i]);
    }

    while(Go < nthreads){
        usleep(10);
    }

    Z.assign_zero();
    for(int i=0; i<nthreads; i++){
        Z+=z[i];
    }

    delete t;
}


template <class T>
void CP::bucket_CExtractZ_basic(const std::vector<std::vector<Share<T>>>& src_fm, Share<T>& z, const int m_start, int m_internal, int thread_num, bool multi_thread){
    /*
    *   PreProcessing data and dataStructure
    */
    int ofs_M = src_fm.size();
    int ofs_w = src_fm[0].size();
    if(!ofs_M){
        return;
    }
    m_internal = ofs_M < (m_start+m_internal) ? ofs_M - m_start : m_internal;
    z.assign_zero();
    int end = ofs_w;
    int t = ceil(ofs_w/2.0);
    int tmp_int = 2*t;
    std::vector<std::vector<Share<T>>> fm_sketch(m_internal);
    for(int m=0; m<m_internal; m++){
        fm_sketch[m].resize(tmp_int);
    }

    // temp variable
    //std::vector<T> a_array(2*t),mac_array(2*t);
    int offset;
    vector<Share<T>> sh_x(m_internal),sh_m(m_internal),R0(m_internal);
    Share<T> sh_tmp;
    T mac,tmp;
    std::vector<T> m(m_internal);

    // bit-wise negate
    Share<T> sh_one(1,mynum,alphai);
    vector<Share<T>> B0(m_internal);
    T score;
    tmp_int = m_start;
    for(int k=0; k<m_internal; k++){
        for(int i=0; i<ofs_w; i++){
            fm_sketch[k][i].sub(sh_one,src_fm[tmp_int][i]);
        }
        for(int i=ofs_w; i<2*t; i++){
            fm_sketch[k][i].assign(sh_one);
        }
        tmp_int++;
    }
    
    vector<vector<Share<T>>> R_group(m_internal);
    for(int k=0; k<m_internal; k++){
        R_group[k].resize(t+1);
        // getR_group(R_group[k],t+1,thread_num);
        // R0[k] = R_group[k][0];
    }

    /*
    *   Start Protocol
    */
    while(end>3){
        const vector<T>& poly_f = get_CExtractZ_PolyNomailFunc(t);  
        for(int k=0; k<m_internal; k++){
            getR_group(R_group[k],t+1,thread_num);
            R0[k] = R_group[k][0];

            sh_x[k].assign_zero();
            for(int i=0; i<t; i++){
                sh_x[k]+=fm_sketch[k][i];
            }
        }
        
        dot_product(sh_x, R0, sh_m,thread_num);  
        delay_open_and_check(sh_m,m,thread_num);

        for(int k=0; k<m_internal; k++){
            B0[k].assign(sh_one);
            B0[k] += sh_x[k]*poly_f[1];
        }

        
        // get the x, x^2,..., x^d
        for(int k=0; k<m_internal; k++){
            T& tmp_m = m[k];
            tmp.assign(tmp_m);
            tmp = tmp*tmp_m;
            for(int i=2; i<t+1; i++){
                sh_tmp = R_group[k][i]*tmp;
                B0[k] += sh_tmp*poly_f[i];
                tmp = tmp*tmp_m;
            } 
        }
        
        score.assign(t);
        for(int k=0; k<m_internal; k++){
            z+=score*B0[k];
        }
        
        _dot_product(B0,fm_sketch,fm_sketch,t,t,t);

        for(int k=0; k<m_internal; k++){
            //_dot_product(B0[k],fm_sketch[k],fm_sketch[k],t,t,t);
            offset = t;
            for(int i=0; i<t; i++){
                fm_sketch[k][i]+=fm_sketch[k][offset++];
            }
        }
        
        end = t;
        t = ceil(end/2.0);

        for(int k=0; k<m_internal; k++){
            for(int i=end; i<2*t; i++){
                fm_sketch[k][i].assign(sh_one);
            }
        }    
    }

    //lookup_x_CExtractZ(fm_sketch,sh_tmp,0,end);
    int len_group = pow(2,end);
    for(int k=0; k<m_internal; k++){
        getR_group(R_group[k],len_group,thread_num);
        R0[k] = R_group[k][0];
        sh_x[k].assign_zero();
        
    }

    tmp_int = 1;
    for(int i=end-1; i>=0; i--){
        for(int k=0; k<m_internal; k++){
            sh_x[k] += fm_sketch[k][i]*tmp_int;
        }
        tmp_int *= 2;
    }


    const vector<T>& poly_f = get_CExtractZ_PolyNomailFunc(len_group);
    // multiply x*R^-1 then reveal it
    dot_product(sh_x, R0, sh_m,thread_num);  
    delay_open_and_check(sh_m,m);

    Share<T> sh_end(end,mynum,alphai);

    for(int k=0; k<m_internal; k++){
        z += sh_end;
        z += sh_x[k]*poly_f[1];

        tmp.assign(m[k]);
        tmp = tmp*m[k];

        for(int i=2; i<len_group; i++){
            sh_tmp = R_group[k][i]*tmp;
            z += sh_tmp*poly_f[i];
            tmp = tmp*m[k];
        }
    }
    
 
    if(multi_thread){
        pthread_mutex_lock(&mutex_go); 
        Go++;
        pthread_mutex_unlock(&mutex_go);
    }
}
template <class T>
void CP::CExtractZ(const std::vector<Share<T>>& src_fm, Share<T>& z, int thread_num ){
    /*
    *   PreProcessing data and dataStructure
    */
    z.assign_zero();
    int ofs_w = src_fm.size();
    int end = ofs_w;
    int t = ceil(ofs_w/2.0);
    std::vector<Share<T>> fm_sketch(2*t);

    // bit-wise negate
    Share<T> sh_one(1,mynum,alphai);
    Share<T> B0;
    T score;
    for(int i=0; i<ofs_w; i++){
        fm_sketch[i].sub(sh_one,src_fm[i]);
    }
    for(int i=ofs_w; i<2*t; i++){
        fm_sketch[i].assign(sh_one);
    }

    std::vector<Share<T>> R_group(t+1);

    // temp variable
    //std::vector<T> a_array(2*t),mac_array(2*t);
    int offset;
    Share<T> sh_x,sh_m,sh_tmp;
    T m,mac,tmp;

    /*
    *   Start Protocol
    */
    while(end>3){
        const vector<T>& poly_f = get_CExtractZ_PolyNomailFunc(t);  
        getR_group(R_group,t+1,thread_num);

        sh_x.assign_zero();
        for(int i=0; i<t; i++){
            sh_x+=fm_sketch[i];
        }

        dot_product(sh_x, R_group[0], sh_m,thread_num);  
        delay_open_and_check(sh_m,m,thread_num);

        B0.assign(sh_one);
        B0 += sh_x*poly_f[1];
        tmp.assign(m);
        tmp = tmp*m;
        // get the x, x^2,..., x^d
       
        for(int i=2; i<t+1; i++){
            sh_tmp = R_group[i]*tmp;
            B0 += sh_tmp*poly_f[i];
            tmp = tmp*m;
        } 


        //bdb_CExtractZ(fm_sketch,B0,t);
        score.assign(t);
        z+=score*B0;

        _dot_product(B0,fm_sketch,fm_sketch,t,t,t);

        offset = t;
        for(int i=0; i<t; i++){
            fm_sketch[i]+=fm_sketch[offset++];
        }
        end = t;
        t = ceil(end/2.0);
        for(int i=end; i<2*t; i++){
            fm_sketch[i].assign(sh_one);
        }
    }

    //lookup_x_CExtractZ(fm_sketch,sh_tmp,0,end);
    int len_group = pow(2,end);
    getR_group(R_group,len_group,thread_num);
    const vector<T>& poly_f = get_CExtractZ_PolyNomailFunc(len_group);

    int tmp_int = 1;
    sh_x.assign_zero();

    for(int i=end-1; i>=0; i--){
        sh_x += fm_sketch[i]*tmp_int;
        tmp_int *= 2;
    } 

    // multiply x*R^-1 then reveal it
    dot_product(sh_x, R_group[0], sh_m,thread_num);  
    delay_open_and_check(sh_m,m);

    Share<T> sh_end(end,mynum,alphai);
    B0.assign(sh_end);
    B0 += sh_x*poly_f[1];
    tmp.assign(m);
    tmp = tmp*m;

    for(int i=2; i<len_group; i++){
        sh_tmp = R_group[i]*tmp;
        B0 += sh_tmp*poly_f[i];
        tmp = tmp*m;
    }
 
    z+=B0;
}


template <class T>
void CP::multi_thread_bucket_bExtractZ(std::vector<std::vector<Share<T>>>& src_fm, Share<T>& Z){
    if(thread_player.size() == 0){
        throw runtime_error("The p2p socket has been expired");
    }

    int ofs_M = src_fm.size();
    if(!ofs_M){
        return;
    }

    std::vector<Share<T>> z(nthreads);
    std::vector<thread_bExtractZ_info<T>> thread_bExtractZ_data(nthreads);

    Go = 0;

    int m_start = 0;
    int m_internal = ceil(((double)ofs_M)/nthreads);

    pthread_t* t = new pthread_t[nthreads];
    for(int i=0; i<nthreads; i++){
        thread_bExtractZ_data[i].obj = this;
        thread_bExtractZ_data[i].src_fm = &src_fm;
        thread_bExtractZ_data[i].z = &z[i];
        thread_bExtractZ_data[i].t_id = i;
        thread_bExtractZ_data[i].m_start = m_start;
        m_start+=m_internal;
        thread_bExtractZ_data[i].m_internal = m_internal;
        pthread_create(&t[i], NULL,thread_bExtractZ, (void*) &thread_bExtractZ_data[i]);
    }

    while(Go < nthreads){
        usleep(10);
    }

    Z.assign_zero();
    for(int i=0; i<nthreads; i++){
        Z+=z[i];
    }

    delete t;
}



void* CP::thread_bExtractZ(void* arg){
    thread_bExtractZ_info<gfp>* data = static_cast<thread_bExtractZ_info<gfp>*>(arg);
    CP* obj = data->obj;
    obj->bucket_bExtractZ(*(data->src_fm),*(data->z),data->m_start,data->m_internal,data->t_id,true);

    pthread_exit(NULL);
}


template <class T>
double CP::fm_approximateCount(const Share<T>& Z, int M){
    const double PHI = 0.77351;
    T a,mac;
    bigint tmp;
    open_and_check(Z,a,mac);
    to_bigint(tmp,a);
    double R = tmp.get_si();
    R /= M;
    return pow(2.0, R)/PHI;

}

void* CP::thread_lipmaa_zeroTest(void* arg){
    thread_lipmaaZero_info<gfp>* data = static_cast<thread_lipmaaZero_info<gfp>*>(arg);
    CP* obj = data->obj;

    obj->bucket_lipmaa_zeroTest_bench(*(data->src_vec),*(data->t0),data->m_start,data->m_internal,data->t_id,true);

    pthread_exit(NULL);
}
template <class T>
void CP::open_and_check(const Share<T>& share_value, T& a, T& mac){
    if(player == 0){
        throw runtime_error("The p2p socket has been expired");
    }
    vector<Share<T>> res_share(nparties);
    vector<octetStream> vec_shares(nparties);
    share_value.pack(vec_shares[mynum]);
    player->Broadcast_Receive(vec_shares);

    for(int k=0;k< nparties;k++){
        res_share[k].unpack(vec_shares[k]);
    }
    check_share(res_share,a,mac,nparties,keyp);
}

template <class T>
void CP::delay_open_and_check(const Share<T>& share_value, T& a, int thread_num){
    if(thread_player[thread_num] == 0){
        throw runtime_error("The p2p socket has been expired");
    }
    vector<Share<T>> res_share(nparties);
    vector<octetStream> vec_shares(nparties);
    share_value.pack(vec_shares[mynum]);
    thread_player[thread_num]->Broadcast_Receive(vec_shares);

    std::vector<PSUCA::W2C_mac>& wait_queue = W2C_mac_queue[thread_num];
    int& count = n_wait2Check[thread_num];

    for(int k=0;k< nparties;k++){
        res_share[k].unpack(vec_shares[k]);
        wait_queue[count].value.add(res_share[k].get_share());
        wait_queue[count].mac.add(res_share[k].get_mac());
    }

    a.assign(wait_queue[count].value);
    count++;
    if(count == PSUCA::max_w2c){
        batch_mac_check(thread_num);
    }
}

void CP::init_wait_to_check_buffer(){
    n_wait2Check.resize(nthreads*communication_multiplier);
    W2C_mac_queue.resize(nthreads*communication_multiplier);

    for(int i=0; i<nthreads*communication_multiplier; i++){
        n_wait2Check[i] = 0;
        W2C_mac_queue[i].resize(PSUCA::max_w2c+10);
    }
}


template <class T>
void CP::delay_open_and_check(const std::vector<Share<T>>& share_value, std::vector<T>& a, int thread_num){
    if(thread_player[thread_num] == 0){
        throw runtime_error("The p2p socket has been expired");
    }

    int size_array = share_value.size();
    if(!size_array){
        return;
    }

    std::vector<PSUCA::W2C_mac>& wait_queue = W2C_mac_queue[thread_num];
    int& count = n_wait2Check[thread_num];

    if((size_array+count) >= PSUCA::max_w2c){
        batch_mac_check(thread_num);
    }

    vector<Share<T>> res_share(nparties);
    vector<octetStream> vec_shares(nparties);

    for(int i=0; i<size_array; i++){
        share_value[i].pack(vec_shares[mynum]);
    } 

    
    if(local_envs){
        thread_player[thread_num]->Broadcast_Receive(vec_shares);
    }else{
        thread_player[thread_num]->sent += vec_shares[mynum].get_length() * (nparties - 1);
        Broadcast_S_and_R(vec_shares, thread_num);
    }
    

    for(int i=0; i<size_array; i++){
        for(int k=0; k<nparties; k++){
            res_share[k].unpack(vec_shares[k]);
            wait_queue[count].value.add(res_share[k].get_share());
            wait_queue[count].mac.add(res_share[k].get_mac());
        }
        a[i].assign(wait_queue[count].value);
        count ++;
    }
}

template <class T>
void CP::_delay_open_and_check(const std::vector<Share<T>>& share_value, std::vector<T>& a, int n_elements, int start_share, int start_a, int thread_num){
    if(thread_player[thread_num] == 0){
        throw runtime_error("The p2p socket has been expired");
    }

    std::vector<PSUCA::W2C_mac>& wait_queue = W2C_mac_queue[thread_num];
    int& count = n_wait2Check[thread_num];

    if((n_elements+count) >= PSUCA::max_w2c){
        batch_mac_check(thread_num);
    }

    vector<Share<T>> res_share(nparties);
    vector<octetStream> vec_shares(nparties);

    for(int i=start_share; i<start_share+n_elements; i++){
        share_value[i].pack(vec_shares[mynum]);
    } 

    if(local_envs){
        thread_player[thread_num]->Broadcast_Receive(vec_shares);
    }else{
        thread_player[thread_num]->sent += vec_shares[mynum].get_length() * (nparties - 1);
        Broadcast_S_and_R(vec_shares, thread_num);
    }
    
    for(int i=start_a; i<start_a+n_elements; i++){
        for(int k=0; k<nparties; k++){
            res_share[k].unpack(vec_shares[k]);
            wait_queue[count].value.add(res_share[k].get_share());
            wait_queue[count].mac.add(res_share[k].get_mac());
        }
        a[i].assign(wait_queue[count].value);
        count ++;
    }
}
void CP::batch_mac_check(int thread_num){
    gfp res;
    if(thread_num == -1){
        for(int k=0; k<nthreads*communication_multiplier; k++){
            for(int i=0; i<n_wait2Check[k]; i++){
                res.mul(W2C_mac_queue[k][i].value,keyp);
                if (!res.equal(W2C_mac_queue[k][i].mac))
                    {
                      cout << "Value:      " << W2C_mac_queue[k][i].value << endl;
                      cout << "Input MAC:  " << W2C_mac_queue[k][i].mac << endl;
                      cout << "Actual MAC: " << res << endl;
                      cout << "MAC key:    " << keyp << endl;
                      throw mac_fail();
                    }
                W2C_mac_queue[k][i].value.assign(0);
                W2C_mac_queue[k][i].mac.assign(0);
            }
            n_wait2Check[k] = 0;
        } 
    }else{
        std::vector<PSUCA::W2C_mac>& wait_queue = W2C_mac_queue[thread_num];
        int& count = n_wait2Check[thread_num];

        for(int i=0; i<count; i++){
            res.mul(wait_queue[i].value,keyp);
            if (!res.equal(wait_queue[i].mac))
                {
                  cout << "Value:      " << wait_queue[i].value << endl;
                  cout << "Input MAC:  " << wait_queue[i].mac << endl;
                  cout << "Actual MAC: " << res << endl;
                  cout << "MAC key:    " << keyp << endl;
                  throw mac_fail();
                }
            wait_queue[i].value.assign(0);
            wait_queue[i].mac.assign(0);
        }
        count = 0;
    }
    
}
template <class T>
void CP::thread_open_and_check(const std::vector<Share<T>>& share_value, std::vector<T>& a, std::vector<T>& mac, int thread_num){
    if(thread_player[thread_num] == 0){
        throw runtime_error("The p2p socket has been expired");
    }

    int size_array = share_value.size();
    if(!size_array){
        return;
    }

    vector<Share<T>> res_share(nparties);
    vector<octetStream> vec_shares(nparties);

    for(int i=0; i<size_array; i++){
        share_value[i].pack(vec_shares[mynum]);
    } 

    thread_player[thread_num]->Broadcast_Receive(vec_shares);
    
    for(int i=0; i<size_array; i++){
        for(int k=0; k<nparties; k++){
            res_share[k].unpack(vec_shares[k]);
        }
        check_share(res_share,a[i],mac[i],nparties,keyp);
    }
}

template <class T>
void CP::open_and_check(const vector<Share<T>>& share_value, vector<T>& a, std::vector<T>& mac){
    if(player == 0){
        throw runtime_error("The p2p socket has been expired");
    }

    int size_array = share_value.size();
    if(!size_array){
        return;
    }

    vector<Share<T>> res_share(nparties);
    vector<octetStream> vec_shares(nparties);

    for(int i=0; i<size_array; i++){
        share_value[i].pack(vec_shares[mynum]);
    } 
    player->Broadcast_Receive(vec_shares);
    for(int i=0; i<size_array; i++){
        for(int k=0; k<nparties; k++){
            res_share[k].unpack(vec_shares[k]);
        }
        check_share(res_share,a[i],mac[i],nparties,keyp);
    }
}



void CP::start(){
    benchmark();
}

DP_connector::DP_connector(int _mynum,int pnb,int my_port, int _num_dp_players, int _oM, int _uN, int _nthreads){
    oM = _oM;
    uN = _uN;
    dataF = 0;
    player = 0;
    nthreads = _nthreads;
    init(_mynum,pnb,my_port,_num_dp_players);
    sent = 0;
}
void DP_connector::init(int _mynum,int pnb,int _my_port, int _num_dp_players)
{  
    nplayers = _num_dp_players;
    mynum=_mynum;
    portnum_base=pnb;
    if(_my_port == DP_connector::DEFAULT_PORT){
      my_port = portnum_base+DP_connector::OFFSET_PORT+mynum;
    }
    setup_server();
}

void DP_connector::key_init(gfp& _alphai,gfp& _keyp){
    alphai.assign(_alphai);
    keyp.assign(_keyp);
}
void DP_connector::player_init(Player* _player){
    player = _player;
}

DP_connector::~DP_connector()
{
    if (server != 0)
        delete server;
}



void DP_connector::setup_server()
{
  server = new ServerSocket(my_port);
  server->init();
}

void DP_connector::share_zero_array(int socket_num, std::vector<Share<gfp>>& share_Z){
    std::vector<std::vector<Share<gfp>>> tmp_gfp_array(1);
    tmp_gfp_array[0].resize(oM);
    share_Z.clear();
    share_Z.resize(oM);
    share_value(tmp_gfp_array,socket_num);
    for(int i=0; i<oM; i++){
        share_Z[i] = tmp_gfp_array[0][i];
    }

}
bool DP_connector::data_collection(int socket_num, std::vector<std::vector<Share<gfp>>>& share_FS){
    bool ret = true;
    int w = ceil(log2 (uN))+4;

    share_FS.clear();

    share_FS.resize(oM);
    for(int m=0; m<oM; m++){
        share_FS[m].resize(w);
    }

    cout<<oM<<"  "<<w<<endl;
    ret = share_value(share_FS,socket_num);

    //ret = check_share_value(tmp_gfp_array,2);

    return ret;    
    
}
void DP_connector::start(Data_Files* _df, std::vector<std::vector<Share<gfp>>>& share_FS){
    dataF = _df;
    /*
    * Collect the data of each data party orderly
    */
    int i=0;
    int socket_num = -1;
    if (server == 0)
        throw runtime_error("The socket communication failure detection");


    /***********************************Start For benchmark*********************************/
    cout << "************** Start data Collection phase **************" << endl;
    //time measure tools 
    struct timeval t1;
    struct timeval t2;
    double cpu_time_used;
    double KB_bytes;
    
    for (i=1; i<=nplayers; i++)
    {
        cerr << "Waiting for Data Party " << i << endl;
        socket_num = server->get_connection_socket(i);
        cerr << "Connected to Data Party " << i << endl;
        
        send(socket_num, GO);

        Share<gfp> test_x; 
        test_x.assign_zero();

        share_value(test_x,socket_num);

        gettimeofday(&t1, NULL);
        data_collection(socket_num, share_FS);
        gettimeofday(&t2, NULL);
        cpu_time_used = (double)(t2.tv_sec-t1.tv_sec)*1000+(double)(t2.tv_usec-t1.tv_usec)/1000;
        printf("2 times Share(x):\ttime:\t%f\t(ms)\t",cpu_time_used);
        
    
        size_t send_bytes = report_size();
        KB_bytes = send_bytes*2/1000;
        printf("\tcommu:\t%f\t(KB)\n\n",KB_bytes);

        close(socket_num);
        socket_num = -1;
    }

    cout << "************** End Data Collection phase **************" << endl;
}

void DP_connector::send_to(int data_player,const octetStream& o,bool donthash) const
{
    TimeScope ts(timer);
    o.Send(data_player);
    if (!donthash)
        { blk_SHA1_Update(&ctx,o.get_data(),o.get_length()); }
    sent += o.get_length();
}

void DP_connector::receive_player(int data_player,octetStream& o,bool donthash) const{
    TimeScope ts(timer);
    o.reset_write_head();
    o.Receive(data_player);
    sent += o.get_length();
    if (!donthash)
        { blk_SHA1_Update(&ctx,o.get_data(),o.get_length()); }
}


template <class T>
bool DP_connector::share_value(Share<T>& x, int socket_num){
    if(socket_num == -1){
        throw runtime_error("The socket has been expired");
    }

    
    Share<T> a;
    T a_value,a_mac,y_value,tmp; //share value of a
    octetStream share_stream;

    if(!dataF->eof<T>(DATA_RAND)){
        dataF->get_one(DATA_MODP, DATA_RAND, a); //this should be changed as rand afterly
    }else{
        throw runtime_error("Cannot read the random share from the data files");
    }


    a_value = a.get_share();
    a_value.pack(share_stream);
    //  reveal a to data party
    send_to(socket_num,share_stream);
    share_stream.reset_write_head();
    //  receive x-a from data party
    receive_player(socket_num,share_stream);
    y_value.unpack(share_stream);

    //compute the data party to-share value x
    a_mac = a.get_mac();
    tmp.assign_zero();
    if(mynum == 0){
        a_value.add(y_value);
        tmp.mul(y_value,alphai);
        a_mac.add(tmp);
    }else{
        tmp.mul(y_value,alphai);
        a_mac.add(tmp);
    }

    x.set_share(a_value);
    x.set_mac(a_mac);

    return true;
}
size_t DP_connector::report_size(){
    size_t result = sent;
    sent = 0;
    return result;
}
template <class T>
bool DP_connector::share_value(std::vector<std::vector<Share<T>>>& des_vec, int socket_num){
    if(socket_num == -1){
        throw runtime_error("The socket has been expired");
    }

    int ofs_M = des_vec.size();
    int ofs_w = des_vec[0].size();
    int elements_size = ofs_M*ofs_w;

    std::vector<Share<T>> a_Array(elements_size);
    T tmp, a_mac, y_value, a_value;

    octetStream share_stream;

    int j = 0;
    while ((j<elements_size) && (!dataF->eof<T>(DATA_RAND)))
    {
        dataF->get_one(DATA_MODP, DATA_RAND, a_Array[j]); //this should be changed as rand afterly
        j++;
    }

    if(j < elements_size)
        throw runtime_error("cannot read enough random share from offline file");

    for(int i=0; i<elements_size; i++){
        a_Array[i].get_share().pack(share_stream);
    }

    //  reveal a to data party
    send_to(socket_num,share_stream);
    share_stream.reset_write_head();
    //  receive x-a from data party
    receive_player(socket_num,share_stream);

    for(int m=0; m<ofs_M; m++){
        for(int i=0; i<ofs_w; i++){
            y_value.unpack(share_stream);
            //compute the data party to-share value x
            a_mac = a_Array[m*ofs_w+i].get_mac();
            a_value = a_Array[m*ofs_w+i].get_share();
            tmp.assign_zero();

            if(mynum == 0){
                a_value.add(y_value);
                tmp.mul(y_value,alphai);
                a_mac.add(tmp);
            }else{
                tmp.mul(y_value,alphai);
                a_mac.add(tmp);
            }

            des_vec[ofs_M-m-1][ofs_w-i-1].set_share(a_value);
            des_vec[ofs_M-m-1][ofs_w-i-1].set_mac(a_mac);
        }
    }
    return true;
}