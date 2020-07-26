// (C) 2018 University of Bristol. See License.txt

/*
 * PairwiseGenerator.cpp
 *
 */

#include "FHEOffline/PairwiseGenerator.h"
#include "FHEOffline/PairwiseMachine.h"
#include "FHEOffline/Producer.h"
#include "Auth/Subroutines.h"

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/GF2EXFactoring.h>
#include <NTL/GF2XFactoring.h>
NTL_CLIENT

#define amortize 8192

template <class FD>
PairwiseGenerator<FD>::PairwiseGenerator(int thread_num,
        PairwiseMachine& machine,
        Dtype data_type) :
    GeneratorBase(thread_num, machine.N),
    data_type(data_type),
    producer(machine.setup<FD>().FieldD, machine.N.my_num(),
            thread_num, (machine.output && data_type == DATA_TRIPLE), machine.setup<FD>().dirname),
    EC(P, machine.other_pks, machine.setup<FD>().FieldD, timers, machine, *this),
    
    // Our extended part
    randProducer(machine.setup<FD>().FieldD, machine.N.my_num(),
            thread_num, (machine.output && data_type == DATA_RAND), machine.setup<FD>().dirname),

    C(machine.sec, machine.setup<FD>().params), volatile_memory(0),
    machine(machine), global_player(machine.N, (1LL << 28) + (thread_num << 16))
{
    for (int i = 1; i < machine.N.num_players(); i++)
        multipliers.push_back(new Multiplier<FD>(i, *this));
    const FD& FieldD = machine.setup<FD>().FieldD;
    a.resize(machine.sec, FieldD);
    b.resize(machine.sec, FieldD);
    c.resize(machine.sec, FieldD);
    a.allocate_slots(FieldD.get_prime());
    b.allocate_slots(FieldD.get_prime());
    // extra limb for addition
    c.allocate_slots((bigint)FieldD.get_prime() << 64);
    b_mod_q.resize(machine.sec,
    { machine.setup<FD>().params, evaluation, evaluation });

    // Our extended part
    r.resize(machine.sec, FieldD);
    r.allocate_slots(FieldD.get_prime());
}

template <class FD>
PairwiseGenerator<FD>::~PairwiseGenerator()
{
    for (auto m : multipliers)
        delete m;
}

template <class FD>
void PairwiseGenerator<FD>::run()
{
    PRNG G;
    G.ReSeed();
    MAC_Check<typename FD::T> MC(machine.setup<FD>().alphai);

    stringstream ss;
    ifstream triple_pool;
    ifstream rand_pool;
    ifstream rand2_pool;

    ofstream outf;

    int left_todo;
    bool playerone = (P.my_num() == 0)? true : false;
    typename FD::T alphai = MC.get_alphai();


    vector<Share<typename FD::T> > Sh_merge(2*amortize);
    vector<typename FD::T> merge(2*amortize);
    vector<Share<typename FD::T> > Sh_batch_a(amortize);
    vector<Share<typename FD::T> > Sh_batch_b(amortize);
    vector<Share<typename FD::T> > Sh_batch_c(amortize);
    vector<Share<typename FD::T> > Sh_temp(amortize);
    vector<typename FD::T> temp(amortize);
    vector<Share<typename FD::T> > Sh_res(amortize*(machine.exp + 1));

    vector<Share<typename FD::T> > Sh_batch_r(amortize);
    vector<Share<typename FD::T> > Sh_batch_r2(amortize);
    vector<typename FD::T> batch_r2(amortize);
    typename FD::T root;
    typename FD::T one(1);
    typename FD::T inv_two(2); inv_two.invert();
    stringstream ss_p;
    ZZ p;
    ss_p << gfp::pr();
    ss_p >> p;

    ZZ z;
    ZZ res;

    switch (data_type)
    {
    case DATA_TRIPLE:
        cout << "[*] Low Gear triple generation ..." << endl;
        while (total < machine.nTriplesPerThread)
        {
            timers["Randomization"].start();
            a.randomize(G);
            b.randomize(G);
            timers["Randomization"].stop();

            size_t prover_memory = EC.generate_proof(C, a, ciphertexts, cleartexts);

            timers["Plaintext multiplication"].start();
            c.mul(a, b);
            timers["Plaintext multiplication"].stop();

            timers["FFT of b"].start();
            for (int i = 0; i < machine.sec; i++)
                b_mod_q.at(i).from_vec(b.at(i).get_poly());
            timers["FFT of b"].stop();

            timers["Proof exchange"].start();
            size_t verifier_memory = EC.create_more(ciphertexts, cleartexts);
            timers["Proof exchange"].stop();

            volatile_memory = max(prover_memory, verifier_memory);

            Rq_Element values({machine.setup<FD>().params, evaluation, evaluation});
            for (int k = 0; k < machine.sec; k++)
            {
                producer.ai = a[k];
                producer.bi = b[k];
                producer.ci = c[k];

                for (int j = 0; j < 3; j++)
                {
                    timers["Plaintext multiplication"].start();
                    producer.macs[j].mul(machine.setup<FD>().alpha, producer.values[j]);
                    timers["Plaintext multiplication"].stop();

                    if (j == 1)  // represent b in FFT form
                        values = b_mod_q[k];
                    else
                    {
                        timers["Plaintext conversion"].start();
                        values.from_vec(producer.values[j].get_poly());
                        timers["Plaintext conversion"].stop();
                    }

                    for (auto m : multipliers)
                        m->multiply_alpha_and_add(producer.macs[j], values);
                }
                producer.reset();
                total += producer.sacrifice(P, MC);
                if (total >= machine.nTriplesPerThread)
                    break;
            }

            timers["MAC check"].start();
            MC.Check(P);
            timers["MAC check"].stop();
        }
        timers.insert(producer.timers.begin(), producer.timers.end());
        break;
    case DATA_RAND:
        cout << "[*] Low Gear Rand() generation ..." << endl;
        while (total < machine.nRandsPerThread)
        {
            timers["Randomization"].start();
            r.randomize(G);
            timers["Randomization"].stop();

            for (int k = 0; k < machine.sec; ++k)
            {
                randProducer.ri = r[k];

                timers["Plaintext multiplication"].start();
                randProducer.mac.mul(machine.setup<FD>().alpha, randProducer.value);
                timers["Plaintext multiplication"].stop();

                // Call F_rand to sample a set of random coefficients
                vector<typename FD::T> rand_coeffs;
                // rand_coeffs.resize(randProducer.num_slots());
                // for (unsigned int i = 0; i < rand_coeffs.size(); ++i)
                // {
                //     Create_Random(rand_coeffs[i], P);
                // }
                Create_Random_Many(rand_coeffs, P, randProducer.num_slots());

                for (auto m : multipliers)
                    m->multiply_alpha_and_add(randProducer.mac, randProducer.value,
                        rand_coeffs,
                        machine.setup<FD>().alphai);

                randProducer.reset();
                timers["Rand generation"].start();
                total += randProducer.sacrifice(P, MC);
                timers["Rand generation"].stop();
                if (total >= machine.nRandsPerThread)
                    break;
            }
        }
        timers.insert(randProducer.timers.begin(), randProducer.timers.end());
        break;
    case DATA_RAND2:
        cout << "[*] Low Gear Rand2() generation ..." << endl;

        ss.str("");
        ss << PREP_DIR << "/" << P.num_players() << "-" << numBits(gfp::pr()) << "-" << machine.sec
            << "/Rands-p-P" << P.my_num();
        if (thread_num)
            ss << "-" << thread_num;
        rand_pool.open(ss.str().c_str(), ios::in | ios::binary);

        ss.str("");
        ss << PREP_DIR << "/" << P.num_players() << "-" << numBits(gfp::pr()) << "-" << machine.sec
            << "/Triples-p-P" << P.my_num();
        if (thread_num)
            ss << "-" << thread_num;
        triple_pool.open(ss.str().c_str(), ios::in | ios::binary);

        if (rand_pool.fail() || triple_pool.fail()) { throw file_error(); }

        if (machine.output) {
            ss.str("");
            ss << PREP_DIR << "/" << P.num_players() << "-" << numBits(gfp::pr()) << "-" << machine.sec
                << "/Rand2s-p-P" << P.my_num();
            if (thread_num)
                ss << "-" << thread_num;
            outf.open(ss.str().c_str(), ios::out | ios::binary);
        }

        left_todo = machine.nRand2sPerThread;

        timers["Generating Rand2()"].start();
        while (left_todo > 0) {
        #ifdef BENCHMARK
            rand_pool.clear();
            rand_pool.seekg(ios::beg);
            triple_pool.clear();
            triple_pool.seekg(ios::beg);
        #endif

            int this_loop = amortize;
            if (this_loop > left_todo) {
                this_loop = left_todo;
                Sh_batch_r.resize(this_loop);
                Sh_batch_a.resize(this_loop);
                Sh_batch_b.resize(this_loop);
                Sh_batch_c.resize(this_loop);

                Sh_merge.resize(2*this_loop);
                merge.resize(2*this_loop);
                Sh_batch_r2.resize(this_loop);
                batch_r2.resize(this_loop);
            }

            for (int i = 0; i < this_loop; i++) {
                Sh_batch_r[i].input(rand_pool, READABLE);
                Sh_batch_a[i].input(triple_pool, READABLE);
                Sh_batch_b[i].input(triple_pool, READABLE);
                Sh_batch_c[i].input(triple_pool, READABLE);

                Sh_merge[i] = Sh_batch_r[i];
                Sh_merge[i].sub(Sh_merge[i], Sh_batch_a[i]);
                Sh_merge[i + this_loop] = Sh_batch_r[i];
                Sh_merge[i + this_loop].sub(Sh_merge[i + this_loop], Sh_batch_b[i]);
            }
            
            MC.POpen_Begin(merge, Sh_merge, P);
            MC.POpen_End(merge, Sh_merge, P);
            
            for (int i = 0; i < this_loop; i++) {
                Sh_batch_r2[i] = Sh_batch_c[i] + merge[i]*Sh_batch_b[i] + merge[i + this_loop]*Sh_batch_a[i];
                Sh_batch_r2[i].add(Sh_batch_r2[i], merge[i]*merge[i + this_loop], playerone, alphai);
            }

            MC.POpen_Begin(batch_r2, Sh_batch_r2, P);
            MC.POpen_End(batch_r2, Sh_batch_r2, P);

            if (machine.output) {
                for (int i = 0; i < this_loop; i++) {
                    if (batch_r2[i].is_zero()) {
                        throw Offline_Check_Error("Meet a bad zero a^2");
                        //continue;
                    }
                    stringstream ss1;
                    ss1 << batch_r2[i];
                    ss1 >> z;
                    z = z%p;

                    //root = batch_r2[i].sqrRoot();     // SPDZ ver., too slow
                    res = SqrRootMod(z, p);             // NTL ver.

                    ZZ temp = p;
                    if (res <= 0) { res = abs(res); }
                    else if (res >= (temp >>= 1)) {
                        res = p - res;
                    }

                    stringstream ss2;
                    ss2 << res;
                    ss2 >> root;

                    root.invert();
                    Sh_batch_r[i] *= root;
                    Sh_batch_r[i].add(Sh_batch_r[i], one, playerone, alphai);
                    Sh_batch_r[i] *= inv_two;

                    Sh_batch_r[i].output(outf, READABLE);
                }
            }
            left_todo -= this_loop;
        }
        timers["Generating Rand2()"].stop();
        if (machine.output)
            outf.close();

        rand_pool.close();
        triple_pool.close();

        total += machine.nRand2sPerThread;
        break;
    case DATA_RANDEXP:
        cout << "[*] Low Gear RandExp(l) generation ..." << endl;
        
        ss.str("");
        ss << PREP_DIR << "/" << P.num_players() << "-" << numBits(gfp::pr()) << "-" << machine.sec
            << "/Rands-p-P" << P.my_num();
        if (thread_num)
        	ss << "-" << thread_num;
        rand_pool.open(ss.str().c_str(), ios::in | ios::binary);

        ss.str("");
        ss << PREP_DIR << "/" << P.num_players() << "-" << numBits(gfp::pr()) << "-" << machine.sec
            << "/Triples-p-P" << P.my_num();
        if (thread_num)
        	ss << "-" << thread_num;
        triple_pool.open(ss.str().c_str(), ios::in | ios::binary);

        if (rand_pool.fail() || triple_pool.fail()) { throw file_error(); }

        if (machine.output) {
            ss.str("");
            ss << PREP_DIR << "/" << P.num_players() << "-" << numBits(gfp::pr()) << "-" << machine.sec
                << "/RandExps-p-" << machine.exp << "-P" << P.my_num();
            if (thread_num)
                ss << "-" << thread_num;
            outf.open(ss.str().c_str(), ios::out | ios::binary);
        }

        left_todo = machine.nRandExpsPerThread;

        timers["Generating RandExp(l)"].start();
        while (left_todo > 0) {
        #ifdef BENCHMARK
            rand_pool.clear();
            rand_pool.seekg(ios::beg);
            triple_pool.clear();
            triple_pool.seekg(ios::beg);
        #endif

            int this_loop = amortize;
            if (this_loop > left_todo) {
                this_loop = left_todo;
                Sh_merge.resize(2*this_loop);
                merge.resize(2*this_loop);
                Sh_batch_a.resize(this_loop);
                Sh_batch_b.resize(this_loop);
                Sh_batch_c.resize(this_loop);
                Sh_temp.resize(this_loop);
                temp.resize(this_loop);
                Sh_res.resize(this_loop*(machine.exp + 1));
            }

            for (int i = 0; i < this_loop; i++) {
                Sh_merge[i].input(rand_pool, READABLE);
                Sh_res[this_loop + i] = Sh_merge[i];    // [R]
                Sh_merge[i + this_loop].input(rand_pool, READABLE);
                Sh_res[i] = Sh_merge[i + this_loop];

                Sh_batch_a[i].input(triple_pool, READABLE);
                Sh_batch_b[i].input(triple_pool, READABLE);
                Sh_batch_c[i].input(triple_pool, READABLE);

                Sh_merge[i].sub(Sh_merge[i], Sh_batch_a[i]);
                Sh_merge[i + this_loop].sub(Sh_merge[i + this_loop], Sh_batch_b[i]);
            }

            MC.POpen_Begin(merge, Sh_merge, P);
            MC.POpen_End(merge, Sh_merge, P);

            for (int i = 0; i < this_loop; i++) {
                Sh_temp[i] = Sh_batch_c[i] + merge[i]*Sh_batch_b[i] + merge[i + this_loop]*Sh_batch_a[i];
                Sh_temp[i].add(Sh_temp[i], merge[i]*merge[i + this_loop], playerone, alphai);
            }

            MC.POpen_Begin(temp, Sh_temp, P);
            MC.POpen_End(temp, Sh_temp, P);

            for (int i = 0; i < this_loop; i++) {
                if (temp[i].is_zero()) {
                    throw Offline_Check_Error("Meet a bad zero a");
                }
                temp[i].invert();
                Sh_res[i] *= temp[i];   // [R^-1]
            }

            for (int i = 2; i <= machine.exp; i++) {
                for (int j = 0; j < this_loop; j++) {
                    Sh_batch_a[j].input(triple_pool, READABLE);
                    Sh_batch_b[j].input(triple_pool, READABLE);
                    Sh_batch_c[j].input(triple_pool, READABLE);

                    Sh_merge[j] = Sh_res[this_loop + j];
                    Sh_merge[j].sub(Sh_merge[j], Sh_batch_a[j]);
                    Sh_merge[j + this_loop] = Sh_res[this_loop*(i - 1) + j];
                    Sh_merge[j + this_loop].sub(Sh_merge[j + this_loop], Sh_batch_b[j]);
                }
                MC.POpen_Begin(merge, Sh_merge, P);
                MC.POpen_End(merge, Sh_merge, P);
                for (int j = 0; j < this_loop; j++) {
                    Sh_res[this_loop*i + j] = Sh_batch_c[j] + merge[j]*Sh_batch_b[j] + merge[j + this_loop]*Sh_batch_a[j];
                    Sh_res[this_loop*i + j].add(Sh_res[this_loop*i + j], merge[j]*merge[j + this_loop], playerone, alphai);
                }
            }

            if (machine.output) {
                for (int i = 0; i < this_loop; i++) {
                    for (int j = 0; j <= machine.exp; j++) {
                        Sh_res[i + this_loop*j].output(outf, READABLE);
                    }
                }
            }
            left_todo -= this_loop;
        }
        timers["Generating RandExp(l)"].stop();
        if (machine.output)
            outf.close();

        rand_pool.close();
        triple_pool.close();

        total += machine.nRandExpsPerThread;
        break;
    default:
        throw not_implemented();
    }

    cout << "[*] Could save " << 1e-9 * a.report_size(CAPACITY) << " GB" << endl;
    timers.insert(EC.timers.begin(), EC.timers.end());
    timers["Networking"] = P.timer;
}

template <class FD>
size_t PairwiseGenerator<FD>::report_size(ReportType type)
{
    size_t res = a.report_size(type) + b.report_size(type) + c.report_size(type);
    for (auto m : multipliers)
        res += m->report_size(type);
    res += producer.report_size(type);
    res += randProducer.report_size(type);

    res += multipliers[0]->report_volatile();
    res += volatile_memory + C.report_size(type);
    res += ciphertexts.get_max_length() + cleartexts.get_max_length();
    res += EC.report_size(type) + EC.volatile_memory;
    res += b_mod_q.report_size(type);
    return res;
}

template <class FD>
size_t PairwiseGenerator<FD>::report_sent()
{
    return P.sent + global_player.sent;
}

template <class FD>
void PairwiseGenerator<FD>::report_size(ReportType type, MemoryUsage& res)
{
    multipliers[0]->report_size(type, res);
    res.add("shares",
            a.report_size(type) + b.report_size(type) + c.report_size(type));
    res.add("producer", producer.report_size(type));

    res.add("Rand() producer", randProducer.report_size(type));
    res.add("rand", r.report_size(type));

    res.add("my ciphertexts", C.report_size(CAPACITY));
    res.add("serialized ciphertexts", ciphertexts.get_max_length());
    res.add("serialized cleartexts", cleartexts.get_max_length());
    res.add("generator volatile", volatile_memory);
    res.add("b mod p", b_mod_q.report_size(type));
    res += EC.memory_usage;
}

template class PairwiseGenerator<FFT_Data>;
template class PairwiseGenerator<P2Data>;
