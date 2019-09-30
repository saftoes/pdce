// (C) 2018 University of Bristol. See License.txt

/*
 * OfflineMachineBase.cpp
 *
 */

#include <Tools/OfflineMachineBase.h>
#include "Tools/int.h"

#include <string>
using namespace std;


OfflineMachineBase::OfflineMachineBase() :
        server(0), my_num(0), nplayers(0), nthreads(0),
        ntriples(0), nTriplesPerThread(0),
        nrands(0), nRandsPerThread(0),
        nrand2s(0), nRand2sPerThread(0),
        output(0)
{
}

OfflineMachineBase::~OfflineMachineBase()
{
    if (server)
        delete server;
}

void OfflineMachineBase::parse_options(int argc, const char** argv)
{
    opt.add(
        "2", // Default.
        0, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "Number of parties (default: 2).", // Help description.
        "-N", // Flag token.
        "--nparties" // Flag token.
    );
    opt.add(
        "", // Default.
        1, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "This player's number, 0/1 (required).", // Help description.
        "-p", // Flag token.
        "--player" // Flag token.
    );
    opt.add(
        "1",
        0,
        1,
        0,
        "Number of threads (default: 1).",
        "-x",
        "--nthreads"
    );
    opt.add(
        "", // Default.
        0, // Required?
        0, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "Write results to files.", // Help description.
        "-o", // Flag token.
        "--output" // Flag token.
    );

    opt.add(
        "131072", // Default.
        0, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "Number of triples (default: 131072).", // Help description.
        "-n", // Flag token.
        "--ntriples" // Flag token.
    );
    opt.add(
        "131072", // Default.
        0, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "Number of Rand() values (default: 131072).", // Help description.
        "-nr", // Flag token.
        "--nrands" // Flag token.
    );
    opt.add(
        "131072", // Default.
        0, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "Number of Rand2() values (default: 131072).", // Help description.
        "-nr2", // Flag token.
        "--nrand2s" // Flag token.
    );
    opt.add(
        "131072", // Default.
        0, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "Number of RandInRange(l) (default: 131072).", // Help description.
        "-nrl", // Flag token.
        "--nranges" // Flag token.
    );
    opt.add(
        "39", // Default.
        0, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "Parameter for RandInRange(l) (default: 39).", // Help description.
        "-L", // Flag token.
        "--rangeLen" // Flag token.
    );
    opt.add(
        "131072", // Default.
        0, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "Number of RandExp(l) (default: 131072).", // Help description.
        "-nrx", // Flag token.
        "--nexps" // Flag token.
    );
    opt.add(
        "15", // Default.
        0, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "Parameter for RandExp(l) (default: 15).", // Help description.
        "-E", // Flag token.
        "--exp" // Flag token.
    );
    opt.add(
        "16384", // Default.
        0, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "Number of RandBits(l) (default: 16384).", // Help description.
        "-nrb", // Flag token.
        "--nbits" // Flag token.
    );

    opt.parse(argc, argv);
    if (!opt.isSet("-p"))
    {
        string usage;
        opt.getUsage(usage);
        cout << usage;
        exit(0);
    }

    opt.get("-p")->getInt(my_num);
    opt.get("-N")->getInt(nplayers);
    opt.get("-x")->getInt(nthreads);

    output = opt.get("-o")->isSet;

    opt.get("-n")->getLongLong(ntriples);
    opt.get("-nr")->getLongLong(nrands);
    opt.get("-nr2")->getLongLong(nrand2s);
    opt.get("-nrl")->getLongLong(nranges);
    opt.get("-L")->getLongLong(l);
    opt.get("-nrx")->getLongLong(nexps);
    opt.get("-E")->getLongLong(exp);
    opt.get("-nrb")->getLongLong(nbits);

#ifdef BENCHMARK
    nTriplesPerThread       = 65536;
    nRandsPerThread         = 65536;
    nRand2sPerThread        = 65536;
    nRandRangesPerThread    = 65536;
    nRandExpsPerThread      = 65536;
    nRandBitsPerThread      = 16384;
    cout << "[*] Warning: working in benchmark mode. See [MY_CFLAGS] in CONFIG ..." << endl;
#else
    nTriplesPerThread       = DIV_CEIL(ntriples, nthreads);
    nRandsPerThread         = DIV_CEIL(nrands, nthreads);
    nRand2sPerThread        = DIV_CEIL(nrand2s, nthreads);
    nRandRangesPerThread    = DIV_CEIL(nranges, nthreads);
    nRandExpsPerThread      = DIV_CEIL(nexps, nthreads);
    nRandBitsPerThread      = DIV_CEIL(nbits, nthreads);
#endif
}

void OfflineMachineBase::start_networking_with_server(string hostname,
        int portnum)
{
    server = Server::start_networking(N, my_num, nplayers, hostname, portnum);
}
