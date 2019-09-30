// (C) 2018 University of Bristol. See License.txt

/*
 * PairwiseSetup.cpp
 *
 */

#include <FHEOffline/PairwiseSetup.h>
#include "FHE/NoiseBounds.h"
#include "FHE/NTL-Subs.h"
#include "Math/Setup.h"
#include "FHEOffline/Proof.h"

template <class FD>
void PairwiseSetup<FD>::init(const Player& P, int sec, int plaintext_length,
    int& extra_slack)
{
    sec = max(sec, 40);
    cout << "Finding parameters for security " << sec << " and field size 2^"
            << plaintext_length << endl;
    PRNG G;
    G.ReSeed();
    //dirname = PREP_DIR;

    stringstream file;
    file << PREP_DIR << "/" << P.num_players() << "-";
    file << plaintext_length << "-" << sec << "/";

    file >> dirname;

    cout << "[*] Current dir: " << endl;
    cout << "\t" << dirname << endl;

    octetStream o;
    if (P.my_num() == 0)
    {
        extra_slack =
            generate_semi_setup(plaintext_length, sec, params, FieldD, true);
        params.pack(o);
        FieldD.pack(o);
        o.store(extra_slack);
        P.send_all(o);
    }
    else
    {
        P.receive_player(0, o);
        params.unpack(o);
        FieldD.unpack(o);
        FieldD.init_field();
        o.get(extra_slack);
    }

    stringstream ss;
    ss << dirname << "/Player-MAC-Keys-P" << P.my_num();

    ifstream test_exist;
    test_exist.open(ss.str().c_str());

    alpha = FieldD;
    if (!test_exist.good()) {
        alpha.randomize(G, Diagonal);
        alphai = alpha.element(0);
    }
    else {
        cout << "[*] MAC key exists, generation ignored ..." << endl;
        int dummy;
        test_exist >> dummy;
        alphai.input(test_exist, true);
        alpha.randomize(G, Diagonal);
        for (unsigned int i = 0; i < alpha.num_slots(); ++i) {
            alpha.set_element(i, alphai);
        }
        test_exist.close();
    }
}

template class PairwiseSetup<FFT_Data>;
template class PairwiseSetup<P2Data>;
