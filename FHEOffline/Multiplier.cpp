// (C) 2018 University of Bristol. See License.txt

/*
 * Multiplier.cpp
 *
 */

#include <FHEOffline/Multiplier.h>
#include "FHEOffline/PairwiseGenerator.h"
#include "FHEOffline/PairwiseMachine.h"

#include <assert.h>

template <class FD>
Multiplier<FD>::Multiplier(int offset, PairwiseGenerator<FD>& generator) :
    generator(generator), machine(generator.machine),
    P(generator.global_player, offset),
    num_players(generator.global_player.num_players()),
    my_num(generator.global_player.my_num()),
    other_pk(machine.other_pks[(my_num + num_players - offset) % num_players]),
    other_enc_alpha(machine.enc_alphas[(my_num + num_players - offset) % num_players]),
    timers(generator.timers),
    C(machine.pk), mask(machine.pk),
    product_share(machine.setup<FD>().FieldD), rc(machine.pk),
    volatile_capacity(0)
{
    product_share.allocate_slots(machine.setup<FD>().params.p0() << 64);
}

template <class FD>
void Multiplier<FD>::multiply_and_add(Plaintext_<FD>& res,
        const Ciphertext& enc_a, const Plaintext_<FD>& b)
{
    Rq_Element bb(enc_a.get_params(), evaluation, evaluation);
    bb.from_vec(b.get_poly());
    multiply_and_add(res, enc_a, bb);
}

template <class FD>
void Multiplier<FD>::multiply_and_add(Plaintext_<FD>& res,
        const Ciphertext& enc_a, const Rq_Element& b)
{
    // Working as P_j
    PRNG G;
    G.ReSeed();
    timers["Ciphertext multiplication"].start();
    C.mul(enc_a, b);
    timers["Ciphertext multiplication"].stop();
    timers["Mask randomization"].start();
    product_share.randomize(G);
    bigint B = 6 * machine.setup<FD>().params.get_R();
    B *= machine.setup<FD>().FieldD.get_prime();
    B <<= machine.sec;
    // slack
    B *= NonInteractiveProof::slack(machine.sec,
            machine.setup<FD>().params.phi_m());
    B <<= machine.extra_slack;
    rc.generateUniform(G, 0, B, B);
    timers["Mask randomization"].stop();
    timers["Encryption"].start();
    other_pk.encrypt(mask, product_share, rc);
    timers["Encryption"].stop();
    timers["Multiplied ciphertext sending"].start();
    octetStream o;
    mask += C;
    mask.pack(o);
    P.reverse_exchange(o);
    C.unpack(o);

    // Working as P_i
    timers["Multiplied ciphertext sending"].stop();
    timers["Decryption"].start();
    res -= product_share;
    machine.sk.decrypt_any(product_share, C);
    res += product_share;
    timers["Decryption"].stop();
    memory_usage.update("multiplied ciphertext", C.report_size(CAPACITY));
    memory_usage.update("mask ciphertext", mask.report_size(CAPACITY));
    memory_usage.update("product shares", product_share.report_size(CAPACITY));
    memory_usage.update("masking random coins", rc.report_size(CAPACITY));
}

template <class FD>
void Multiplier<FD>::multiply_alpha_and_add(Plaintext_<FD>& res,
        const Rq_Element& b)
{
    multiply_and_add(res, other_enc_alpha, b);
}

template <class FD>
size_t Multiplier<FD>::report_size(ReportType type)
{
    return C.report_size(type) + mask.report_size(type)
            + product_share.report_size(type) + rc.report_size(type);
}

template <class FD>
void Multiplier<FD>::report_size(ReportType type, MemoryUsage& res)
{
    (void)type;
    res += memory_usage;
}


template class Multiplier<FFT_Data>;
template class Multiplier<P2Data>;



template <class FD>
void Multiplier<FD>::multiply_and_add(Plaintext_<FD>& res, const Ciphertext& enc_a,
        const Plaintext_<FD>& b,
        const vector<typename FD::T>& rand_coeffs,
        const typename FD::T& alphai)
{
    // Working as P_j
    PRNG G;
    G.ReSeed();
    timers["Ciphertext multiplication"].start();
    C.mul(enc_a, b);
    timers["Ciphertext multiplication"].stop();
    timers["Mask randomization"].start();
    product_share.randomize(G);
    bigint B = 6 * machine.setup<FD>().params.get_R();
    B *= machine.setup<FD>().FieldD.get_prime();
    B <<= machine.sec;
    // slack
    B *= NonInteractiveProof::slack(machine.sec,
            machine.setup<FD>().params.phi_m());
    B <<= machine.extra_slack;
    rc.generateUniform(G, 0, B, B);
    timers["Mask randomization"].stop();
    timers["Encryption"].start();
    other_pk.encrypt(mask, product_share, rc);
    timers["Encryption"].stop();
    timers["Multiplied ciphertext sending"].start();
    octetStream o;
    mask += C;

    // Step 5/6 in the Input protocol of Overdrive
    typename FD::T rho;
    rho.assign_zero();
    typename FD::T sigma_i;
    sigma_i.assign_zero();
    typename FD::T temp;

    for (unsigned int i = 0; i < product_share.num_slots(); ++i)
    {
        temp.mul(rand_coeffs[i], b.element(i));
        rho += temp;
        temp.mul(rand_coeffs[i], product_share.element(i));
        sigma_i += temp;
    }

    mask.pack(o);
    rho.pack(o);
    sigma_i.pack(o);
    P.reverse_exchange(o);
    C.unpack(o);
    rho.unpack(o);
    sigma_i.unpack(o);

    // Working as P_i
    timers["Multiplied ciphertext sending"].stop();
    timers["Decryption"].start();
    res -= product_share;
    machine.sk.decrypt_any(product_share, C);
    res += product_share;
    timers["Decryption"].stop();

    timers["Rand checking"].start();
    typename FD::T check;
    check.assign_zero();
    for (unsigned int i = 0; i < product_share.num_slots(); ++i)
    {
        temp.mul(rand_coeffs[i], product_share.element(i));
        check += temp;
    }
    check -= sigma_i;
    temp.mul(alphai, rho);
    check.sub(temp, check);

    assert( check.is_zero() == true );

    timers["Rand checking"].stop();


    memory_usage.update("multiplied ciphertext", C.report_size(CAPACITY));
    memory_usage.update("mask ciphertext", mask.report_size(CAPACITY));
    memory_usage.update("product shares", product_share.report_size(CAPACITY));
    memory_usage.update("masking random coins", rc.report_size(CAPACITY));
}

template <class FD>
void Multiplier<FD>::multiply_alpha_and_add(Plaintext_<FD>& res,
        const Plaintext_<FD>& b,
        const vector<typename FD::T>& rand_coeffs,
        const typename FD::T& alphai)
{
    multiply_and_add(res, other_enc_alpha, b, rand_coeffs, alphai);
}
