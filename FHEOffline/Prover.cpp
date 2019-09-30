// (C) 2018 University of Bristol. See License.txt


#include "Prover.h"

#include "Tools/random.h"


template <class FD, class U>
Prover<FD,U>::Prover(Proof& proof, const FD& FieldD) :
  volatile_memory(0)
{
  s.resize(proof.V, proof.pk->get_params());
  y.resize(proof.V, FieldD);
#ifdef LESS_ALLOC_MORE_MEM
  s.allocate_slots(bigint(1) << proof.B_rand_length);
  y.allocate_slots(bigint(1) << proof.B_plain_length);
  t = s[0];
  z = y[0];
  // extra limb to prevent reallocation
  t.allocate_slots(bigint(1) << (proof.B_rand_length + 64));
  z.allocate_slots(bigint(1) << (proof.B_plain_length + 64));
#endif
}

template <class FD, class U>
void Prover<FD,U>::Stage_1(const Proof& P, octetStream& ciphertexts,
    const AddableVector<Ciphertext>& c,
    const FHE_PK& pk, bool Diag, bool binary)
{
    size_t allocate = 3 * c.size() * c[0].report_size(USED);
    ciphertexts.resize_precise(allocate);
    ciphertexts.reset_write_head();
    c.pack(ciphertexts);

    int V=P.V;

//  AElement<T> AE;
//  ZZX rd;
//  ZZ pr=(*AE.A).prime();
//  ZZ bd=B_plain/(pr+1);
    PRNG G;
    G.ReSeed();
    Random_Coins rc(pk.get_params());
    Ciphertext ciphertext(pk.get_params());
    ciphertexts.store(V);
    for (int i=0; i<V; i++)
    {
//      AE.randomize(Diag,binary);
//      rd=RandPoly(phim,bd<<1);
//      y[i]=AE.plaintext()+pr*rd;


        y[i].randomize(G, P.B_plain_length, Diag, binary);
        s[i].resize(3, P.phim);
        s[i].generateUniform(G, P.B_rand_length);
        rc.assign(s[i][0], s[i][1], s[i][2]);
        pk.encrypt(ciphertext,y[i],rc);   // a_i = Enc_pk(y_i, s_i)
        ciphertext.pack(ciphertexts);
    }
}


template <class FD, class U>
bool Prover<FD,U>::Stage_2(Proof& P, octetStream& cleartexts,
                        const vector<U>& x,
                        const Proof::Randomness& r,
                        const vector<int>& e)
{
    size_t allocate = P.V * P.phim
        * (5 + numBytes(P.plain_check) + 3 * (5 + numBytes(P.rand_check)));
    cleartexts.resize_precise(allocate);
    cleartexts.reset_write_head();

    unsigned int i,k;
    int j,ee;
#ifndef LESS_ALLOC_MORE_MEM
    AddableVector<bigint> z;
    AddableMatrix<bigint> t;
#endif
    cleartexts.reset_write_head();
    cleartexts.store(P.V);

    // NOTE: indices of matrix M_e are controlled by i, k
    for (i=0; i<P.V; i++)
    {   z=y[i];
        t=s[i];
        for (k=0; k<P.sec; k++)
        {   j=(i+1)-(k+1);
            // NOTE: compressed form of the matrix M_e
            //       elements in M_e is either 0 or 1
            if (j<0 || j>=(int) P.sec) { ee=0; }
            else                       { ee=e[j]; }

            if (ee!=0)
	          {
              z += x[j];
	            t += r[j];
            }
        }
        if (not P.check_bounds(z, t, i))
            return false;
        z.pack(cleartexts);
        t.pack(cleartexts);
    }
#ifndef LESS_ALLOC_MORE_MEM
    volatile_memory = t.report_size(CAPACITY) + z.report_size(CAPACITY);
#endif
#ifdef PRINT_MIN_DIST
    cout << "Minimal distance (log) " << log2(P.dist) << ", compare to " <<
        log2(P.plain_check.get_d() / pow(2, P.B_plain_length))  << endl;
#endif
    return true;
}



/* This is the non-interactive version using the ROM
*/
template <class FD, class U>
size_t Prover<FD,U>::NIZKPoK(Proof& P, octetStream& ciphertexts, octetStream& cleartexts,
                        const FHE_PK& pk,
                        const AddableVector<Ciphertext>& c,
                        const vector<U>& x,
                        const Proof::Randomness& r,
                        bool Diag,bool binary)
{
  vector<int> e(P.sec);

//  AElement<T> AE;
//  for (i=0; i<P.sec; i++)
//    { AE.assign(x.at(i));
//      if (!AE.to_type(0))
//         { cout << "Error in making x[i]" << endl;
//           cout << i << endl;
//         }
//    }

  bool ok=false;
  int cnt=0;
  while (!ok)
    { cnt++;
      Stage_1(P,ciphertexts,c,pk,Diag,binary);
      P.get_challenge(e, ciphertexts);
      // Check check whether we are OK, or whether we should abort
      ok = Stage_2(P,cleartexts,x,r,e);
    }
  if (cnt > 1)
      cout << "\t\tNumber iterations of prover = " << cnt << endl;
  return report_size(CAPACITY) + volatile_memory;
}


template<class FD, class U>
size_t Prover<FD,U>::report_size(ReportType type)
{
  size_t res = 0;
  for (unsigned int i = 0; i < s.size(); i++)
    res += s[i].report_size(type);
  for (unsigned int i = 0; i < y.size(); i++)
    res += y[i].report_size(type);
#ifdef LESS_ALLOC_MORE_MEM
  res += z.report_size(type) + t.report_size(type);
#endif
  return res;
}


template<class FD, class U>
void Prover<FD, U>::report_size(ReportType type, MemoryUsage& res)
{
  res.update("prover s", s.report_size(type));
  res.update("prover y", y.report_size(type));
#ifdef LESS_ALLOC_MORE_MEM
  res.update("prover z", z.report_size(type));
  res.update("prover t", t.report_size(type));
#endif
  res.update("prover volatile", volatile_memory);
}


template class Prover<FFT_Data, Plaintext_<FFT_Data> >;
template class Prover<FFT_Data, AddableVector<bigint> >;

template class Prover<P2Data, Plaintext_<P2Data> >;
template class Prover<P2Data, AddableVector<bigint> >;
