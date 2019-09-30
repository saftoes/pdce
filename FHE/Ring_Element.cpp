// (C) 2018 University of Bristol. See License.txt


#include "FHE/Ring_Element.h"
#include "Exceptions/Exceptions.h"
#include "FHE/FFT.h"

void reduce_step(vector<modp>& aa,int i,const FFT_Data& FFTD)
{ modp temp=aa[i];
  for (int j=0; j<FFTD.phi_m(); j++)
    { switch (FFTD.Phi()[j])
        { case 0:
             break;
          case 1:
             Sub(aa[i-FFTD.phi_m()+j],aa[i-FFTD.phi_m()+j],temp,FFTD.get_prD());
             break;
          case -1:
             Add(aa[i-FFTD.phi_m()+j],aa[i-FFTD.phi_m()+j],temp,FFTD.get_prD());
             break;
          case 2:
             Sub(aa[i-FFTD.phi_m()+j],aa[i-FFTD.phi_m()+j],temp,FFTD.get_prD());
             Sub(aa[i-FFTD.phi_m()+j],aa[i-FFTD.phi_m()+j],temp,FFTD.get_prD());
             break;
          case -2:
             Add(aa[i-FFTD.phi_m()+j],aa[i-FFTD.phi_m()+j],temp,FFTD.get_prD());
             Add(aa[i-FFTD.phi_m()+j],aa[i-FFTD.phi_m()+j],temp,FFTD.get_prD());
             break;
          case 3:
             Sub(aa[i-FFTD.phi_m()+j],aa[i-FFTD.phi_m()+j],temp,FFTD.get_prD());
             Sub(aa[i-FFTD.phi_m()+j],aa[i-FFTD.phi_m()+j],temp,FFTD.get_prD());
             Sub(aa[i-FFTD.phi_m()+j],aa[i-FFTD.phi_m()+j],temp,FFTD.get_prD());
             break;
          case -3:
             Add(aa[i-FFTD.phi_m()+j],aa[i-FFTD.phi_m()+j],temp,FFTD.get_prD());
             Add(aa[i-FFTD.phi_m()+j],aa[i-FFTD.phi_m()+j],temp,FFTD.get_prD());
             Add(aa[i-FFTD.phi_m()+j],aa[i-FFTD.phi_m()+j],temp,FFTD.get_prD());
             break;
          default:
             throw not_implemented();
        }
    }
}


Ring_Element::Ring_Element(const FFT_Data& fftd,RepType r)
{ 
  FFTD=&fftd;
  rep=r;          
  assign_zero();
}


void Ring_Element::assign_zero()
{
  element.resize((*FFTD).phi_m());
  for (int i=0; i<(*FFTD).phi_m(); i++)
    { assignZero(element[i],(*FFTD).get_prD()); }
}


void Ring_Element::assign_one()
{
  element.resize((*FFTD).phi_m());
  modp fill;
  if (rep==polynomial) { assignZero(fill,(*FFTD).get_prD()); }
  else                 { assignOne(fill,(*FFTD).get_prD()); }
  for (int i=1; i<(*FFTD).phi_m(); i++)
    { element[i]=fill; }
  assignOne(element[0],(*FFTD).get_prD());
}



void Ring_Element::negate()
{
  for (int i=0; i<(*FFTD).phi_m(); i++)
    { Negate(element[i],element[i],(*FFTD).get_prD()); }
}



void add(Ring_Element& ans,const Ring_Element& a,const Ring_Element& b)
{
  if (a.rep!=b.rep)   { throw rep_mismatch(); }
  if (a.FFTD!=b.FFTD) { throw pr_mismatch();  }  
  ans.partial_assign(a);
  for (int i=0; i<(*ans.FFTD).phi_m(); i++)
    { Add(ans.element[i],a.element[i],b.element[i],(*a.FFTD).get_prD()); }
}



void sub(Ring_Element& ans,const Ring_Element& a,const Ring_Element& b)
{
  if (a.rep!=b.rep)   { throw rep_mismatch(); }
  if (a.FFTD!=b.FFTD) { throw pr_mismatch();  }
  ans.partial_assign(a);
  for (int i=0; i<(*ans.FFTD).phi_m(); i++)
    { Sub(ans.element[i],a.element[i],b.element[i],(*a.FFTD).get_prD()); }
}



void mul(Ring_Element& ans,const Ring_Element& a,const Ring_Element& b)
{
  if (a.rep!=b.rep)   { throw rep_mismatch(); }
  if (a.FFTD!=b.FFTD) { throw pr_mismatch();  }
  ans.partial_assign(a);
  if (ans.rep==evaluation)
    { // In evaluation representation, so we can just multiply componentwise
      for (int i=0; i<(*ans.FFTD).phi_m(); i++)
        { Mul(ans.element[i],a.element[i],b.element[i],(*a.FFTD).get_prD()); }
    }
  else if ((*ans.FFTD).get_twop()!=0)
    { // This is the case where m is not a power of two

      // Here we have to do a poly mult followed by a reduction
      // We could be clever (e.g. use Karatsuba etc), but instead
      // we do the school book method followed by term re-writing

      // School book mult
      vector<modp> aa(2*(*a.FFTD).phi_m());
      for (int i=0; i<2*(*a.FFTD).phi_m(); i++)
        { assignZero(aa[i],(*a.FFTD).get_prD()); }
      modp temp;
      for (int i=0; i<(*a.FFTD).phi_m(); i++)
        { for (int j=0; j<(*a.FFTD).phi_m(); j++)
	    { Mul(temp,a.element[i],b.element[j],(*a.FFTD).get_prD()); 
              int k=i+j;
              Add(aa[k],aa[k],temp,(*a.FFTD).get_prD());
            }
        }
      // Now apply reduction, assumes Ring.poly is monic
      for (int i=2*(*a.FFTD).phi_m()-1; i>=(*a.FFTD).phi_m(); i--)
        { reduce_step(aa,i,*a.FFTD);
          assignZero(aa[i],(*a.FFTD).get_prD()); 
        }
     // Now stick into answer
     for (int i=0; i<(*ans.FFTD).phi_m(); i++)
       { ans.element[i]=aa[i]; }
    }
  else if ((*ans.FFTD).get_twop()==0)
    { // m a power of two case
      Ring_Element aa(*ans.FFTD,ans.rep);
      modp temp;
      for (int i=0; i<(*ans.FFTD).phi_m(); i++)
        { for (int j=0; j<(*ans.FFTD).phi_m(); j++)
            { Mul(temp,a.element[i],b.element[j],(*a.FFTD).get_prD());
              int k=i+j;
              if (k>=(*ans.FFTD).phi_m())
                 { k-=(*ans.FFTD).phi_m();
                   Negate(temp,temp,(*a.FFTD).get_prD());
                 }
              Add(aa.element[k],aa.element[k],temp,(*a.FFTD).get_prD());
            }
        }
      ans=aa;
    }
  else
    { throw not_implemented(); }
}


void mul(Ring_Element& ans,const Ring_Element& a,const modp& b)
{
  ans.partial_assign(a);
  for (int i=0; i<(*ans.FFTD).phi_m(); i++)
    { Mul(ans.element[i],a.element[i],b,(*a.FFTD).get_prD()); }
}



void Ring_Element::randomize(PRNG& G,bool Diag)
{
  if (Diag==false)
    { for (int i=0; i<(*FFTD).phi_m(); i++) 
       { element[i].randomize(G,(*FFTD).get_prD()); }
    }
  else 
    { element[0].randomize(G,(*FFTD).get_prD());
      if (rep==polynomial)
        { for (int i=1; i<(*FFTD).phi_m(); i++) 
            { assignZero(element[i],(*FFTD).get_prD()); }
        }
      else
        { for (int i=1; i<(*FFTD).phi_m(); i++) 
            { element[i]=element[0]; }  
        }
    }
}


ostream& operator<<(ostream& s,const Ring_Element& e)
{
  s << "[";
  if (e.rep==polynomial) { s << "P "; }
  else                   { s << "E "; }
  bigint te;
  for (int i=0; i<(*e.FFTD).phi_m(); i++) 
    { to_bigint(te,e.element[i],(*e.FFTD).get_prD());
      s << te << " ";
    }
  s << "]";
  return s;
}

/*
 * Important: e must already have been initialised with the correct
 * Zp_data, or this may fail!
 */
istream& operator>>(istream& s, Ring_Element& e)
{
  vector<modp> elem;
  bigint cur;
  RepType rep;
  int ch = s.get();
  while (isspace(ch))
  {
      ch = s.get();
  }
  
  if (ch != '[')
     { throw IO_Error("Bad Ring_Element input: no '['"); }
  ch = s.get();
  if (ch == 'P')
      rep = polynomial;
  else if (ch == 'E')
      rep = evaluation;
  else
    { throw IO_Error("Error reading Ring_Element : ch="+to_string((char)(ch))); }
  e.rep = rep;
  ch = s.peek();
  while (isspace(ch))
  {
      s.get();
      ch = s.peek();
  }
  modp te;
  while (ch != ']')
  {
      if (!(s >> cur))
          { throw IO_Error("Bad Ring_Element input"); }
      
      to_modp(te,cur,(*e.FFTD).get_prD());
      elem.push_back(te);
      
      ch = s.peek();
      while (isspace(ch))
      {
          s.get();
          ch = s.peek();
      }
  }
  s.get();
  
  e.element = elem;
  return s;
}


void Ring_Element::change_rep(RepType r)
{ 
  if (rep==r) { return; }
  if (r==evaluation)
    { rep=evaluation;
      if ((*FFTD).get_twop()==0)
        { // m a power of two variant
          FFT_Iter2(element,(*FFTD).phi_m(),(*FFTD).get_root(0),(*FFTD).get_prD());
	}
      else
        { // Non m power of two variant and FFT enabled
          vector<modp> fft((*FFTD).m());
          BFFT(fft,element,*FFTD);
          for (int i=0; i<(*FFTD).phi_m(); i++)
	    { element[i]=fft[(*FFTD).p(i)]; }
	}
    }
  else
    { rep=polynomial;
      if ((*FFTD).get_twop()==0)
	{ // m a power of two variant
          modp root2;
          Sqr(root2,(*FFTD).get_root(1),(*FFTD).get_prD());
          FFT_Iter(element, (*FFTD).phi_m(),root2,(*FFTD).get_prD());
          modp w;
          w = (*FFTD).get_iphi();
          for (int i=0; i<(*FFTD).phi_m(); i++)
            { Mul(element[i], element[i], w, (*FFTD).get_prD());
              Mul(w, w, (*FFTD).get_root(1),(*FFTD).get_prD());
            }
        }
      else
        { // Non power of 2 m variant and FFT enabled
          vector<modp> fft((*FFTD).m());
          for (int i=0; i<(*FFTD).m(); i++) 
            { assignZero(fft[i],(*FFTD).get_prD()); }
          for (int i=0; i<(*FFTD).phi_m(); i++) 
	    { fft[(*FFTD).p(i)]=element[i]; }
          BFFT(fft,fft,*FFTD,false);
          // Need to reduce fft mod Phi_m
          for (int i=(*FFTD).m()-1; i>=(*FFTD).phi_m(); i--)
            { reduce_step(fft,i,*FFTD); }
          for (int i=0; i<(*FFTD).phi_m(); i++) 
	    { element[i]=fft[i]; }
	}
    }
}


bool Ring_Element::equals(const Ring_Element& a) const
{
  if (rep!=a.rep)   { throw rep_mismatch(); }
  if (*FFTD!=*a.FFTD) { throw pr_mismatch();  }
  for (int i=0; i<(*FFTD).phi_m(); i++)
    { if (!areEqual(element[i],a.element[i],(*FFTD).get_prD())) { return false; } }
  return true;
}


void Ring_Element::from_vec(const vector<bigint>& v)
{
  RepType t=rep;
  rep=polynomial;
  bigint tmp;
  for (int i=0; i<(*FFTD).phi_m(); i++)
    {
      tmp = v[i];
      element[i].convert_destroy(tmp, FFTD->get_prD());
    }
  change_rep(t);
//  cout << "RE:from_vec<bigint>:: " << *this << endl;
}

void Ring_Element::from_vec(const vector<int>& v)
{
  RepType t=rep;
  rep=polynomial;
  for (int i=0; i<(*FFTD).phi_m(); i++)
    { to_modp(element[i],v[i],(*FFTD).get_prD()); }
  change_rep(t);
//  cout << "RE:from_vec<int>:: " << *this << endl;
}

template <class T>
void Ring_Element::from(const Generator<T>& generator)
{
  RepType t=rep;
  rep=polynomial;
  T tmp;
  for (int i=0; i<(*FFTD).phi_m(); i++)
    {
      generator.get(tmp);
      element[i].convert_destroy(tmp, (*FFTD).get_prD());
    }
  change_rep(t);
}

ConversionIterator Ring_Element::get_iterator() const
{
  if (rep != polynomial)
    throw runtime_error("simple iterator only available in polynomial represention");
  return {element, (*FFTD).get_prD()};
}

RingReadIterator Ring_Element::get_copy_iterator() const
{
  return *this;
}

RingWriteIterator Ring_Element::get_write_iterator()
{
  return *this;
}

vector<bigint>  Ring_Element::to_vec_bigint() const
{
  vector<bigint> v;
  to_vec_bigint(v);
  return v;
}


void Ring_Element::to_vec_bigint(vector<bigint>& v) const
{
  v.resize(FFTD->phi_m());
  if (rep==polynomial)
     { for (int i=0; i<(*FFTD).phi_m(); i++)
         { to_bigint(v[i],element[i],(*FFTD).get_prD()); }
     }
  else
     { Ring_Element a=*this;
       a.change_rep(polynomial);
       for (int i=0; i<(*FFTD).phi_m(); i++)
         { to_bigint(v[i],a.element[i],(*FFTD).get_prD()); }
     }
}





modp Ring_Element::get_constant() const
{
  if (rep==polynomial)
     { return element[0]; }
  Ring_Element a=*this;
  a.change_rep(polynomial);
  return a.element[0]; 
}



void store(octetStream& o,const vector<modp>& v,const Zp_Data& ZpD)
{
  o.store((int)v.size());
  for (unsigned int i=0; i<v.size(); i++)
     { v[i].pack(o,ZpD); }
}


void get(octetStream& o,vector<modp>& v,const Zp_Data& ZpD)
{
  unsigned int length;
  o.get(length);
  v.resize(length);
  for (unsigned int i=0; i<length; i++)
    { v[i].unpack(o,ZpD); }
}


void Ring_Element::pack(octetStream& o) const
{
  check_size();
  o.store(rep);
  store(o,element,(*FFTD).get_prD());
}


void Ring_Element::unpack(octetStream& o)
{
  unsigned int a;
  o.get(a);
  rep=(RepType) a;
  check_rep();
  get(o,element,(*FFTD).get_prD());
  check_size();
}


void Ring_Element::check_rep()
{
  if (rep != evaluation and rep != polynomial)
    throw runtime_error("invalid representation");
}


void Ring_Element::check_size() const
{
  if ((int)element.size() != FFTD->phi_m())
    throw runtime_error("invalid element size");
}

void Ring_Element::output(ostream& s) const
{
  s.write((char*)&rep, sizeof(rep));
  auto size = element.size();
  s.write((char*)&size, sizeof(size));
  for (auto& x : element)
    x.output(s, FFTD->get_prD(), false);
}


void Ring_Element::input(istream& s)
{
  s.read((char*)&rep, sizeof(rep));
  check_rep();
  auto size = element.size();
  s.read((char*)&size, sizeof(size));
  element.resize(size);
  for (auto& x : element)
    x.input(s, FFTD->get_prD(), false);
}


void Ring_Element::check(const FFT_Data& FFTD) const
{
  if (&FFTD != this->FFTD)
    throw params_mismatch();
}


size_t Ring_Element::report_size(ReportType type) const
{
  if (type == CAPACITY)
    return sizeof(modp) * element.capacity();
  else
    return sizeof(mp_limb_t) * (*FFTD).get_prD().get_t() * element.size();
}

template void Ring_Element::from(const Generator<bigint>& generator);
template void Ring_Element::from(const Generator<int>& generator);
