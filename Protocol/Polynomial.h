// (C) 2018 University of NKU. Free for used

/*
 * Polynomial.h
 *
 */

#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_


/*
 * @Params of the Class DP
 *
 */
#include "Math/bigint.h"
#include "Math/gfp.h"

#include <NTL/ZZ_p.h>
#include <mpirxx.h>

#include <utility>
#include <vector>



namespace PSUCA{
  namespace Polynomial{
  	enum Poly_type
  	{
  		poly_zero = 1,
  		poly_cdb = 2,
  		poly_lookup = 3

  	};

    class Interpolatation_Polynomial{
    private:
    	bool set;
      	int npoints;

      	bigint p;
      	std::vector<std::pair<NTL::ZZ,NTL::ZZ>> points; // first for x, second for y
      	std::vector<NTL::ZZ_p> ZZ_coeff;
          
    public:
    	Interpolatation_Polynomial();
    	// load coeff from the file
    	Interpolatation_Polynomial(const bigint& _p, const std::string& folder_name, int tau, Poly_type type);
	   	Interpolatation_Polynomial(const bigint& _p, int _npoints, const std::vector<std::pair<int,int>>& _points);

	    void interpolate();
	    NTL::ZZ_p evaluate(const int& x);
	    void output(std::vector<gfp>& polyf);
 
    	/*
    	*	Conversion Utility function
    	*/ 
    	void mpz_to_ZZ(const mpz_t& mpz_var, NTL::ZZ& zz_var);
    	void bigint_to_ZZ(const bigint& bigint_var, NTL::ZZ& zz_var);
    	void ZZ_to_bigint(const NTL::ZZ& zz_var, bigint& bigint_var);
    	void ZZ_to_mpz(const NTL::ZZ& zz_var, mpz_t& mpz_var);
    	void string_to_ZZ(const char* str, int len, NTL::ZZ& number);

    	/*
    	*	File Utility function
    	*/ 
    	static bool is_directory(const std::string& path);
    	static bool create_directory(const std::string& path, mode_t mode);

    	static void load_points(std::string folder_name, int _npoints, Poly_type type, std::vector<std::pair<int,int>>& _points);
    	static void write_points(std::string folder_name, int tau, Poly_type type);

    	void write_coeff(std::string folder_name, Poly_type type);

    };
  }
}






#endif /* POLYNOMIAL_H_ */