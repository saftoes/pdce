// (C) 2018 University of NKU. Free for used

/*
 * Polynomial.cpp
 *
 */
#include "Polynomial.h"
#include "Exceptions/Exceptions.h"

#include <sstream> 
#include <string> 
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <math.h> 

#include <NTL/mat_ZZ_p.h>
#include <NTL/vec_ZZ_p.h>

using namespace std;

namespace PSUCA{
  namespace Polynomial{

  	using namespace NTL;

  	Interpolatation_Polynomial::Interpolatation_Polynomial(){
  		set = false;
  	}

  	Interpolatation_Polynomial::Interpolatation_Polynomial(const bigint& _p, int _npoints, const std::vector<std::pair<int,int>>& _points){
  		p = _p;
  		npoints = _npoints;
  		points.resize(npoints);
  		for(int i=0; i<npoints; i++){
  			points[i].first = _points[i].first;
  			points[i].second = _points[i].second;
  		}
      set = true;
      
  		interpolate();
  	}

  	Interpolatation_Polynomial::Interpolatation_Polynomial(const bigint& _p, const std::string& folder_name, int tau, Poly_type type){
  		p = _p;
  		ZZ zz_tmp;
  		bigint_to_ZZ(p, zz_tmp);
  		ZZ_p::init(zz_tmp); 

  		int n;
  		std::string file_path;


  		switch(type){
  			case poly_zero:
  				file_path = folder_name + "coeff_polyzero_" + to_string(tau) + ".txt";
  				n = tau + 1;
  				break;
  			case poly_cdb:
  				file_path = folder_name + "coeff_polycdb_" + to_string(tau) + ".txt";
  				n = tau + 1;
  				break;
  			case poly_lookup:
  				file_path = folder_name + "coeff_polylookup_" + to_string(tau) + ".txt";
  				n = pow(2,tau);
  				break;
  			default:
  				throw bad_value();
  		}

  		ZZ_coeff.resize(n);
  		ifstream file(file_path.c_str());
  		if(!file.good()){
  			throw file_missing(file_path);
  		}

  		npoints = n;
  		string tmp;
  		ZZ tmp_z;
  		ZZ_p tmp_zp;
  		for(int i=0; i<n; i++){
  			tmp.clear();
  			std::getline(file,tmp);

  			string_to_ZZ(tmp.c_str(), tmp.length(), tmp_z);
  			tmp_zp = to_ZZ_p(tmp_z);
  			ZZ_coeff[i] = tmp_zp;
  		}    

  		set = true;
  	}


  	void Interpolatation_Polynomial::interpolate(){
  		if(!set){
  			throw runtime_error("the polynomial parameter has not been initialized\n");
  		}
  		// get the prime in the form of ZZ_p
  		ZZ zz_tmp;
  		bigint_to_ZZ(p, zz_tmp);
  		ZZ_p::init(zz_tmp); 

  		// Generate a vandermonde Matrix
  		int n = npoints;
  		mat_ZZ_p Vander;
  		ZZ_p dimension;
  		Vander.SetDims(n, n);
  		for(int i=0; i<n; i++){
  			ZZ_p tmp(1);
  			ZZ_p multiplier;
  			multiplier = to_ZZ_p(points[i].first);
  			for(int j=0; j<n; j++){
  				Vander[i][j] = tmp;
  				tmp = tmp*multiplier;
  			}
  		}

  		// Generate the Y vector
  		vec_ZZ_p Y;
  		Y.SetLength(n);
  		for(int i=0; i<n; i++){
  			Y[i] = to_ZZ_p(points[i].second);
  		}

  		mat_ZZ_p inverse_Vander = inv(Vander);

  		vec_ZZ_p coeff = inverse_Vander * Y;
  		ZZ_coeff.resize(n);
  		for(int i=0; i<n; i++){
  			ZZ_coeff[i] = coeff[i];
  		}
  	}

  	ZZ_p Interpolatation_Polynomial::evaluate(const int& x){
  		if(!set){
  			throw runtime_error("the polynomial parameter has not been initialized\n");
  		}
  		int n = npoints;
		ZZ_p tmp(1);
		ZZ_p result(0);
		ZZ_p multiplier(x);

		for(int i=0; i<n; i++){
			result += tmp*ZZ_coeff[i];
			tmp = tmp*multiplier;
		}

		return result;
  	}

  	void Interpolatation_Polynomial::output(std::vector<gfp>& polyf){
  		if(!set){
  			throw runtime_error("the polynomial parameter has not been initialized\n");
  		}
  		polyf.clear();
  		polyf.resize(npoints);

  		ZZ tmp_z;
  		bigint tmp_b;
  		for(int i=0; i<npoints; i++){
  			tmp_z = rep(ZZ_coeff[i]);
  			ZZ_to_bigint(tmp_z,tmp_b);
  			to_gfp(polyf[i],tmp_b);
  		}
  	}


  	void Interpolatation_Polynomial::mpz_to_ZZ(const mpz_t& mpz_var, NTL::ZZ& zz_var){
  		char* mpz_num = mpz_get_str(NULL,10,mpz_var);
  		std::stringstream ss;
  		ss << mpz_num;
  		ss >> zz_var;
  	}

  	void Interpolatation_Polynomial::ZZ_to_mpz(const NTL::ZZ& zz_var, mpz_t& mpz_var){
  		std::stringstream ss;
  		ss << zz_var;
  		mpz_set_str(mpz_var, ss.str().c_str(),10);
  	}

  	void Interpolatation_Polynomial::bigint_to_ZZ(const bigint& bigint_var, NTL::ZZ& zz_var){
  		std::string bigint_num = bigint_var.get_str(10);
  		std::stringstream ss;
  		ss << bigint_num.c_str();
  		ss >> zz_var;
  	}

  	void Interpolatation_Polynomial::ZZ_to_bigint(const NTL::ZZ& zz_var, bigint& bigint_var){
  		std::stringstream ss;
  		ss << zz_var;
  		bigint_var.set_str(ss.str().c_str(),10);
  	}

  	void Interpolatation_Polynomial::string_to_ZZ(const char* str, int len, NTL::ZZ& number){
  		number = conv<ZZ>(str[0]-48); //convert from ascii
	    for(int i = 1; i < len; i++)
	    {
	        number *= 10; // base
	        number += conv<ZZ>(str[i]-48);
	    }
  	}

  	void Interpolatation_Polynomial::load_points(std::string folder_name, int _npoints, Poly_type type, std::vector<std::pair<int,int>>& _points){
  		std::string file_path;
  		int n;
  		switch(type){
  			case poly_zero:
  				file_path = folder_name + "points_polyzero_" + to_string(_npoints) + ".txt";
  				n = _npoints + 1;
  				break;
  			case poly_cdb:
  				file_path = folder_name + "points_polycdb_" + to_string(_npoints) + ".txt";
  				n = _npoints + 1;
  				break;
  			case poly_lookup:
  				file_path = folder_name + "points_polylookup_" + to_string(_npoints) + ".txt";
  				n = pow(2,_npoints);
  				break;
  			default:
  				throw bad_value();

  		}
  		_points.clear();
  		_points.resize(n);
  		ifstream file(file_path.c_str());
  		if(!file.good()){
  			throw file_missing(file_path);
  		}
  		for(int i=0; i<n; i++){
  			file >> _points[i].first;
  			file >> _points[i].second; 
  		}
  	}

  	bool Interpolatation_Polynomial::is_directory(const std::string& path)
	{
	    struct stat sb;
	    
	    if (stat(path.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))
	    {
	        return true;
	    }
	    return false;
	}

	// named pricinple points_{functionality}_{tau}
  	void Interpolatation_Polynomial::write_points(std::string folder_name, int tau, Poly_type type){
  		if(!is_directory(folder_name)){
  			if(!create_directory(folder_name,S_IRWXU)){
    			throw file_error("could not mkdir " + folder_name);
    		}
	    }
	    int n;
	    if(type == poly_zero){
	    	std::string file_path = folder_name + "points_polyzero_" + to_string(tau) + ".txt";
	    	n = tau+1;
	  		ofstream file(file_path.c_str());
	  		file << "0\t0\n";
	  		for(int i=1; i<n; i++){
	  			file<<i<<"\t1\n";
	  		}
	  		file.close();
	    }

	   	if(type == poly_cdb){
	    	std::string file_path = folder_name + "points_polycdb_" + to_string(tau) + ".txt";
	    	n = tau+1;
	  		ofstream file(file_path.c_str());
	  		file << "0\t1\n";
	  		for(int i=1; i<n; i++){
	  			file<<i<<"\t0\n";
	  		}
	  		file.close();
	    }

	    if(type == poly_lookup){
	    	std::string file_path = folder_name + "points_polylookup_" + to_string(tau) + ".txt";
	    	n = pow(2,tau);
	  		ofstream file(file_path.c_str());
	  		file << "0\t"<<tau<<"\n";
	  		int ptr = 1;
	  		for(int i=0; i<tau; i++){
	  			int size = pow(2,i);
	  			for(int j=0; j<size; j++){
	  				file<<ptr<<"\t"<<(tau-i-1)<<"\n";
	  				ptr++;
	  			}
	  		}
	  		file.close();
	    }
  		
  	}


  	void Interpolatation_Polynomial::write_coeff(std::string folder_name, Poly_type type){
  		if(!is_directory(folder_name)){
  			if(!create_directory(folder_name,S_IRWXU)){
    			throw file_error("could not mkdir " + folder_name);
    		}
	    }
	    if(!set){
  			throw runtime_error("the polynomial parameter has not been initialized\n");
  		}


  		std::string file_path;
  		int n = ZZ_coeff.size(); 
  		switch(type){
  			case poly_zero:
  				n = n -1;
  				file_path = folder_name + "coeff_polyzero_" + to_string(n) + ".txt";
  				break;
  			case poly_cdb:
  				n = n -1;
  				file_path = folder_name + "coeff_polycdb_" + to_string(n) + ".txt";
  				break;
  			case poly_lookup:
  				n = log2(n);
  				file_path = folder_name + "coeff_polylookup_" + to_string(n) + ".txt";
  				break;
  			default:
  				throw bad_value();
  		}

  		ofstream file(file_path.c_str());
  		for(size_t i=0; i<ZZ_coeff.size(); i++){
  			file<<ZZ_coeff[i]<<"\n";
  		}
  		file.close();
  	}


  	bool Interpolatation_Polynomial::create_directory(const std::string& path, mode_t mode){
  		if (mkdir(path.data(),mode) != 0) {
	        return false;
	    }
	    return true;
  	}
  }
}