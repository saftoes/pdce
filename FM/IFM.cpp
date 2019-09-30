// (C) 2018 University of NKU. Free for used
// Author: stoneboat@mail.nankai.edu.cn

/*
 * IFM.cpp
 *
 */

#include "IFM.h"
#include <iostream>

void IntegerFlajoletMartin::init_fake_IFM(){
	std::vector<uint64_t> fm_sketch(M);
	generate_fake_FM(N, M, fm_sketch);

	FlajoletMartin::randomSketchGenerator rsg;
	unsigned char* rnd_num = new unsigned char[M*w];
	rsg.random_bytes(M*w,rnd_num);
	  
	int k = 0;
	uint64_t tmp = 0;
	bool is_one;
	for(int m=0; m<M; m++){
		tmp = fm_sketch[m];
		k = w;
		while(k--){
			is_one = (bool)(tmp & 1);
			tmp = tmp >> 1;
			if(is_one){
				IFS_group[m][w-k-1] = 1 + (int)rnd_num[m*w+w-k-1];
			}else{
				IFS_group[m][w-k-1] = 0;
			}
		}

	}
	fill_IFS = true;
}

void IntegerFlajoletMartin::init_IFS(){
	IFS_group.resize(M);

	for(int i=0; i<M; i++){
		IFS_group[i].resize(w);
	}
}

void IntegerFlajoletMartin::check_fake_IFM(){
	init_fake_IFM();
	if(!fill_IFS){
		std::cout<<" have not fill the IFS structure, exit with no action\n";
	}

	double R = 0;
	double ret = 0;

	uint64_t tmp = 0;
	bool bit;
	for(int m=0; m<M; m++){
		tmp = 0;
		for(int i=w-1; i>=0; i--){
			bit = (bool)IFS_group[m][i] ; 
			tmp = tmp | ((int)bit);
			tmp = tmp << 1;
		}
		tmp = tmp >> 1;
		R += ls0(tmp);
	}

	R /= M;
	ret = pow(2.0, R)/PHI;
	std::cout<<" number of unique items is "<<N<<"\t estimate size is "<< ret<<" with "<< M <<"\t trails\n";

}




