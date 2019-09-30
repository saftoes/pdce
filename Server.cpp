// (C) 2018 University of Bristol. See License.txt

/*
 * Server.cpp
 *
 */

#include "Networking/Server.h"
#include <iostream>
using namespace std;

int main(int argc, char** argv)
{
	if (argc != 4) {
    	cout << "Usage: ./Server.x n_computation_parties portbase n_threads\n";
    	return 0;
  	}

  	constexpr int port_increase = 75;
    constexpr int communication_multiplier = 2;
  	int nmachines = atoi(argv[1]);
  	int PortnumBase = atoi(argv[2]);
  	int n_threads = atoi(argv[3]);

  	for(int i=0; i<n_threads*communication_multiplier; i++){
  		Server(nmachines, PortnumBase+i*port_increase).start();
  	}
}
