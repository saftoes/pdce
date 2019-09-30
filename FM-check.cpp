// (C) 2018 University of NKU. Free for used
// Author: stoneboat@mail.nankai.edu.cn

/*
 * FM-check.cpp
 *
 */

#include "FM/PFM.h"
#include "FM/IFM.h"

#include <iostream>
#include <string>
#include <vector>

#include <sys/time.h>

#include "Tools/ezOptionParser.h"




/******************************************************************************/

int main(int argc, char **argv){
  int check_PFM = 0;
  int check_IFM = 0;


	ez::ezOptionParser opt;

  opt.syntax = "./FM-check.x [OPTIONS]\n";
  opt.example = "./FM-check.x  \n";

  opt.add(
          "0", // Default.
          0, // Required?
          1, // Number of args expected.
          0, // Delimiter if expecting multiple args.
          "check the plaintext FlajoletMartin", // Help description.
          "-p", // Flag token.
          "--plaintext FlajoletMartin" // Flag token. bit-structured plaintext FlajoletMartin
  );
  opt.add(
          "0", // Default.
          0, // Required?
          1, // Number of args expected.
          0, // Delimiter if expecting multiple args.
          "check the integer FlajoletMartin", // Help description.
          "-I", // Flag token.
          "--integer FlajoletMartin" // Flag token. bit-structured plaintext FlajoletMartin
  );
  opt.add(
          "1", // Default.
          0, // Required?
          1, // Number of args expected.
          0, // Delimiter if expecting multiple args.
          "number of threads", // Help description.
          "-x", // Flag token.
          "--nthreads" // Flag token. bit-structured plaintext FlajoletMartin
  );

  opt.parse(argc, (const char **) argv);

  int nthreads;
  opt.get("--nthreads")->getInt(nthreads);
  opt.get("--plaintext FlajoletMartin")->getInt(check_PFM);
  opt.get("--integer FlajoletMartin")->getInt(check_IFM);

  std::string usage;

  if (check_PFM)
  {
      std::cout<< "check the bit-structured plaintext FlajoletMartin\n";
      // read the data
      
      struct timeval t1;
      struct timeval t2;
      double cpu_time_used;

      std::cout<<"test random FM generating\n";
      for (int n = 0; n < 2; n ++){
        for (int m = 0; m < 2; m ++){
          std::cout<<"Experiment parameter\tM:\t"<<FlajoletMartin::M[m]<<"\tN:\t"<<FlajoletMartin::N[n]<<"\tw:\t"<<FlajoletMartin::w[n]<<std::endl;

          gettimeofday(&t1, NULL);
          BitsFlajoletMartin bfm(FlajoletMartin::M[m],FlajoletMartin::N[n],nthreads);
          double napprox = bfm.randomFM_approximateCount();
          gettimeofday(&t2, NULL);
          cpu_time_used = (double)(t2.tv_sec-t1.tv_sec)*1000+(double)(t2.tv_usec-t1.tv_usec)/1000;
          printf("\tcount = %lf \ttime %f (ms) \n",napprox,cpu_time_used);

        }
      }

      std::cout<<"test fake FM generating\n";
      for (int n = 0; n < 3; n ++){
        for (int m = 0; m < 3; m ++){
          gettimeofday(&t1, NULL);
          BitsFlajoletMartin bfm(FlajoletMartin::M[m],FlajoletMartin::N[n],nthreads);
          bfm.check_fake_FM();
          gettimeofday(&t2, NULL);
          cpu_time_used = (double)(t2.tv_sec-t1.tv_sec)*1000+(double)(t2.tv_usec-t1.tv_usec)/1000;
          printf("\tgenerate  %d FlajoletMartin sketch need \ttime %f (ms) \n",FlajoletMartin::M[m],cpu_time_used);

        }
      }     
      std::cout<< "check done\n";
  }


  if (check_IFM)
  {
      std::cout<< "integer-structed oblivious FlajoletMartin sketch\n";
      // read the data
      
      struct timeval t1;
      struct timeval t2;
      double cpu_time_used;

      std::cout<<"test fake OFM generating\n";
      for (int n = 0; n < 3; n ++){
        for (int m = 0; m < 3; m ++){
          gettimeofday(&t1, NULL);
          IntegerFlajoletMartin ifm(FlajoletMartin::M[m],FlajoletMartin::N[n],nthreads);
          ifm.check_fake_IFM();
          gettimeofday(&t2, NULL);
          cpu_time_used = (double)(t2.tv_sec-t1.tv_sec)*1000+(double)(t2.tv_usec-t1.tv_usec)/1000;
          printf("\tgenerate  %d FlajoletMartin sketch need \ttime %f (ms) \n",FlajoletMartin::M[m],cpu_time_used);

        }
      }     
      std::cout<< "check done\n";
  }





}


/******************************************************************************/
