# PDCE

## Preface

This repository provides the implementation of [our PDCE protocol](https://). This implementation is based on [SPDZ-2](https://github.com/bristolcrypto/SPDZ-2). To fulfill our additional functionalities, we also realized several types of offline tuples (**Rand**, **Rand2** and **RandExp**). 

We also implements command line interfaces related to our differential privacy expriments in this repository. See [README.md](https://github.com/saftoes/pdce/tree/master/fms/README.md) for more details.

## Requirements
- OS: Ubuntu 16.04 LTS
- GCC, version 7.4.0
- MPIR library, version 2.7.2
- libsodium library, version 1.0.11
- NTL library, **version 11.0.0 (exact)**
- valgrind, version 3.14.0
- [farmhash](https://github.com/google/farmhash) (which are used as the implementation of no-crypto hash)
- OPENSSL, version 1.1.1b


## To compile PDCE
Download this repository and run the following commands in shell.
```bash
cd pdce/
make clean && make all
```

For further compilation details about SPDZ, see [README.md](https://github.com/bristolcrypto/SPDZ-2) in [SPDZ-2](https://github.com/bristolcrypto/SPDZ-2).

## To benchmark PDCE

### Offline phase of PDCE
Besides **Triple** generation realized by Low Gear protocol in [SPDZ-2](https://github.com/bristolcrypto/SPDZ-2), additional three types of tuples (**Rand**, **Rand2** and **RandExp**) are also generated under the framework of Low Gear. To benchmark their generation, one can run `pairwise-offline.x` with the following parameters:

| Parameters                    | Description
| ---------                     | -------- 
| -N [N_PARTIES]                | Number of parties (default: 2)
| -p [ID]                       | ID of this party (from 0 to N_PARTIES - 1)
| -x [N_THREADS]                | Number of threads used (default: 1)
| -o                            | Write generated tuples to files
| -h [HOST]                     | IPv4 host of party 0 (default: localhost)
| -pn [PORT_NUM]                | Base port number (default: 5000)
| -f [FIELD_SIZE]               | Field size of finite field (default: 55)
| -s [STATISTIC_SEC]            | Statistical security parameter for underlying Low Gear protocol (default: 40)
| -n [N_TRIPLES]                | Number of **Triple**s to be generated (default: 131072)
| -r                            | To generate **Rand** tuples instead of default **Triple**
| -nr [N_RANDS]                 | Number of **Rand**s to be generated (default: 131072)
| -r2                           | To generate **Rand2** tuples instead of default **Triple**
| -nr2 [N_RAND2S]               | Number of **Rand2**s to be generated (default: 131072)
| -rx                           | To generate **RandExp** tuples instead of default **Triple**
| -nrx [N_RANDEXPS]             | Number of **RandExp**s to be generated (default: 131072)
| -E [TAU]                      | Value of tau, the parameter for **RandExp** generation

Shares of output tuples, FHE parameters and MAC key share can be found in folders named `[N_PARTIES]-[FIELD_SIZE]-[STATISTIC_SEC]` in `pdce/Player-Data/`.

#### Some examples
To generate 1024 **Triple**s in 5-parties setting (LAN/WAN), in which each party works with 16 threads and the field size of 70, run
```bash
# On party 0
./pairwise-offline.x -o -h [Party_0_HOST] -x 16 -f 70 -N 5 -p 0 -n 1024
# On party 1
./pairwise-offline.x -o -h [Party_0_HOST] -x 16 -f 70 -N 5 -p 1 -n 1024
# On party 2
./pairwise-offline.x -o -h [Party_0_HOST] -x 16 -f 70 -N 5 -p 2 -n 1024
# On party 3
./pairwise-offline.x -o -h [Party_0_HOST] -x 16 -f 70 -N 5 -p 3 -n 1024
# On party 4
./pairwise-offline.x -o -h [Party_0_HOST] -x 16 -f 70 -N 5 -p 4 -n 1024
```

To generate 1024 **Rand**s in 5-parties setting (LAN/WAN), in which each party works with 16 threads and the field size of 70, run
```bash
# On party 0
./pairwise-offline.x -o -h [Party_0_HOST] -x 16 -f 70 -N 5 -p 0 -nr 1024 -r
# On party 1
./pairwise-offline.x -o -h [Party_0_HOST] -x 16 -f 70 -N 5 -p 1 -nr 1024 -r
# On party 2
./pairwise-offline.x -o -h [Party_0_HOST] -x 16 -f 70 -N 5 -p 2 -nr 1024 -r
# On party 3
./pairwise-offline.x -o -h [Party_0_HOST] -x 16 -f 70 -N 5 -p 3 -nr 1024 -r
# On party 4
./pairwise-offline.x -o -h [Party_0_HOST] -x 16 -f 70 -N 5 -p 4 -nr 1024 -r
```

To generate 1024 **RandExp**s (with tau = 30) in 5-parties setting (LAN/WAN), in which each party works with 16 threads and the field size of 70, run
```bash
# On party 0
./pairwise-offline.x -o -h [Party_0_HOST] -x 16 -f 70 -N 5 -p 0 -nrx 1024 -rx -E 30
# On party 1
./pairwise-offline.x -o -h [Party_0_HOST] -x 16 -f 70 -N 5 -p 1 -nrx 1024 -rx -E 30
# On party 2
./pairwise-offline.x -o -h [Party_0_HOST] -x 16 -f 70 -N 5 -p 2 -nrx 1024 -rx -E 30
# On party 3
./pairwise-offline.x -o -h [Party_0_HOST] -x 16 -f 70 -N 5 -p 3 -nrx 1024 -rx -E 30
# On party 4
./pairwise-offline.x -o -h [Party_0_HOST] -x 16 -f 70 -N 5 -p 4 -nrx 1024 -rx -E 30
```
Note that **Rand2** and **RandExp** rely on pre-computed **Triple** and **Rand**. So before generating **Rand2** and **RandExp**, you should always generate sufficient **Triple**s and **Rand**s previously.


### Online phase of PDCE
#### Simple local example
To test the functionality of each protocol of online phase, say, one data party with around $10^9$ unique items and two computation parties securely approximating the unique items with $2000$ trails in one local machine.

Run the offline phase to generate the necessary offline data first. Then open a new cmd and start the Names server which establishs the communication network between Computation parties.
```bash
./Server.x 2 5000 1
```

Open two new cmd and run the computation parties
```bash
./CP.x  -lgp 70 -np 2 -p 0 -sec 40 -ndp 1 -oM 2000 -uN 1000000000 -lan 1 -nt 1 
./CP.x  -lgp 70 -np 2 -p 1 -sec 40 -ndp 1 -oM 2000 -uN 1000000000 -lan 1 -nt 1
```

In the cmd run Names server, run the data party
```bash
./DP.x -ndp 1 -ncp 2 -p 1 -lgp 70 -oM 2000 -uN 1000000000 -nt 1
```

#### parameters for Names server
```bash
./Server.x {number of computation parties} 5000 {number of threads}
```

#### parameters for Computation party
```bash
./CP.x  -lgp {length of the prime} -np {number of computation parties} -p {party number} -sec {security level} -ndp {number of data parties} -oM {number of flajolet-martin trails} -uN {around number of unique items} -lan {whether in lan environment} -nt {number of threads} -h {ip of Data party}
```

#### parameters for Data party
```bash
./DP.x -ndp {number of data parties} -ncp {number of computation parties} -p {party number} -lgp {length of the prime} -oM {number of flajolet-martin trails} -uN {around number of unique items} -ip {ip files}
```

#### Simple remote 5 computation party, 1 data party example, running 16 threads
Note, you should generate enough offline data with "-x 16" parameters first.

Create the file ,e.g., called ip_addr.data. Write the computation parties' ip addresses in it from party 0 to party 4, seperating them pressing **Enter** key. Note, you need to press **Enter** after the party 4's ip address also.

In machine 0, start the Names server which establishs the communication network between Computation parties.
```bash
./Server.x 5 5000 16
```

In machine 1-6 run the computation party. 
```bash
./CP.x  -lgp 70 -np 5 -p {party number: from 0 to 4} -sec 40 -ndp 1 -oM 2000 -uN 1000000000 -lan 0 -nt 16 -h {ip where the data party resides}
```
Note if the machines are in local network, set the -lan to 1, otherwise to 0. This is because in WAN mode, we parallerize the sending and receiving of the sockets to reduce the communication delay. For details, see the comments in source file Protocol/ComputationParty.h .

In machine 0, run the data party.
```bash
./DP.x -ndp 1 -ncp 5 -p 1 -lgp 70 -oM 2000 -uN 1000000000 -ip ip_addr.data
```

#### network parameter tuning
To full use of the bandwith which can make better your online benchmark performance in WAN mode. You can larger the tcp buffer with the following commands
```bash
sudo su 
echo 'net.core.wmem_max=12582912' >> /etc/sysctl.conf &&
echo 'net.core.rmem_max=12582912' >> /etc/sysctl.conf &&
echo 'net.ipv4.tcp_rmem= 40960 807380 12582912' >> /etc/sysctl.conf &&
echo 'net.ipv4.tcp_wmem= 40960 807380 12582912' >> /etc/sysctl.conf &&
echo 'net.ipv4.tcp_window_scaling = 1' >> /etc/sysctl.conf &&
echo 'net.ipv4.tcp_timestamps = 1' >> /etc/sysctl.conf &&
echo 'net.ipv4.tcp_sack = 1' >> /etc/sysctl.conf &&
echo 'net.core.netdev_max_backlog = 5000' >> /etc/sysctl.conf &&
sysctl -p
su - {the username of the machine}
```
