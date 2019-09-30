include CONFIG

USE_NTL = 1  # spdz-2 offline
MOD = -DMAX_MOD_SZ=6 # spdz-2 offline
LDLIBS := -lntl $(LDLIBS)

MATH = $(patsubst %.cpp,%.o,$(wildcard Math/*.cpp))

TOOLS = $(patsubst %.cpp,%.o,$(wildcard Tools/*.cpp))

NETWORK = $(patsubst %.cpp,%.o,$(wildcard Networking/*.cpp))

AUTH = $(patsubst %.cpp,%.o,$(wildcard Auth/*.cpp))

PROCESSOR = $(patsubst %.cpp,%.o,$(wildcard Processor/*.cpp))

FM = $(patsubst %.cpp,%.o,$(wildcard FM/*.cpp))

PROTOCOL = $(patsubst %.cpp,%.o,$(wildcard Protocol/*.cpp)) $(FM)

#ifeq ($(USE_NTL),1)
FHEOFFLINE = $(patsubst %.cpp,%.o,$(wildcard FHEOffline/*.cpp FHE/*.cpp))
#endif

# OT stuff needs GF2N_LONG, so only compile if this is enabled
ifeq ($(USE_GF2N_LONG),1)
OT = $(patsubst %.cpp,%.o,$(filter-out OT/OText_main.cpp,$(wildcard OT/*.cpp)))
OT_EXE = ot.x ot-offline.x
endif

COMMON = $(MATH) $(TOOLS) $(NETWORK) $(AUTH)
COMPLETE = $(COMMON) $(PROCESSOR) $(FHEOFFLINE) $(TINYOTOFFLINE) $(OT)

LIB = libSPDZ.a
LIBSIMPLEOT = SimpleOT/libsimpleot.a

# used for dependency generation
OBJS = $(COMPLETE)
DEPS := $(OBJS:.o=.d)

##################################################################### PDCE-Protocol Start###############################################
FM_LDLIBS := -lfarmhash -lcrypto $(LDLIBS)

## our compilation
all: PDCE

PDCE: PDCE-offline PDCE-online

PDCE-offline: pairwise-offline.x
PDCE-online: CP.x DP.x Server.x

CP.x: CP.cpp $(COMMON) $(PROCESSOR) $(PROTOCOL)
	$(CXX) $(CFLAGS) CP.cpp -o CP.x $(COMMON) $(PROCESSOR) $(PROTOCOL) $(FM_LDLIBS)

DP.x: DP.cpp $(COMMON) $(PROCESSOR) $(PROTOCOL)
	$(CXX) $(CFLAGS) DP.cpp -o DP.x $(COMMON) $(PROCESSOR) $(PROTOCOL) $(FM_LDLIBS)

FM-check.x: FM-check.cpp $(COMMON) $(FM) 
	$(CXX) $(CFLAGS) FM-check.cpp -o FM-check.x $(COMMON) $(FM) $(FM_LDLIBS) 

##################################################################### PDCE-Protocol End###############################################



spdz-all: gen_input online offline externalIO

ifeq ($(USE_NTL),1)
spdz-all: overdrive she-offline
endif

-include $(DEPS)

%.o: %.cpp
	$(CXX) $(CFLAGS) -MMD -c -o $@ $<

online: Fake-Offline.x Server.x Player-Online.x Check-Offline.x

offline: $(OT_EXE) Check-Offline.x

gen_input: gen_input_f2n.x gen_input_fp.x 

externalIO: client-setup.x bankers-bonus-client.x bankers-bonus-commsec-client.x

she-offline: Check-Offline.x spdz2-offline.x

overdrive: simple-offline.x pairwise-offline.x cnc-offline.x

Fake-Offline.x: Fake-Offline.cpp $(COMMON) $(PROCESSOR)
	$(CXX) $(CFLAGS) -o $@ $^ $(LDLIBS)

Check-Offline.x: Check-Offline.cpp $(COMMON) $(PROCESSOR)
	$(CXX) $(CFLAGS) Check-Offline.cpp -o Check-Offline.x $(COMMON) $(PROCESSOR) $(LDLIBS)

Server.x: Server.cpp $(COMMON)
	$(CXX) $(CFLAGS) Server.cpp -o Server.x $(COMMON) $(LDLIBS)

Player-Online.x: Player-Online.cpp $(COMMON) $(PROCESSOR)
	$(CXX) $(CFLAGS) Player-Online.cpp -o Player-Online.x $(COMMON) $(PROCESSOR) $(LDLIBS)

ifeq ($(USE_GF2N_LONG),1)
ot.x: $(OT) $(COMMON) OT/OText_main.cpp
	$(CXX) $(CFLAGS) -o $@ $^ $(LDLIBS) $(LIBSIMPLEOT)

ot-check.x: $(OT) $(COMMON)
	$(CXX) $(CFLAGS) -o ot-check.x OT/BitVector.o OT/OutputCheck.cpp $(COMMON) $(LDLIBS)

ot-bitmatrix.x: $(OT) $(COMMON) OT/BitMatrixTest.cpp
	$(CXX) $(CFLAGS) -o ot-bitmatrix.x OT/BitMatrixTest.cpp OT/BitMatrix.o OT/BitVector.o $(COMMON) $(LDLIBS)

ot-offline.x: $(OT) $(COMMON) ot-offline.cpp
	$(CXX) $(CFLAGS) -o $@ $^ $(LDLIBS) $(LIBSIMPLEOT)
endif

check-passive.x: $(COMMON) check-passive.cpp
	$(CXX) $(CFLAGS) -o $@ $^ $(LDLIBS)

gen_input_f2n.x: Scripts/gen_input_f2n.cpp $(COMMON)
	$(CXX) $(CFLAGS) Scripts/gen_input_f2n.cpp	-o gen_input_f2n.x $(COMMON) $(LDLIBS)

gen_input_fp.x: Scripts/gen_input_fp.cpp $(COMMON)
	$(CXX) $(CFLAGS) Scripts/gen_input_fp.cpp	-o gen_input_fp.x $(COMMON) $(LDLIBS)

read_fp.x: Scripts/gen_input_fp.cpp $(COMMON)
	$(CXX) $(CFLAGS) Scripts/read_fp.cpp	-o read_fp.x $(COMMON) $(LDLIBS)

client-setup.x: client-setup.cpp $(COMMON) $(PROCESSOR)
	$(CXX) $(CFLAGS) -o $@ $^ $(LDLIBS)

bankers-bonus-client.x: ExternalIO/bankers-bonus-client.cpp $(COMMON) $(PROCESSOR)
	$(CXX) $(CFLAGS) -o $@ $^ $(LDLIBS)

bankers-bonus-commsec-client.x: ExternalIO/bankers-bonus-commsec-client.cpp $(COMMON) $(PROCESSOR)
	$(CXX) $(CFLAGS) -o $@ $^ $(LDLIBS)

	
#ifeq ($(USE_NTL),1)
simple-offline.x: $(COMMON) $(FHEOFFLINE) simple-offline.cpp
	$(CXX) $(CFLAGS) -o $@ $^ $(LDLIBS)

pairwise-offline.x: $(COMMON) $(FHEOFFLINE) pairwise-offline.cpp
	$(CXX) $(CFLAGS) -o $@ $^ $(LDLIBS)

cnc-offline.x: $(COMMON) $(FHEOFFLINE) cnc-offline.cpp
	$(CXX) $(CFLAGS) -o $@ $^ $(LDLIBS)

spdz2-offline.x: $(COMMON) $(FHEOFFLINE) spdz2-offline.cpp
	$(CXX) $(CFLAGS) -o $@ $^ $(LDLIBS)
#endif


clean:
	-rm */*.o */*.key */*.x *.o */*.d *.d *.x core.* *.a gmon.out

clean-log:
	-rm */*.log *.log



