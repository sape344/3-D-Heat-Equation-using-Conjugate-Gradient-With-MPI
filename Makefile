


CC = g++
CFLAGS  = -O3 -time -std=c++2a

all: serial_heat_map paralel_heat_map

serial_heat_map:  serial_heat_map.cpp 
	$(CC) serial_heat_map.cpp $(CFLAGS) -o serial_heat_map.o 

paralel_heat_map: paralel_heat_map.cpp 
	mpic++ paralel_heat_map.cpp $(CFLAGS) -o paralel_heat_map.o

clean: 
	$(RM) count *.o *~
