CC=g++  -std=c++14
#CFLAGS+=-Wall -I/usr/local/Cellar/opencv/4.1.0_2/include/opencv4/opencv2 -I/usr/local/Cellar/opencv/4.1.0_2/include/opencv4 -O3 -I/usr/local/include
CFLAGSLINK+= -lgmp -lpthread -L/usr/local/lib/ -lssl -lcrypto
CFLAGSLINK +=-L/usr/local/opt/openssl@1.1/lib
CFLAGSLINK +=-I/usr/local/opt/openssl@1.1/include

sloth : sloth_core.o
	$(CC) $(CFLAGS) $(CFLAGSLINK) sloth_core.cpp -o sloth
sloth_core.o : sloth_core.cpp
	$(CC) sloth_core.cpp $(CFLAGS) -c