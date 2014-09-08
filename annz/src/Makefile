CPP = g++

all: annz_net annz_train annz_test

annz_net: annz_net.cpp
	$(CPP) -g -o annz_net annz_net.cpp

network.o: network.h network.cpp util.o
	$(CPP) -g -c network.cpp

train_network.o: train_network.h train_network.cpp
	$(CPP) -g -c train_network.cpp

test_network.o: test_network.h test_network.cpp
	$(CPP) -g -c test_network.cpp

util.o: util.cpp util.h
	$(CPP) -g -c util.cpp

train.o: train.cpp train.h
	$(CPP) -g -c train.cpp

test.o: test.cpp test.h
	$(CPP) -g -c test.cpp

annz_test: annz_test.cpp network.o util.o test.o test_network.o
	$(CPP) -g -o annz_test annz_test.cpp test_network.o network.o util.o test.o

annz_train: annz_train.cpp network.o util.o train.o train_network.o
	$(CPP) -g -o annz_train annz_train.cpp train_network.o network.o util.o train.o

clean:
	rm *.o *~
