HEADERS =-I./ -I/usr/local/include/ -I/${boost_1_62_0_path}/boost_1_62_0
LIBS =-L./ -L/usr/local/lib -L/${boost_1_62_0_path}/boost_1_62_0/stage/lib

Trajquery:  trajquery.cpp
	mpicxx $(HEADERS) $(LIBS) -O3 -std=c++11 -o Trajquery trajquery.cpp 
clean:
	rm  Trajquery
