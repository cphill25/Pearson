all: pearson3

pearson3 : pearson.3.3.cpp
	g++ -O3 -o pearson3 pearson.3.3.cpp

clean:
	/bin/rm pearson3 *.temp
