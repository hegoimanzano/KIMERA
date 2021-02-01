# To compile kimera. 
KIMERA: main.cpp
	g++     -Wall -fexceptions -fopenmp -O3 -std=c++11 include/*.h main.cpp src/*.cpp -Iinclude -o Kimera.exe 

# To compile kimera, Other example. COMPLETE PATHS. AT LEAST c++4.9.2 COMPILER IS NECESSARY. SELECT ONE, REMOVE THE REST.
#KIMERA: main.cpp
#	(path-to-compilers)/g++ -L/path-to-compilers/lib64/  -Wl,-rpath=(path-to-compilers)/lib64/   -Wall -fexceptions -fopenmp -std=c++11 -O3 include/*.h  main.cpp src/*.cpp -Iinclude -o Kimera.exe 
	 
# To compile kimera. Other example
#KIMERA: main.cpp 
#	icpc  -Wall -fexceptions -qopenmp -fast -std=c++11   include/*.h main.cpp src/*.cpp -Iinclude -o Kimera.exe

# To compile kimera on Mac
#KIMERA: main.cpp 
#	clang++ -Wall -fexceptions -fopenmp -std=c++11 -O3 main.cpp src/*.cpp -Iinclude -o Kimera.exe -lomp
