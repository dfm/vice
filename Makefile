demo: src/demo.cc include/vice/*.h
	g++ -O3 -march=native -std=c++14 -march=native -Iinclude -Ivendor/starry -Ivendor/starry/lib/eigen_3.3.3 -Ivendor/starry/lib/boost_1_66_0 -Ivendor/starry/lib/LBFGSpp/include src/demo.cc -o bin/demo
