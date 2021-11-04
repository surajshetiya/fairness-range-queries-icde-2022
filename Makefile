CXX=clang++
CFLAGS=-I.

$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c

fair_1d: fairness_1d.o
	$(CXX) -ggdb -o fairness_1d.o fairness_1d.cpp