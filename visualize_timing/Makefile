CXX = g++
CXXFLAGS = -std=c++11

main: visualize.cpp matplotlibcpp.h
	$(CXX) $(CXXFLAGS) visualize.cpp -I/usr/include/python2.7 -lpython2.7 -o visualize

clean:
	rm visualize
