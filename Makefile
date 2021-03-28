CXX = g++
CXXFLAGS = -std=c++11

main: main.o gradient_planner_lib3.o
	$(CXX) $(CXXFLAGS) main.o gradient_planner_lib3.o -I/usr/include/python2.7 -lpython2.7 -o main
main.o: main.cpp gradient_planner_lib3.h matplotlibcpp.h
	$(CXX) $(CXXFLAGS) -I/usr/include/python2.7 -lpython2.7 -c main.cpp
	# $(CXX) $(CXXFLAGS) -c main.cpp
gradient_planner_lib3.o: gradient_planner_lib3.h

clean:
	rm *.o
	rm main
