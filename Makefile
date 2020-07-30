CXX = g++
CXXFLAGS = -std=c++11

main: main.o gradient_planner_lib.o
	$(CXX) $(CXXFLAGS) -o main main.o gradient_planner_lib.o
main.o: main.cpp gradient_planner_lib.h
	$(CXX) $(CXXFLAGS) -c main.cpp
gradient_planner_lib.o: gradient_planner_lib.h

