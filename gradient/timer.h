//na zacatku mereneho scopu spustit pomoci "Timer timer"
#include <chrono>
#include <iostream>
#include <string>

struct Timer{
    std::chrono::steady_clock::time_point start,end;
    std::chrono::duration<float> duration;
    std::string name;
    
    //konstruktor - zavolano pri vytvoreni
    Timer() : name("Default timer"){
	start = std::chrono::steady_clock::now();
    }
    Timer(const char* str) : name(str){
	start = std::chrono::steady_clock::now();
    }

    //destruktor - zavolano pri zniceni
    ~Timer(){
	end = std::chrono::steady_clock::now();
	duration = end - start;
	float ms = duration.count() * 1000.0f;
	std::cout << name << " took " << ms << "ms\n";
    }
};
