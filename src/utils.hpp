# include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <time.h>
#include <sys/stat.h>
#include <unistd.h>

#ifndef UTILS_HPP
#define UTILS_HPP

using namespace std;

inline vector<string> split(string str, char delimiter){
	vector<string> internal;
	stringstream ss(str);
	string temp;
	while (getline(ss, temp, delimiter)){
		internal.push_back(temp);
	}
	return internal;
}

inline string convertToString(char* a, int size) 
{ 
    int i; 
    string s = ""; 
    for (i = 0; i < size; i++) { 
        s = s + a[i]; 
    } 
    return s; 
} 

inline void currentDateTime(string &date_str, string &time_str) {
    time_t     now = time(0); 
    struct tm  tstruct;
    char       date_buf[11];
    char       time_buf[9];
    tstruct = *localtime(&now);
    strftime(date_buf, sizeof(date_buf), "%Y-%m-%d", &tstruct);
    strftime(time_buf, sizeof(time_buf), "%H-%M-%S", &tstruct);

    date_str = convertToString(date_buf, 10);
    time_str = convertToString(time_buf, 8);

    return;
}

inline bool file_exist (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}
inline string make_directory (vector<string> args){
    string outputdir = "./results/";
    for (int i = 0 ; i < (int)args.size() ; i++){
        outputdir += args[i] + "/";
        mkdir(outputdir.c_str(), 0776);
    }
    return outputdir;
}
#endif