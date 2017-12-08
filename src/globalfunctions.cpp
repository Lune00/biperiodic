#include"globalfunctions.hpp"

using namespace std;

//Format output to current tic step k
string formatfile(const string file, const int k){
	std::stringstream ss;
	ss << std::setw(wi)<< std::setfill('0') << k ;
	string newfile = ss.str() + file;
	return newfile;
}


