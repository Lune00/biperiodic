#ifndef Config_hpp
#define Config_hpp

#include<fstream>
#include<iostream>
#include<string>


class Cell;
class Algo;
class Sample;
class Interactions;

class Config{
	private:

		//Paths to write outputs:
		std::string folder_spl_;
		std::string folder_cell_;
		std::string folder_Interactions_;
		std::string folder_analyse_;

	public:
		 Config();
		 ~Config();
		 //Interface d'initialistion
		 //Return 1 pb, return 0 ok
		 //void plug(Algo&,Cell&,Sample&);
		 int init(std::ifstream&,Algo&,Cell&,Sample&,Interactions&);
};

#endif
