#ifndef Config_hpp
#define Config_hpp

#include<fstream>
#include<iostream>
#include<string>

#include"Algo.hpp"
#include"Cell.hpp"
#include"Sample.hpp"
#include"Tenseur.hpp"

class Config{
	private:

		std::string folder_spl_;
		std::string folder_cell_;
		std::string folder_analyse_;

		Sample* sample_;
		Algo* algo_;
		Cell* cell_;

	public:
		 Config();
		 ~Config();
		 //Interface d'initialistion
		 //Return 1 pb, return 0 ok
		 void plug(Algo&,Cell&,Sample&);
		 int init(std::ifstream&,Algo&,Cell&,Sample&);
};

#endif
