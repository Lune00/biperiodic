#include"Sample.hpp"

using namespace std;

Sample::Sample(){

}

Sample::~Sample(){
}


void Sample::init(ifstream& is){

	string token;
	is >> token;
	while(is){
		if(token=="includeFile") is >> fichsample_;
		if(token=="rho") is >> rho_;
		is >> token ;
	}
}

void Sample::loadSample(ifstream& is){

}

void Sample::write(ofstream& os){

}

//Fill particles with mass and inertia
void Sample::attributeMass(){


}

void Sample::setminmax(){



}
