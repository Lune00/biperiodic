#include"Interactions.hpp"

using namespace std;

Interactions::Interactions(){

	dv_ = 2. ;
	dsv_ = 3. ;
	nv_ = 1 ;
	nsv_ = 1;


}

Interactions::~Interactions(){

}

void Interactions::plug(Sample& spl){
	spl_ = & spl;
}

void Interactions::updatevlist(){


}
