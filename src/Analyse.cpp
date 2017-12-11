#include"Analyse.hpp"
#include"Cell.hpp"
#include"Sample.hpp"


using namespace std;

Analyse::Analyse(){}


Analyse::~Analyse(){}


void Analyse::init(ifstream& is){



}

void Analyse::plug(Sample& spl, Cell& cell){
	spl_ = &spl;
	cell_ = &cell;
}

//Call for different analyses asked by the user
void Analyse::analyse(int tic, double t){
	cout<<"tic = "<<tic<<" - analyse..."<<endl;
	printSample(tic);
}

//Print sample coordonées absolues avec une couche d'epaisseur e
//de particules periodiques
void Analyse::printSample(int tic){

	//Epaisseur couche particules en coordonées reduites
	double e = 1.1 ;
	string localfolder = folder_+"/frames";
	string makepath = "mkdir -p " +localfolder;
	system(makepath.c_str());
	string frame = formatfile(localfolder,"frame.ps",tic);
	//Get vector of image particles within the range e near boundaries
	std::vector<Particle> images = spl_->getimages(e);
	cout<<"Nombre de particules images:"<<images.size()<<endl;
	//Write PS file: particles + images
	writePS(frame,images);

	return ;
}

void Analyse::writePS(const string frame, const vector<Particle>& images){




}
