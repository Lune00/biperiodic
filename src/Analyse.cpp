#include"Analyse.hpp"
#include"Cell.hpp"
#include"Sample.hpp"


using namespace std;

Analyse::Analyse(){ e_ = 4.1;}


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

	//Epaisseur couche particules en rayon max
	string localfolder = folder_+"/frames";
	string makepath = "mkdir -p " +localfolder;
	system(makepath.c_str());
	string frame = formatfile(localfolder,"frame.ps",tic);
	//Get vector of image particles within the range e near boundaries
	std::vector<Particle> images = spl_->getimages(e_);
	cout<<"Nombre de particules images:"<<images.size()<<endl;
	//Write PS file: particles + images
	writePS(frame,images);

	return ;
}

//Representation of the samle in absolute coordinates with
//images particle, and without ModularTransformation
//show cell deformation
//With ModularTransformation (constant box) later... (but to do)
void Analyse::writePS(const string frame, const vector<Particle>& images){

	ofstream ps(frame.c_str());
	Tensor2x2 h = cell_->geth();
	//Bounding box based on e_ and original cell geometry:
	//Warning: the cell can get out of bounds

	double scaledotradius = 0.05 ;
	double margin = 2. * e_ ;
	double zoom = 1. ;

	double xmin = 0. ;
	double xmax = cell_->getLx() + margin ;
	double ymin = 0.  ;
	double ymax = cell_->getLy() + margin;


//	double xc = xmin + width * 0.5 - spl_->getrmax() ;
//	double yc = ymin + height * 0.5 ;
	double xc = cell_->getxc();
	double yc = cell_->getyc();

	double width = 2. * (xmax - xc) - margin + spl_->getrmax();
	double height = 2. * (ymax - yc) - margin;


	ps<<"%!PS-Adobe-3.0 EPSF-3.0"<<endl;
	ps<<"%%BoundingBox:"<<" "<<xmin<<" "<<ymin<<" "<<width<<" "<<height<<endl;
	ps<<"%%Pages:1"<<endl;
	ps<<"0.1 setlinewidth 0. setgray "<<endl;
	ps <<"0. 0. .23 setrgbcolor clippath fill"<<endl;
	ps << "/colordisk {0.3 0.7 1.0} def"<< endl;
	ps << "/colordot {0. 0. 0.} def" <<endl;
	ps << "/colorimage {1. 1. 1.} def" <<endl;
	//Draw sample
	for(vector<Particle>::const_iterator it = spl_->inspectSample().begin(); it != spl_->inspectSample().end(); it++){

		Vecteur rabs = h * it->getR();
		double x = xc + rabs.getx();
		double y = yc + rabs.gety();
		double r = it->getRadius();
		double theta = it->getRot();
		double xrcostheta = x + r * cos(theta) * 0.8 ;
		double yrcostheta = y + r * sin(theta) * 0.8 ;
		double radiusrot = r * scaledotradius;

		ps <<" newpath "<<endl ;
		ps <<x<<" "<<y<<" "<<r<<" colordisk setrgbcolor 0.0 setlinewidth 0 360 arc gsave fill grestore "<<endl; 
		ps <<"stroke"<<endl;
		ps << "newpath "<<endl;
		ps <<xrcostheta<<" "<<yrcostheta<<" "<<radiusrot<<" colordot setrgbcolor 0.0 setlinewidth 0 360 arc gsave fill grestore"<<endl;
		ps<<"stroke"<<endl;
	}

	//Draw image particles
	for(vector<Particle>::const_iterator it = images.begin(); it != images.end(); it++){
		Vecteur rabs = h * it->getR();
		double x = xc + rabs.getx();
		double y = yc + rabs.gety();
		double r = it->getRadius();
		double theta = it->getRot();
		double xrcostheta = x + r * cos(theta) * 0.8 ;
		double yrcostheta = y + r * sin(theta) * 0.8 ;
		double radiusrot = r * scaledotradius;

		ps <<" newpath "<<endl ;
		ps <<x<<" "<<y<<" "<<r<<" colorimage setrgbcolor 0.0 setlinewidth 0 360 arc gsave fill grestore "<<endl; 
		ps <<"stroke"<<endl;
		ps << "newpath "<<endl;
		ps <<xrcostheta<<" "<<yrcostheta<<" "<<radiusrot<<" colordot setrgbcolor 0.0 setlinewidth 0 360 arc gsave fill grestore"<<endl;
		ps<<"stroke"<<endl;
	}

	//Draw periodic cell

	double ux=h.getxx();
	double uy=h.getyx();
	double vx=h.getxy();
	double vy=h.getyy();

	double lw = 0.1 * spl_->getrmax();

	//a1 base vector
	ps<<"newpath"<<endl;
	ps<<xc<<" "<<yc<<" moveto " <<endl;
	ps<<ux<<" "<<uy<<" rlineto"<<endl;
	ps<<vx<<" "<<vy<<" rlineto"<<endl;
	ps<<-ux<<" "<<-uy<<" rlineto"<<endl;
	ps<<-vx<<" "<<-vy<<" rlineto"<<endl;
	ps<<"closepath"<<endl;
	//yellow 1, 0.937, 0.078
	//light gray 0.819 0.819 0.8191
	//green 0.435, 1, 0.078
	//darker gray 0.549, 0.549, 0.549
	ps<<"0.549 0.549 0.549 setrgbcolor"<<endl;
	ps<<lw<<" setlinewidth "<<endl;
	ps<<"stroke"<<endl;


	return ;

}
