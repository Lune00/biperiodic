#include"Analyse.hpp"
#include"Cell.hpp"
#include"Sample.hpp"
#include"Interactions.hpp"

using namespace std;

//TODO: path to analysis and choices
Analyse::Analyse(){

	allFalse();
	//Thickness of periodic band around sample in printSample
	e_ = 4.;
}

void Analyse::allFalse(){
	printSample_ = false;
	strain_ = false;
	stress_ = false;
	compacity_ = false;
	SP_ = false;
}



Analyse::~Analyse(){}

//TODO
void Analyse::init(ifstream& is){

	string token;
	is >> token;
	while(is){

		if(token=="printSample"){
			printSample_ = true ;
		}

		if(token=="strain"){
			strain_ = true ;
		}

		if(token=="stress"){
			stress_ = true ;
		}

		if(token=="energy"){
			energy_ = true ;
		}
		if(token=="compacity"){
			compacity_ = true;
		}
		if(token=="SP"){
			SP_ = true;
			is >> nbinsSP_ ;
		}
		if(token=="interp"){
			interpenetration_ = true;
		}

		if(token=="}") break;

		is >> token;

	}

	//File managment (empty files used for new analysis)
	cleanFiles();
}


void Analyse::cleanFiles(){

	if(energy_){
		string filename = folder_ + "/energy.txt";
		ofstream o(filename.c_str());
		o.close();
	}

	if(strain_){
		string filename = folder_ + "/strain.txt";
		ofstream o(filename.c_str()); 
		o.close();
	}
	if(stress_){
		string filename = folder_ + "/stress_int.txt";
		string filename2 = folder_ + "/stress_ext.txt";
		ofstream o(filename.c_str());
		o.close();
		o.open(filename2.c_str());
		o.close();
	}
	if(compacity_){
		string filename = folder_ + "/compacity.txt";
		ofstream o(filename.c_str());
		o.close();
	}
	if(SP_){

		string filename = folder_ + "/SProfile.txt";
		ofstream o(filename.c_str());
		o.close();
	}
	if(interpenetration_){
		string filename = folder_ + "/interpenetration.txt";
		ofstream o(filename.c_str());
		o.close();
	}

}

void Analyse::plug(Sample& spl, Cell& cell,Interactions& Int){
	spl_ = &spl;
	cell_ = &cell;
	Int_ = &Int;
}




//Call for different analyses asked by the user
void Analyse::analyse(int tic, double t){
	if(printSample_) printSample(tic);
	if(energy_) computeEnergy(t);
	if(strain_) strain(t);
	if(stress_) stress(t);
	if(compacity_) compacity(t);
	if(SP_) ProfileVelocity(t);
	if(interpenetration_) Interpenetration(t);
}


void Analyse::Interpenetration(const double t) const{
	vector<double> av_max = Int_->getAverageMaxPenetration();
	string file = folder_ + "/interpenetration.txt";
	ofstream os(file.c_str(),ios::app);
	os<<t<<" "<<av_max[0]<<" "<<av_max[1]<<endl;
	os.close();
}



//On utilise les coordonees reduites
//meme si c'est pas une bonne ref car ca depend de la cellule
//On suppose que hyy ne change pas beaucoup...
//Sinon il faut travailler en coordonees absolues et h0(cel ref)
void Analyse::ProfileVelocity(const double t)const {

	double Ly = 1.;
	double ampProbe = Ly / (double)nbinsSP_;
	vector<Probe*> lprobe (nbinsSP_);
	vector<double > Xprofile(nbinsSP_,0.);
	vector<unsigned int> Nbod(nbinsSP_,0);


	Tensor2x2 hd = cell_->gethd();
	Tensor2x2 h = cell_->geth();
	Tensor2x2 hinv = h.getInverse();


	//init probes
	for( unsigned int i = 0 ; i < nbinsSP_ ; i++){
		double yc = 0.5 * ampProbe + i * ampProbe;
		lprobe[i] = new Probe (yc,ampProbe);
	}

	for(vector<Particle>::const_iterator it = spl_->inspectSample().begin(); it != spl_->inspectSample().end(); it++){

		unsigned int j = 0 ;
		while( j < nbinsSP_){
			if(lprobe[j]->containCenter(*it))
			{

				//To do in global function: get v from sd
				//Vecteur v = 
				Xprofile[j] += absVelocity(*it,h,hinv,hd).getx();
				Nbod[j]++;

			}
			if(lprobe[j]->intersection(*it))
			{


			}
			else{
				j++;
			}


		}
	}

	//Writing & Normalise profile:
	string file = folder_ + "/SProfile.txt";
	ofstream os(file.c_str(),ios::app);

	for(unsigned int i = 0; i<nbinsSP_;i++){
		if(Nbod[i]==0) Nbod[i]=1;
		Xprofile[i] /= (double)(Nbod[i]);
		os<<t<<" "<< lprobe[i]->gety()<<" "<<Xprofile[i]<<endl;
	}

	os.close();

	return ;
}

void Analyse::compacity(const double t) const{

	string file = folder_ + "/compacity.txt";
	ofstream os(file.c_str(),ios::app);
	double Vp = 0. ;
	for(vector<Particle>::const_iterator it = spl_->inspectSample().begin(); it != spl_->inspectSample().end(); it++){

		Vp += it->getVolume();
	}
	double sf = Vp / cell_->getVolume();
	os<<t<<" "<<sf<<endl;
	os.close();
}

//TMP
//Add stress_int directement au meme fichier ici
// 3:($7/$9)
void Analyse::strain(const double t) const{
	const Tensor2x2 stress_int = Int_->getStressInt();
	const Tensor2x2 strain = cell_->getStrainTensor();
	string strain_file = folder_ + "/strain.txt";
	ofstream os(strain_file.c_str(),ios::app);
	os<<t<<" ";
	strain.write(os);
	os<<" ";
	stress_int.write(os);
	os<<endl;
	os.close();
}

void Analyse::stress(const double t) const{

	const Tensor2x2 stress_int = Int_->getStressInt();
	const Tensor2x2 stress_ext = cell_->getStressExt();
	string stressI_file = folder_ + "/stress_int.txt";
	string stressE_file = folder_ + "/stress_ext.txt";
	ofstream os(stressI_file.c_str(),ios::app);
	ofstream os2(stressE_file.c_str(),ios::app);
	os<<t<<" ";
	stress_int.write(os);
	os<<endl;
	os2<<t<<" ";
	stress_ext.write(os2);
	os2<<endl;
	os.close();
	os2.close();
}

//Add injected energy, friction and dissipation in inelastic collision
void Analyse::computeEnergy(const double t) {

	string file_energy = folder_ + "/energy.txt";
	ofstream os(file_energy.c_str(),ios::app);

	double TKE = spl_->getTKE();
	double RKE = spl_->getRKE();
	double Elastic = Int_->getElasticEnergy();
	double TotalE = TKE + RKE + Elastic;

	os <<t<<" "<<TKE<<" "<<RKE<<" "<<Elastic<<" "<<TotalE<<endl;

	os.close();
}

//Print sample in absolute coordinates, centered on the center of
//the cell, with a thickness e (in Rmax) of periodic particles around
void Analyse::printSample(int tic){

	//Epaisseur couche particules en rayon max
	string localfolder = folder_+"/frames";
	string makepath = "mkdir -p " +localfolder;
	system(makepath.c_str());
	string frame = formatfile(localfolder,"frame.ps",tic);
	//Get vector of image particles within the range e near boundaries
	std::vector<Particle> images = spl_->getimages(e_);
	//Write PS file: particles + images
	writePS(frame,images);

	return ;
}

//Representation of the sample in absolute coordinates with
//images particle, and without ModularTransformation
//With ModularTransformation (constant box) later... (but to do)
void Analyse::writePS(const string frame, const vector<Particle>& images){


	bool label = true ;
	bool forcenetwork = true ;

	ofstream ps(frame.c_str());
	//Is h used??
	Tensor2x2 h = cell_->geth();
	//Bounding box based on e_ and original cell geometry:
	//Warning: the cell can get out of bounds
	double scaledotradius = 0.05 ;
	double margin = e_ * spl_->getrmax() ;

	//Not yet used
	double zoom = 1. ;

	double xc = 0.;
	double yc = 0. ;

	double xcframe = cell_->getxc();
	double ycframe = cell_->getyc();
	double lx2 = cell_->get_width() * 0.5 ;
	double ly2 = cell_->get_height() * 0.5 ;


	double scalefont = 0.3 * spl_->getrmax();

	double rmax = spl_->getrmax();

	ps<<"%!PS-Adobe-3.0 EPSF-3.0"<<endl;
	ps<<"%%BoundingBox:"<<" "<<xcframe-lx2-margin<<" "<<ycframe-ly2-margin<<" "<<xcframe + lx2+margin<<" "<<ycframe + ly2+margin<<endl;
	ps<<"%%Pages:1"<<endl;
	ps<<"0.1 setlinewidth 0. setgray "<<endl;
	ps <<"0. 0. .23 setrgbcolor clippath fill"<<endl;
	ps << "/colordisk {0.3 0.7 1.0} def"<< endl;
	ps << "/colordot {0. 0. 0.} def" <<endl;
	ps << "/colorimage {1. 1. 1.} def" <<endl;
	ps<<"/Times-Roman findfont"<<endl;
	ps<<scalefont<<" scalefont"<<endl;
	ps<<"setfont"<<endl;

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


		//If label (Id particles)
		if(label){
			ps<<"newpath"<<endl;
			ps<<x<<" "<<y<<" moveto"<<endl;
			ps<<"("<<it->getId()<<") true charpath"<<endl;
			ps<<"0.5 setlinewidth"<<endl;
			ps<<"1. 0. 0. setrgbcolor"<<endl;
			ps<<"fill"<<endl;
			ps<<"stroke"<<endl;
		}
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


	//Draw network
	if(forcenetwork){

		int N = Int_->getnc();

		if(N==0) return ;

		double Fmean = 0. ;
		double Fmin = Int_->inspectContact(0)->getfn(); 
		double Fmax = Int_->inspectContact(0)->getfn(); 

		for(int i = 0 ; i < N ; i++){
			const Contact * c = Int_->inspectContact(i);
			Fmean += c->getfn() ;
			Fmax = max(c->getfn(),Fmax);
			Fmin = min(c->getfn(),Fmin);
		}

		Fmean /= (double)N;
		//cerr<<"fmean = "<<Fmean<<" fmin = "<<Fmin<<" fmax = "<<Fmax<<endl;

		for(int i = 0 ; i < N ; i++){
			const Contact * c = Int_->inspectContact(i);
			if(!c->isActif()) continue;
			double fn = c->getfn();
			double fnres = (fn - Fmin)/(Fmax+Fmin);
			double lw = (fn / Fmean) * 0.1 * rmax ; 
			//double lw = 0.2 * rmax ;

			Vecteur ri = h * c->geti()->getR();
			Vecteur rj = h * c->getj()->getR();

			if(fnres != fnres){
			  cerr<<"@nan : Fmin = "<<Fmin<<" Fmax = "<<Fmax<<" fn = "<<fn<<endl;
			  continue;
			}

			double xi = ri.getx();
			double yi = ri.gety();
			double xj = xi + c->getbranch().getx();
			double yj = yi + c->getbranch().gety();
			ps<<"/coul_force {1 setlinecap 1 "<<1. - fnres<<" "<<1. -fnres<<" setrgbcolor} def"<<endl;
			ps<<lw<<" setlinewidth coul_force"<<endl;
			ps<<xi<<" "<<yi<<" moveto "<<xj<<" "<<yj<<" lineto stroke"<<endl;
		}
	}

	return ;

}

