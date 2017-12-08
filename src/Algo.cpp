#include"Algo.hpp"
#include"Cell.hpp"
#include"Sample.hpp"
#include"Interactions.hpp"

using namespace std;

void Algo::init(ifstream& is){
	string token;
	is >> token;
	while(is){
		if(token=="dt") is >> dt_;
		if(token=="ns") is >> ns_;
		if(token=="nrecord") is >> nrecord_;
		if(token=="}") break;
		is >> token;
	}
}

//A adapter comme on souhaite en terme de tests
bool Algo::initcheck(){
	if( t_ > 0. ) return false;
	if (dt_ < 1. && ns_ != 0 && nrecord_ != 0 ) return true;
	else return false;
}

void Algo::plug(Cell& cell, Sample& spl, Interactions& Int){
	Int_ = &Int;
	cell_ = &cell;
	spl_ = &spl;
}


void Algo::run(){

	double tfinal = ns_ * dt_ ;

	cout<<"Simulation:"<<endl;
	cout<<"dt = "<<dt_<<endl;
	cout<<"ns = "<<ns_<<endl;
	cout<<"tfinal = "<<tfinal<<endl;

	//Be precautious, restart clock for a new run:
	t_ = 0. ;
	tic_ = 0;
	ticw_ = 0 ;

	while(t_<tfinal){
		//Update verlet list
		Int_->updateverlet(tic_);

		//Time step: integration & periodicity
		verletalgo2();

		//Temp: Analyse, writing:
		if( tic_ % nrecord_ == 0) write();

		t_+=dt_;
		tic_++;
	}
}


//Pour l'instant on l'implemente de maniere naive
//et non optimale (ecriture condensee a l'aide tenseurs/vecteurs)
//On verra apres comment rendre ca plus compacte
void Algo::verletalgo2(){

  double dt2_2 = 0.5 * dt_ * dt_ ;
  double dt_2 = 0.5 * dt_ ;

  vector<Particle>* ps = spl_->getSample();

  // ------------- FIRST STEP VERLET ALGO STARTS HERE

  for(std::vector<Particle>::iterator it = ps->begin(); it != ps->end(); it++){

	  Vecteur a = it->geta();
	  Vecteur r = it->getR();
	  Vecteur v = it->getV();
	  //Rotations, to add
	  r = r + v * dt_ + a * dt2_2;
	  v = v + a * dt_2;
	  //update position et acceleration debut pas temps
	  it->setRV(r,v);

  }

  //Periodicite en position des particules
  //Peut etre a bouger dans Sample plutot
  //Ca me parait plus etre un taff de sample de modifier les positions
  cell_->PeriodicBoundaries2(ps);

  //Integration du mvt de la cellule:
  Tensor2x2 h = cell_->geth();
  Tensor2x2 hd = cell_->gethd();
  Tensor2x2 hdd = cell_->gethdd();

  //h = h + hd * dt_ + hdd * dt2_2 ;
  //Controle en force ou controle en vitesse
  //Temporaire:
  double hxx, hxy , hyx , hyy ;
  double hdxx, hdxy , hdyx , hdyy ;

  if(cell_->getControlxx() == 'f' ) {
	  hxx = h.getxx() + hd.getxx() * dt_ + hdd.getxx() * dt2_2 ;
	  hdxx = hd.getxx() + hdd.getxx() * dt_2 ;
  }
  else{
	  //Vitesse imposee reste la meme (hdd, definit par Ld au debut)
	  //hdd par definition nulle si control en vitesse sur hdd
	  hxx = h.getxx() + hd.getxx() * dt_ ; 
	  hdxx = hd.getxx();
  }

  if(cell_->getControlxy() == 'f' ) {
	  hxy = h.getxy() + hd.getxy() * dt_  + hdd.getxy() * dt2_2 ;
	  hdxy = hd.getxy() + hdd.getxy() * dt_2 ;
  }
  else{
	  hxy = h.getxy() + hd.getxy() * dt_ ; 
	  hdxy = hd.getxy();
  }

  if(cell_->getControlyx() == 'f' ) {
	  hyx = h.getyx() + hd.getyx() * dt_  + hdd.getyx() * dt2_2 ;
	  hdyx = hd.getyx() + hdd.getyx() * dt_2 ;
  }
  else{
	  hyx = h.getyx() + hd.getyx() * dt_ ; 
	  hdyx = hd.getyx();
  }

  if(cell_->getControlyy() == 'f' ) {
	  hyy = h.getyy() + hd.getyy() * dt_  + hdd.getyy() * dt2_2 ;
	  hdyy = hd.getyy() + hdd.getyy() * dt_2 ;
  }
  else{
	  hyy = h.getyy() + hd.getyy() * dt_ ; 
	  hdyy = hd.getyy();
  }

  //Set new h
  h.set(hxx,hxy,hyx,hyy);
  hd.set(hdxx,hdxy,hdyx,hdyy);
  cell_->update(h,hd);


  // ------------- FIRST STEP VERLET ALGO END HERE

  Int_->detectContacts();

  //Calcul des forces entre particules a la nouvelle position fin du pas de temps
  Int_->computeForces();

  //Calcul du tenseur de contraintes internes: sigma_int

  //------------- SECOND STEP VERLET ALGO STARTS HERE
  //Calcul des vitesses a la fin du pas de temps:

  //On retransforme la vitesse en coordonnee reduite
  Tensor2x2 hinv = h.getInverse();
  for(std::vector<Particle>::iterator it = ps->begin(); it != ps->end(); it++){
	  Vecteur a = it->geta();
	  Vecteur v = it->getV();
	  a = hinv * a ;
	  v = v + a * dt_2 ;
	  it->setV(v);
  }

  //Apply stress_ext: si controle en force
  //stress ext est egal a celui impose
  //Sinon il est egal a -stressint et acc nulle
  cell_->ApplyBC();
  double V = cell_->getVolume();
  double mh = cell_->getMasse();
  Tensor2x2 TotalStress = cell_->getStressInt() + cell_->getStressExt();
  hdd = hinv * (V/mh) * (TotalStress);
  //On a l'accleration en fin de pas, on peut integrer l'espace en vitesse a la fin du pas
  hd = hd + hdd * dt_2 ;
  cell_->updatehd(hd);

  //Cell deformation
  cell_->CalculStrainTensor();
}


void Algo::write(){
	cout<<"Write outputs"<<endl;
	spl_->writeAbsolute(ticw_);
	Int_->writeContacts(ticw_);
	//cell_->write(testcell,t);
	//cell_->writeStrainTensor(strain,t);
	ticw_++;
}
