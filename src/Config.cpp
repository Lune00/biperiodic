#include"Config.hpp"

using namespace std;

//Essai mecanique: 2 fichiers input parametres pour Config, echantillon initial (xmin,xmax...)
//Prendra input eventuellement dans un fichier
Config::Config(){
	//Pour le moment tout en vitesse (cinematique)
	for(int i=0;i<4;i++){
		BCU[i] = 'v';
	}
	//BCU[3]='f';
	//Par souci de convention, les DL non controles sont init a 0
	//Par definition, le taux de cisaillement pur est egale a la moitie de LdUser.yx(ou xy)
	LdUser_.set(1.,0.,0.,0.);
	StressUser_.set(0.,0.,0.,0.);
}
