#include<iostream>//cout,cin
using namespace std;
double periodic(double xi,double xj,double lo,double hi){
	double rij=xi-xj;
	if(rij>(hi-lo)/2){
		rij=hi-lo-rij;
	}else if((-rij)>(hi-lo)/2){
		rij=lo-hi-rij;
	}return rij;
}
