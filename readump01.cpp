#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#include<iostream>//string
#include<vector>//vector
#include<fstream>//ifstream
#include<sstream>//istringstream
using namespace std;
void readump01(string file_name,vector<int> &steps,int &N,vector<double> &xlo,vector<double> &xhi,vector<double> &ylo,vector<double> &yhi,vector<double> &zlo,vector<double> &zhi,vector<vector<int> > &id,vector<vector<int> > &mol,vector<vector<int> > &type,vector<vector<double> > &q,vector<vector<double> > &x,vector<vector<double> > &y,vector<vector<double> > &z,vector<vector<double> > &ix,vector<vector<double> > &iy,vector<vector<double> > &iz,vector<vector<double> > &masses,vector<vector<double> > &vx,vector<vector<double> > &vy,vector<vector<double> > &vz,vector<vector<double> > &fx,vector<vector<double> > &fy,vector<vector<double> > &fz){
	steps.clear();
	xlo.clear();
	xhi.clear();
	ylo.clear();
	yhi.clear();
	zlo.clear();
	zhi.clear();
	id.clear();
	mol.clear();
	type.clear();
	q.clear();
	x.clear();
	y.clear();
	z.clear();
	ix.clear();
	iy.clear();
	iz.clear();
	masses.clear();
	vx.clear();
	vy.clear();
	vz.clear();
	fx.clear();
	fy.clear();
	fz.clear();
	int step,tag;
	double Xlo,Xhi,Ylo,Yhi,Zlo,Zhi;
	vector<int> atm_id,mol_id,atm_type;
	vector<double> Q,X,Y,Z,Ix,Iy,Iz,mass,Vx,Vy,Vz,Fx,Fy,Fz;
	ifstream dump;
    dump.open(file_name.c_str());
    string line;
    while(getline(dump,line)){
    	istringstream iss(line);
    	string item;
    	iss>>item;
    	string sub;
    	iss>>sub;
    	if(sub=="TIMESTEP"){
    		getline(dump,line);
    		sscanf(line.c_str(),"%i",&step);
    		steps.push_back(step);
    	}else if(sub=="NUMBER"){
    		getline(dump,line);
    		sscanf(line.c_str(),"%i",&N);
    		atm_id.resize(N);
    		mol_id.resize(N);
    		atm_type.resize(N);
    		Q.resize(N);
    		X.resize(N);
    		Y.resize(N);
    		Z.resize(N);
			Ix.resize(N);
			Iy.resize(N);
			Iz.resize(N);
    		mass.resize(N);
    		Vx.resize(N);
    		Vy.resize(N);
    		Vz.resize(N);
			Fx.resize(N);
			Fy.resize(N);
			Fz.resize(N);
    	}else if(sub=="BOX"){
    		getline(dump,line);
    		sscanf(line.c_str(),"%lf %lf",&Xlo,&Xhi);
			getline(dump,line);
    		sscanf(line.c_str(),"%lf %lf",&Ylo,&Yhi);
			getline(dump,line);
    		sscanf(line.c_str(),"%lf %lf",&Zlo,&Zhi);
    		xlo.push_back(Xlo);
    		xhi.push_back(Xhi);
    		ylo.push_back(Ylo);
    		yhi.push_back(Yhi);
    		zlo.push_back(Zlo);
    		zhi.push_back(Zhi);
    	}else if(sub=="ATOMS"){
    		for(int n=0;n<N;n++){
    			getline(dump,line);
    			sscanf(line.c_str(),"%i %*i %*i %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f",&tag);
    			sscanf(line.c_str(),"%i %i %i %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&atm_id[tag-1],&mol_id[tag-1],&atm_type[tag-1],&Q[tag-1],&X[tag-1],&Y[tag-1],&Z[tag-1],&Ix[tag-1],&Iy[tag-1],&Iz[tag-1],&mass[tag-1],&Vx[tag-1],&Vy[tag-1],&Vz[tag-1],&Fx[tag-1],&Fy[tag-1],&Fz[tag-1]);
    		}id.push_back(atm_id);
    		mol.push_back(mol_id);
    		type.push_back(atm_type);
    		q.push_back(Q);
    		x.push_back(X);
    		y.push_back(Y);
    		z.push_back(Z);
			ix.push_back(Ix);
			iy.push_back(Iy);
			iz.push_back(Iz);
    		masses.push_back(mass);
    		vx.push_back(Vx);
    		vy.push_back(Vy);
    		vz.push_back(Vz);
			fx.push_back(Fx);
			fy.push_back(Fy);
			fz.push_back(Fz);
    	}
    }dump.close();
}
