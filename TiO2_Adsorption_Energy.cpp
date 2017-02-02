#include<iostream>//cout
#include<string>//string
#include<vector>//vector
#define _USE_MATH_DEFINES
#include<math.h>//fabs
#include<fstream>//ofstream
#include<algorithm>//find
using namespace std;

void readump01(string file_name,vector<int> &steps,int &N,vector<double> &xlo,vector<double> &xhi,vector<double> &ylo,vector<double> &yhi,vector<double> &zlo,vector<double> &zhi,vector<vector<int> > &id,vector<vector<int> > &mol,vector<vector<int> > &type,vector<vector<double> > &q,vector<vector<double> > &x,vector<vector<double> > &y,vector<vector<double> > &z,vector<vector<double> > &ix,vector<vector<double> > &iy,vector<vector<double> > &iz,vector<vector<double> > &masses,vector<vector<double> > &vx,vector<vector<double> > &vy,vector<vector<double> > &vz,vector<vector<double> > &fx,vector<vector<double> > &fy,vector<vector<double> > &fz);
double periodic(double init,double fin,double lo,double hi);

int main(){
	int N,files;
	vector<string> file_names;
    vector<int> steps,tmp_steps;
    vector<double> xlo,xhi,ylo,yhi,zlo,zhi,tmp_xlo,tmp_xhi,tmp_ylo,tmp_yhi,tmp_zlo,tmp_zhi;
    vector<vector<int> > id,mol,type,tmp_id,tmp_mol,tmp_type;
    vector<vector<double> > q,x,y,z,ix,iy,iz,masses,vx,vy,vz,fx,fy,fz,tmp_q,tmp_x,tmp_y,tmp_z,tmp_ix,tmp_iy,tmp_iz,tmp_masses,tmp_vx,tmp_vy,tmp_vz,tmp_fx,tmp_fy,tmp_fz;
	int fluid,surfs,closest;
	double cutoff,dx,dy,dz,R,closeR,closedx,closedy,closedz,xHH,yHH,zHH,dipole,theta;
	vector<int> surf_types,fluids,sites;
	vector<double> angZ_avg,ang_avg,angZ_std,ang_std;
	vector<vector<double> > angZ,ang;
	cout<<"number of dump files = ";
	cin>>files;
	file_names.resize(files);
	for(int file=0;file<files;file++){
		cout<<"name of dump file "<<file+1<<": ";
		cin>>file_names[file];
	}cout<<"fluid atom type = ";
	cin>>fluid;
	cout<<"number of surface bead types = ";
	cin>>surfs;
	for(int surf=0;surf<surfs;surf++){
		surf_types.push_back(0);
		cout<<"surface bead type "<<surf+1<<" = ";
		cin>>surf_types[surf];
	}cout<<"cutoff = ";
	cin>>cutoff;
	readump01(file_names[0],steps,N,xlo,xhi,ylo,yhi,zlo,zhi,id,mol,type,q,x,y,z,ix,iy,iz,masses,vx,vy,vz,fx,fy,fz);
	for(int file=1;file<files;file++){
		readump01(file_names[file],tmp_steps,N,tmp_xlo,tmp_xhi,tmp_ylo,tmp_yhi,tmp_zlo,tmp_zhi,tmp_id,tmp_mol,tmp_type,tmp_q,tmp_x,tmp_y,tmp_z,tmp_ix,tmp_iy,tmp_iz,tmp_masses,tmp_vx,tmp_vy,tmp_vz,tmp_fx,tmp_fy,tmp_fz);
		steps.insert(steps.end(),tmp_steps.begin()+1,tmp_steps.end());
		xlo.insert(xlo.end(),tmp_xlo.begin()+1,tmp_xlo.end());
		xhi.insert(xhi.end(),tmp_xhi.begin()+1,tmp_xhi.end());
		ylo.insert(ylo.end(),tmp_ylo.begin()+1,tmp_ylo.end());
		yhi.insert(yhi.end(),tmp_yhi.begin()+1,tmp_yhi.end());
		zlo.insert(zlo.end(),tmp_zlo.begin()+1,tmp_zlo.end());
		zhi.insert(zhi.end(),tmp_zhi.begin()+1,tmp_zhi.end());
		id.insert(id.end(),tmp_id.begin()+1,tmp_id.end());
		mol.insert(mol.end(),tmp_mol.begin()+1,tmp_mol.end());
		type.insert(type.end(),tmp_type.begin()+1,tmp_type.end());
		q.insert(q.end(),tmp_q.begin()+1,tmp_q.end());
		x.insert(x.end(),tmp_x.begin()+1,tmp_x.end());
		y.insert(y.end(),tmp_y.begin()+1,tmp_y.end());
		z.insert(z.end(),tmp_z.begin()+1,tmp_z.end());
		ix.insert(ix.end(),tmp_ix.begin()+1,tmp_ix.end());
		iy.insert(iy.end(),tmp_iy.begin()+1,tmp_iy.end());
		iz.insert(iz.end(),tmp_iz.begin()+1,tmp_iz.end());
		masses.insert(masses.end(),tmp_masses.begin()+1,tmp_masses.end());
		vx.insert(vx.end(),tmp_vx.begin()+1,tmp_vx.end());
		vy.insert(vy.end(),tmp_vy.begin()+1,tmp_vy.end());
		vz.insert(vz.end(),tmp_vz.begin()+1,tmp_vz.end());
		fx.insert(fx.end(),tmp_fx.begin()+1,tmp_fx.end());
		fy.insert(fy.end(),tmp_fy.begin()+1,tmp_fy.end());
		fz.insert(fz.end(),tmp_fz.begin()+1,tmp_fz.end());
	}for(int n=0;n<N;n++){
		if(type[0][n]==fluid){
			fluids.push_back(id[0][n]);
		}else if(find(surf_types.begin(),surf_types.end(),type[0][n])!=surf_types.end()){
			sites.push_back(id[0][n]);
		}
	}



           for(int t=0;t<steps.size();t++)
                {
		angZ.push_back(vector<double>(181,0));
		ang.push_back(vector<double>(181,0));
		for(int i=0;i<fluids.size();i++){
			closeR=xhi[t]-xlo[t]+yhi[t]-ylo[t]+zhi[t]-zlo[t];
			for(int j=0;j<sites.size();j++){
				dx=periodic(x[t][sites[j]-1],x[t][fluids[i]-1],xlo[t],xhi[t]);
				dy=periodic(y[t][sites[j]-1],y[t][fluids[i]-1],ylo[t],yhi[t]);
				dz=periodic(z[t][sites[j]-1],z[t][fluids[i]-1],zlo[t],zhi[t]);
				R=pow((dx*dx)+(dy*dy)+(dz*dz),.5);
				if(R<closeR){
					closedx=dx;
					closedy=dy;
					closedz=dz;
					closeR=R;
					closest=sites[j];
                                        cout << sites[j] << endl;
				}
		}


                        xHH=periodic(x[t][fluids[i]-1],x[t][fluids[i]]+(periodic(x[t][fluids[i]],x[t][fluids[i]+1],xlo[t],xhi[t])/2),xlo[t],xhi[t]);
			
                        yHH=periodic(y[t][fluids[i]-1],y[t][fluids[i]]+(periodic(y[t][fluids[i]],y[t][fluids[i]+1],ylo[t],yhi[t])/2),ylo[t],yhi[t]);
			 
                        zHH=periodic(z[t][fluids[i]-1],z[t][fluids[i]]+(periodic(z[t][fluids[i]],z[t][fluids[i]+1],zlo[t],zhi[t])/2),zlo[t],zhi[t]);
			
                        dipole=pow((xHH*xHH)+(yHH*yHH)+(zHH*zHH),.5);
			
                        theta=acos(-zHH/dipole)*(180/M_PI);
			
                        angZ.back()[(int)round(theta)]++;
			
                        theta=acos(((-xHH*closedx)+(-yHH*closedy)+(-zHH*closedz))/(closeR*dipole))*(180/M_PI);
			
                        ang.back()[(int)round(theta)]++;
		        
                        }
	}ang_avg.resize(181,0);
	angZ_avg.resize(181,0);
	ang_std.resize(181,0);
	angZ_std.resize(181,0);
	for(int a=0;a<ang_avg.size();a++){
		for(int u=0;u<ang.size();u++){
			ang_avg[a]+=ang[u][a]/ang.size();
			angZ_avg[a]+=angZ[u][a]/angZ.size();
		}for(int u=0;u<ang.size();u++){
			ang_std[a]+=(ang[u][a]-ang_avg[a])*(ang[u][a]-ang_avg[a]);
			angZ_std[a]+=(angZ[u][a]-angZ_avg[a])*(angZ[u][a]-angZ_avg[a]);
		}ang_std[a]=pow(ang_std[a]/ang.size(),.5);
		angZ_std[a]=pow(angZ_std[a]/angZ.size(),.5);
	}ofstream ang_file;
	ang_file.open("angles.txt",ios::trunc);
	ang_file<<"dipole angle"<<"\t"<<"Z (avg)"<<"\t"<<"pair vector (avg)"<<"\t"<<"Z (stdev)"<<"\t"<<"pair vector (stdev)"<<"\n";
	for(int a=0;a<ang_avg.size();a++){
		ang_file<<a<<"\t"<<angZ_avg[a]<<"\t"<<ang_avg[a]<<"\t"<<angZ_std[a]<<"\t"<<ang_std[a]<<"\n";
	}ang_file.close();
	cout<<"press any key & enter to close";
	cin>>file_names[0];
	return 0;
}
