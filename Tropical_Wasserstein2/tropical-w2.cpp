#include <iostream>
#include <fftw3.h>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <time.h>

using namespace std;

class poisson_solver{
public:
    fftw_plan planIn;
    fftw_plan planOut;
    double *fourierWorkspace;
    double *u;
    double *kernel;
    double *tridiagonalWorkspace;

    int n1;
    int n2;
    int nt;
    double dx;
    double dy;
    double dt;

    poisson_solver(int n1, int n2, int nt, double dx, double dy, double dt) {
    	this->n1=n1;
    	this->n2=n2;
    	this->nt=nt;
    	this->dx=dx;
    	this->dy=dy;
    	this->dt=dt;

        fourierWorkspace =(double*) fftw_malloc(n1*n2*sizeof(double));

		planIn = fftw_plan_r2r_2d(n2,n1, fourierWorkspace, fourierWorkspace, FFTW_REDFT10,FFTW_REDFT10, FFTW_MEASURE);
		planOut = fftw_plan_r2r_2d(n2,n1, fourierWorkspace, fourierWorkspace, FFTW_REDFT01,FFTW_REDFT01, FFTW_MEASURE);

		u=new double[n1*n2*nt];
        kernel=new double[n1*n2];
        tridiagonalWorkspace=new double[n1*n2*nt];

        create_negative_laplacian_kernel_2d();
    }

    void create_negative_laplacian_kernel_2d(){
	    int pcount=n1*n2;
	    
	    for(int i=0;i<n2;i++){
	        for(int j=0;j<n1;j++){
	            double negativeLaplacian=2/(dx*dx)*(1-cos(M_PI*(j)*dx)) + 2/(dy*dy)*(1-cos(M_PI*(i)*dy));
	            kernel[i*n1+j]=negativeLaplacian;
	        }
	    }
	}

	void forward_tridiagonal_sweep(){

	    for(int i=0;i<n1*n2;i++){
	        tridiagonalWorkspace[i]=0.0/1.0;
	        u[i]=0;
	    }
	     
	    for(int k=1;k<nt-1;k++){
	        for(int i=0;i<n1*n2;i++){
	            double alpha=kernel[i]/(nt*nt);
	            tridiagonalWorkspace[k*n1*n2+i]=-1/(2+alpha+tridiagonalWorkspace[(k-1)*n1*n2+i]);
	            u[k*n1*n2+i]=(u[k*n1*n2+i]+u[(k-1)*n1*n2+i])/(2+alpha+tridiagonalWorkspace[(k-1)*n1*n2+i]);
	        }
	    }
	    
	    for(int i=0;i<n1*n2;i++){
	        u[(nt-1)*n1*n2+i]=0;
	    }
	    
	}

	void backward_tridiagonal_sweep(){
	    for(int k=nt-2;k>=0;k--){
	        for(int i=0;i<n1*n2;i++){
	            u[k*n1*n2+i]=u[k*n1*n2+i]-tridiagonalWorkspace[k*n1*n2+i]*u[(k+1)*n1*n2+i];
	        }
	        
	    }
	}

	void perform_inverse_laplacian(){

		fftw_execute(planIn);
		fourierWorkspace[0]=0;

		for(int i=1;i<n1*n2;++i){
			fourierWorkspace[i]/=4*(n1)*(n2)*kernel[i];
		}

		fftw_execute(planOut);
	}

	void get_fourier_coefficients(double* u){

		for(int i=0;i<n1*n2;++i){
			fourierWorkspace[i]/=1.0*nt*nt;	
		}

		fftw_execute(planIn);

		for(int i=0;i<n1*n2;++i){
			u[i]=fourierWorkspace[i]/(4.0*n1*n2);
		}
	}
	void back_to_real_space(){
		fftw_execute(planOut);
	}

	void destroy_all_fftps(){
		free(u);
	    free(fourierWorkspace);
	    free(kernel);
	    fftw_destroy_plan(planIn);
	    fftw_destroy_plan(planOut);
	}
};


void create_csv_file_for_parameters(int n1,int n2,int nt){
	ofstream outfile;
	outfile.open("./data/parameters.csv");
	outfile<<n1<<","<<n2<<","<<nt;
	outfile.close();
}

void create_csv_file(const double* A,string filename,int n1,int n2,int nt){
	ofstream outfile;
	outfile.open(filename);
	for(int i=0;i<n1*n2*nt;++i){
		outfile<<A[i]<<",";
	}
	outfile.close();
}

void create_csv_file_rho(const double* rho,const double* mu,const double* nu,string filename,int n1,int n2,int nt){
	ofstream outfile;
	outfile.open(filename);
	for(int i=0;i<n1*n2;++i){
		outfile<<mu[i]<<",";
	}
	for(int i=0;i<n1*n2*(nt-2);++i){
		outfile<<rho[n1*n2+i]<<",";
	}
	for(int i=0;i<n1*n2;++i){
		outfile<<nu[i]<<",";
	}
	outfile.close();
}



void create_initial_densities_one_to_four(double* rho,double base,int n1, int n2, int nt, double dx, double dy, double dt){
	double sum_rho0=0;
	double sum_rho1=0;
	for(int i=0;i<n2;++i){
		for(int j=0;j<n1;++j){
			double x=j*dx;
			double y=i*dy;
			if(abs(x-0.5)<0.1 && abs(y-0.5)<0.1){
				rho[i*n1+j]+=1;
			}

			if(abs(x-0.2)<0.05 && abs(y-0.2)<0.05){
				rho[(nt-1)*n1*n2+i*n1+j]+=1;
			}
			if(abs(x-0.2)<0.05 && abs(y-0.8)<0.05){
				rho[(nt-1)*n1*n2+i*n1+j]+=1;
			}
			if(abs(x-0.8)<0.05 && abs(y-0.2)<0.05){
				rho[(nt-1)*n1*n2+i*n1+j]+=1;
			}
			if(abs(x-0.8)<0.05 && abs(y-0.8)<0.05){
				rho[(nt-1)*n1*n2+i*n1+j]+=1;
			}
			rho[i*n1+j]+=base;
			rho[(nt-1)*n1*n2+i*n1+j]+=base;

			sum_rho0+=rho[i*n1+j];
			sum_rho1+=rho[(nt-1)*n1*n2+i*n1+j];
		}
	}
	for(int i=0;i<n1*n2;++i){
		rho[i]/=sum_rho0*dx*dy;
		rho[(nt-1)*n1*n2+i]/=sum_rho1*dx*dy;
	}

	for(int n=1;n<nt-1;++n){
		for(int i=0;i<n1*n2;++i){
			rho[n*n1*n2+i]=1;	
		}
	}
}


void create_initial_densities_one_to_one_positive(double* rho,double base,int n1, int n2, int nt, double dx, double dy, double dt){
	double sum_rho0=0;
	double sum_rho1=0;
	for(int i=0;i<n2;++i){
		for(int j=0;j<n1;++j){
			double x=j*dx;
			double y=i*dy;
			if(abs(x-1.0/3)<0.1 && abs(y-1.0/3)<0.1){
				rho[i*n1+j]+=1;
			}

			if(abs(x-2.0/3)<0.1 && abs(y-2.0/3)<0.1){
				rho[(nt-1)*n1*n2+i*n1+j]+=1;
			}
			rho[i*n1+j]+=base;
			rho[(nt-1)*n1*n2+i*n1+j]+=base;

			sum_rho0+=rho[i*n1+j];
			sum_rho1+=rho[(nt-1)*n1*n2+i*n1+j];
		}
	}
	for(int i=0;i<n1*n2;++i){
		rho[i]/=sum_rho0*dx*dy;
		rho[(nt-1)*n1*n2+i]/=sum_rho1*dx*dy;
	}

	for(int n=1;n<nt-1;++n){
		for(int i=0;i<n1*n2;++i){
			rho[n*n1*n2+i]=1;	
		}
	}
}

void create_initial_densities_one_to_one_positive_circle(double* rho,double base,int n1, int n2, int nt, double dx, double dy, double dt){
	double sum_rho0=0;
	double sum_rho1=0;
	for(int i=0;i<n2;++i){
		for(int j=0;j<n1;++j){
			double x=j*dx;
			double y=i*dy;
			if(pow(x-0.3,2)+pow(y-0.3,2)<pow(0.1,2)){
				rho[i*n1+j]+=1;
			}

			if(pow(x-0.7,2)+pow(y-0.7,2)<pow(0.1,2)){
				rho[(nt-1)*n1*n2+i*n1+j]+=1;
			}
			rho[i*n1+j]+=base;
			rho[(nt-1)*n1*n2+i*n1+j]+=base;

			sum_rho0+=rho[i*n1+j];
			sum_rho1+=rho[(nt-1)*n1*n2+i*n1+j];
		}
	}
	for(int i=0;i<n1*n2;++i){
		rho[i]/=sum_rho0*dx*dy;
		rho[(nt-1)*n1*n2+i]/=sum_rho1*dx*dy;
	}

	for(int n=1;n<nt-1;++n){
		for(int i=0;i<n1*n2;++i){
			rho[n*n1*n2+i]=1;	
		}
	}
}

void create_initial_densities_one_to_one_negative(double* rho,double base,int n1, int n2, int nt, double dx, double dy, double dt){
	double sum_rho0=0;
	double sum_rho1=0;
	for(int i=0;i<n2;++i){
		for(int j=0;j<n1;++j){
			double x=j*dx;
			double y=i*dy;
			if(abs(x-1.0/3)<0.1 && abs(y-2.0/3)<0.1){
				rho[i*n1+j]+=1;
			}

			if(abs(x-2.0/3)<0.1 && abs(y-1.0/3)<0.1){
				rho[(nt-1)*n1*n2+i*n1+j]+=1;
			}
			rho[i*n1+j]+=base;
			rho[(nt-1)*n1*n2+i*n1+j]+=base;

			sum_rho0+=rho[i*n1+j];
			sum_rho1+=rho[(nt-1)*n1*n2+i*n1+j];
		}
	}
	for(int i=0;i<n1*n2;++i){
		rho[i]/=sum_rho0*dx*dy;
		rho[(nt-1)*n1*n2+i]/=sum_rho1*dx*dy;
	}

	for(int n=1;n<nt-1;++n){
		for(int i=0;i<n1*n2;++i){
			rho[n*n1*n2+i]=1;	
		}
	}
}

double calculate_tropical_norm(double mx, double my){
    if(mx>=my && my>=0){
        return mx;
    }
    else if(my>=mx && mx>=0){
        return my;
    }
    else if(mx>=0 && 0 >= my){
        return mx - my;
    }
    else if(my>=0 && 0 >= mx){
        return my - mx;
    }
    else if(0>=mx && mx >= my){
        return - my;
    }
    else if(0>=my && my >= mx){
        return - mx;
    }

    cout << "wrong answer";

    exit(0);
    return 0;
}


double calculateEnergy(const double* rho, const double* mx, const double* my,int n1, int n2, int nt, double dx, double dy, double dt){
    double sum = 0;

    for(int n=0; n<nt; ++n){
        for(int i=0; i<n2; ++i){
            for(int j=0; j<n1; ++j){
            	double mxvalue=mx[n*n1*n2+i*n1+j];
            	double myvalue=my[n*n1*n2+i*n1+j];
            	double rhovalue=rho[n*n1*n2+i*n1+j];
                double tropNorm=calculate_tropical_norm(mxvalue,myvalue);
                sum += tropNorm * tropNorm / (2.0 * rhovalue);
            }
        }
    }
    return sum * dx * dy * dt;
}

double calculate_dual(const double* Phi, const double* rho, int n1, int n2, int nt, double dx, double dy, double dt){
	double sum = 0;
	for(int i=0;i<n2;++i){
		for(int j=0;j<n1;++j){
			sum += Phi[(nt-1)*n1*n2+i*n1+j]*rho[(nt-1)*n1*n2+i*n1+j] - Phi[i*n1+j]*rho[i*n1+j];
		}
	}
	return sum*dx*dy;
}

double calculate_grad_mx(const double* mxTmp,int n,int i,int j,int n1,int n2,int nt){
	double mxvalue;
	if(j==0){
		mxvalue=1.0*n1*(mxTmp[n*n1*n2+i*n1+j]);
	}else{
		mxvalue=1.0*n1*(mxTmp[n*n1*n2+i*n1+j]-mxTmp[n*n1*n2+i*n1+j-1]);
	}
	return mxvalue;
}

double calculate_grad_my(const double* myTmp,int n,int i,int j,int n1,int n2,int nt){
	double myvalue;
	if(i==0){
		myvalue=1.0*n2*(myTmp[n*n1*n2+i*n1+j]);
	}else{
		myvalue=1.0*n2*(myTmp[n*n1*n2+i*n1+j]-myTmp[n*n1*n2+(i-1)*n1+j]);
	}
	return myvalue;
}
	

void update_Phi(poisson_solver& fftps, double* Phi,const double* rhoTmp,const double* mxTmp,const double* myTmp,double sigma_PDHG,int n1,int n2,int nt,double dx,double dy,double dt){


    // Update t = 0

	for(int i=0;i<n2;++i){
    	for(int j=0;j<n1;++j){
    		int n = 0;
	    	double partial_t_rho=1.0/dt*(rhoTmp[(n+1)*n1*n2+i*n1+j]-rhoTmp[n*n1*n2+i*n1+j]);
			double mxvalue=calculate_grad_mx(mxTmp,n,i,j,n1,n2,nt);
			double myvalue=calculate_grad_my(myTmp,n,i,j,n1,n2,nt);
			fftps.fourierWorkspace[i*n1+j]=(partial_t_rho+mxvalue+myvalue);	
    	}
    }

    fftps.perform_inverse_laplacian();
    for(int i=0;i<n1*n2;++i){
    	Phi[i] += sigma_PDHG*fftps.fourierWorkspace[i];
    }

    // Update t = 1

    for(int i=0;i<n2;++i){
    	for(int j=0;j<n1;++j){
    		int n = nt-1;
	    	double partial_t_rho=1.0/dt*(rhoTmp[n*n1*n2+i*n1+j]-rhoTmp[(n-1)*n1*n2+i*n1+j]);
			double mxvalue=calculate_grad_mx(mxTmp,n,i,j,n1,n2,nt);
			double myvalue=calculate_grad_my(myTmp,n,i,j,n1,n2,nt);
			fftps.fourierWorkspace[i*n1+j]=(partial_t_rho+mxvalue+myvalue);	
    	}
    }

    fftps.perform_inverse_laplacian();
    for(int i=0;i<n1*n2;++i){
    	Phi[(nt-1)*n1*n2+i] += sigma_PDHG*fftps.fourierWorkspace[i];
    }

	// Update 0 < t < 1

	for(int n=0;n<nt;++n){	
		for(int i=0;i<n2;++i){
			for(int j=0;j<n1;++j){
				if(n==0){
					fftps.fourierWorkspace[i*n1+j]=0;
				}else if(n==nt-1){
					fftps.fourierWorkspace[i*n1+j]=0;
				}else{
					double partial_t_rho=0.5/dt*(rhoTmp[(n+1)*n1*n2+i*n1+j]-rhoTmp[(n-1)*n1*n2+i*n1+j]);
					double mxvalue=calculate_grad_mx(mxTmp,n,i,j,n1,n2,nt);
					double myvalue=calculate_grad_my(myTmp,n,i,j,n1,n2,nt);

					fftps.fourierWorkspace[i*n1+j]=(partial_t_rho+mxvalue+myvalue);	
				}
			}
		}

		fftps.get_fourier_coefficients(&fftps.u[n*n1*n2]);
		
	}

    fftps.forward_tridiagonal_sweep();    
    fftps.backward_tridiagonal_sweep();
    
    for(int n=0;n<nt;n++){
        for(int i=0;i<n1*n2;i++){
            fftps.fourierWorkspace[i]=fftps.u[n*n1*n2+i];
        }
        
        fftps.back_to_real_space();
        
        for(int i=0;i<n1*n2;i++){
            Phi[n*n1*n2+i]+=sigma_PDHG*fftps.fourierWorkspace[i];
        }
    } 

}

void update_mx(poisson_solver& fftps, double* mx,const double* mxTmp,const double* rho,const double* my,const double* Phi,double tau_PDHG,int n1,int n2,int nt,double dx,double dy,double dt){

	for(int n=0;n<nt;++n){
		for(int i=0;i<n2;++i){
			for(int j=0;j<n1-1;++j){
				double mxvalue,myvalue;

				mxvalue=mx[n*n1*n2+i*n1+j];
				if(i==0){
					myvalue=0.25*(my[n*n1*n2+i*n1+j]+my[n*n1*n2+i*n1+j+1]);
				}else{
					myvalue=0.25*(my[n*n1*n2+i*n1+j]+my[n*n1*n2+(i-1)*n1+j]+my[n*n1*n2+(i-1)*n1+j+1]+my[n*n1*n2+i*n1+j+1]);	
				}

				double rhovalue=0.5*(rho[n*n1*n2+i*n1+j]+rho[n*n1*n2+i*n1+j+1]);
				double grad_x_Phi=1.0/dx*(Phi[n*n1*n2+i*n1+j+1]-Phi[n*n1*n2+i*n1+j]);
				double grad_y_Phi=0.25/dy*(Phi[n*n1*n2+int(fmin(i+1,n2-1))*n1+j]-Phi[n*n1*n2+int(fmax(i-1,0))*n1+j])+0.25/dy*(Phi[n*n1*n2+int(fmin(i+1,n2-1))*n1+j+1]-Phi[n*n1*n2+int(fmax(i-1,0))*n1+j+1]);

				double newmxvalue;

				double c1=mxvalue+tau_PDHG*grad_x_Phi;
				double c2=myvalue+tau_PDHG*grad_y_Phi;
				double mu=tau_PDHG/rhovalue;

				if((c2/(1.0+mu)>=c1 && c1>=0)||(c2/(1.0+mu)<=c1 && c1<=0)){
					newmxvalue=c1;
				}else if((c1/(1.0+mu)>=c2 && c2>=0)||(c1/(1.0+mu)<=c2 && c2<=0)){
					newmxvalue=c1/(1+mu);
				}else if((c2>-(1+mu)/mu*c1 && c2<-mu/(1+mu)*c1)||(c2>-mu/(1+mu)*c1 && c2<-(1+mu)/mu*c1)){
					newmxvalue=((1+mu)*c1+mu*c2)/(1+2*mu);
				}else if((c2>=1.0/(1+mu)*c1 && c2<=(1+mu)*c1)||(c2<=1.0/(1+mu)*c1 && c2>=(1+mu)*c1)){
					newmxvalue=(c1+c2)/(2+mu);
				}else if((c2<=0 && c2>=-mu/(1+mu)*c1)||(c2>=0 && c2<=-mu/(1+mu)*c1)){
					newmxvalue=c1/(1+mu);
				}else if((c1<=0 && c2>=-(1+mu)/mu*c1)||(c1>=0 && c2<=-(1+mu)/mu*c1)){
					newmxvalue=0;
				}else{
					cout<<"wrong answer"<<endl;
					exit(1);
				}
				mx[n*n1*n2+i*n1+j]=newmxvalue;
			}
		}
	}
	
}

void update_my(poisson_solver& fftps, double* my,const double* myTmp,const double* rho,const double* mx,const double* Phi,double tau_PDHG,int n1,int n2,int nt,double dx,double dy,double dt){

	for(int n=0;n<nt;++n){
		for(int i=0;i<n2-1;++i){
			for(int j=0;j<n1;++j){
				double mxvalue,myvalue;

				myvalue=my[n*n1*n2+i*n1+j];
				if(j==0){
					mxvalue=0.25*(mx[n*n1*n2+i*n1+j]+mx[n*n1*n2+(i+1)*n1+j]);	
				}else{
					mxvalue=0.25*(mx[n*n1*n2+i*n1+j]+mx[n*n1*n2+(i+1)*n1+j]+mx[n*n1*n2+(i+1)*n1+j-1]+mx[n*n1*n2+i*n1+j-1]);	
				}
				

				double rhovalue=0.5*(rho[n*n1*n2+i*n1+j]+rho[n*n1*n2+(i+1)*n1+j]);
				double grad_y_Phi=1.0/dy*(Phi[n*n1*n2+(i+1)*n1+j]-Phi[n*n1*n2+i*n1+j]);
				double grad_x_Phi=0.25/dx*(Phi[n*n1*n2+i*n1+int(fmin(j+1,n1-1))]-Phi[n*n1*n2+i*n1+int(fmax(j-1,0))])+0.25/dx*(Phi[n*n1*n2+(i+1)*n1+int(fmin(j+1,n1-1))]-Phi[n*n1*n2+(i+1)*n1+int(fmax(j-1,0))]);

				double newmyvalue;

				double c1=mxvalue+tau_PDHG*grad_x_Phi;
				double c2=myvalue+tau_PDHG*grad_y_Phi;
				double mu=tau_PDHG/rhovalue;

				if((c2/(1.0+mu)>=c1 && c1>=0)||(c2/(1.0+mu)<=c1 && c1<=0)){
					newmyvalue=c2/(1+mu);
				}else if((c1/(1.0+mu)>=c2 && c2>=0)||(c1/(1.0+mu)<=c2 && c2<=0)){
					newmyvalue=c2;
				}else if((c2>-(1+mu)/mu*c1 && c2<-mu/(1+mu)*c1)||(c2>-mu/(1+mu)*c1 && c2<-(1+mu)/mu*c1)){
					newmyvalue=((1+mu)*c2+mu*c1)/(1+2*mu);
				}else if((c2>=1.0/(1+mu)*c1 && c2<=(1+mu)*c1)||(c2<=1.0/(1+mu)*c1 && c2>=(1+mu)*c1)){
					newmyvalue=(c1+c2)/(2+mu);
				}else if((c2<=0 && c2>=-mu/(1+mu)*c1)||(c2>=0 && c2<=-mu/(1+mu)*c1)){
					newmyvalue=0;
				}else if((c1<=0 && c2>=-(1+mu)/mu*c1)||(c1>=0 && c2<=-(1+mu)/mu*c1)){
					newmyvalue=c2/(1+mu);
				}else{
					cout<<"wrong answer"<<endl;
					exit(1);
				}

				my[n*n1*n2+i*n1+j]=newmyvalue;
			}
		}
	}
	
}

double sign(double x){
    
    double s= (x>0) - (x<0);
    
    return s;
    
}

double real3rdRoot1=-.5;   // equals cos(2*M_PI/3);
double im3rdRoot1=0.86602540378;   //equals sin(2*M_PI/3);
double real3rdRoot2=-.5;  //  equals cos(4*M_PI/3)=real3rdRoot1;
double im3rdRoot2=-0.86602540378;  //equals sin(4*M_PI/3)=-im3rdRoot1;


double cubic_solve(double b, double c, double d){
    
    double b3over3=(b/3)*(b/3)*(b/3);
    
    double p=c-b*(b/3);
    double q=d+2*b3over3-b*(c/3);
    double solution=0;
    
    if(p==0){
        
        solution=-sign(q)*exp(log(fabs(q))/3.0);
        
    }else{
        double discrim=(q/2)*(q/2)+(p/3)*(p/3)*(p/3);
        
        double s=sqrt(fabs(discrim));
        
        if(discrim<0){
            
            double theta=atan2(s,-q/2);
            
            double x=s*s+q*q/4;
            double rc=exp(log(x)/6);
            
            double thetac=theta/3;
            
            double real=rc*cos(thetac);
            double im=rc*sin(thetac);
            
            double solution1=2*real;
            
            
            double solution2=2*(real*real3rdRoot1-im*im3rdRoot1);
            double solution3=2*(real*real3rdRoot2-im*im3rdRoot2);
            
            solution=fmax(solution1,fmax(solution2,solution3));
            
            
        }else if(discrim>0){
            
            double u3=-q/2+s;
            double v3=-q/2-s;
            
            double u=sign(u3)*exp(log(fabs(u3))/3);
            double v=sign(v3)*exp(log(fabs(v3))/3);
            
            solution=u+v;
            
        }else{
            solution=fmax(3*q/p, -3*q/(2*p));
            
        }
    }
    
    return solution-b/3;
    
}


void update_rho(poisson_solver& fftps, double* rho,const double* rhoTmp,const double* mx,const double* my,const double* Phi,double tau_PDHG,int n1,int n2,int nt,double dx,double dy,double dt){

	for(int n=1;n<nt-1;++n){
		for(int i=0;i<n2;++i){
			for(int j=0;j<n1;++j){
				double mxvalue;
				double myvalue;
				if(j==0){
					mxvalue=0.5*mx[n*n1*n2+i*n1+j];
				}else{
					mxvalue=0.5*(mx[n*n1*n2+i*n1+j]+mx[n*n1*n2+i*n1+j-1]);
				}
				if(i==0){
					myvalue=0.5*my[n*n1*n2+i*n1+j];
				}else{
					myvalue=0.5*(my[n*n1*n2+i*n1+j]+my[n*n1*n2+(i-1)*n1+j]);
				}
				
				double tropical_norm=calculate_tropical_norm(mxvalue,myvalue);
				double rhovalue=rho[n*n1*n2+i*n1+j];
				double rhoTmpvalue=rhoTmp[n*n1*n2+i*n1+j];
				double partial_t_Phi=0.5/dt*(Phi[(n+1)*n1*n2+i*n1+j]-Phi[(n-1)*n1*n2+i*n1+j]);

				double newrhovalue=cubic_solve(-rhoTmpvalue-tau_PDHG*partial_t_Phi, 0, -0.5*tau_PDHG*tropical_norm*tropical_norm);
				rho[n*n1*n2+i*n1+j]=fmax(0,newrhovalue);
			}
		}
	}
}

int main(int argc, char **argv)
{

    if(argc!=8){
        cout << "Need to do the following : " << endl;
        cout << "./main.exe [n1] [n2] [nt] [tau] [sigma] [tolerance] [max iteration]" << endl;
        exit(1);
    }

    // Parameters for Grids
    
    int n1 = stoi(argv[1]);  // x panels
    int n2 = stoi(argv[2]);  // y panels
    int nt = stoi(argv[3]);  // t panels
    double tau_PDHG=stod(argv[4]);
    double sigma_PDHG=stod(argv[5]);
    double tolerance = stod(argv[6]);
    int max_iteration_PDHG=stoi(argv[7]);
    
    create_csv_file_for_parameters(n1,n2,nt);

    double base=0.5;

    double dx=1.0/n1;
    double dy=1.0/n2;
    double dt=1.0/nt;

    cout<<"XXX G-Prox PDHG with Tropical Norm XXX"<<endl;
    cout<<endl;
    cout<<"n1 : "<<n1<<" n2 : "<<n2<<" nt : "<<nt<<" base : "<<base<<endl;
    cout<<"dx : "<<scientific<<dx<<endl;
    cout<<"dy : "<<scientific<<dy<<endl;
    cout<<"dt : "<<scientific<<dt<<endl;
    cout<<fixed;
    cout<<"tau_PDHG   : "<<tau_PDHG<<endl;
    cout<<"sigma_PDHG : "<<sigma_PDHG<<endl;
    cout<<"max_iteration_PDHG : "<<max_iteration_PDHG<<endl;
    cout<<"tolerance          : "<<scientific<<tolerance<<endl;


    double* rho=new double[n1*n2*nt];
    double* rhoTmp=new double[n1*n2*nt];
    double* mx=new double[n1*n2*nt];
    double* my=new double[n1*n2*nt];
    double* mxTmp=new double[n1*n2*nt];
    double* myTmp=new double[n1*n2*nt];
    double* Phi=new double[n1*n2*nt];

    clock_t t;
    t = clock();
    poisson_solver fftps(n1,n2,nt,dx,dy,dt);
    t = clock() - t;
    printf ("\nCPU time for setting up FFT: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);
  

	create_initial_densities_one_to_four(rho,base,n1,n2,nt,dx,dy,dt);	
	// create_initial_densities_one_to_one_negative(rho,base,n1,n2,nt,dx,dy,dt);
	// create_initial_densities_one_to_one_positive(rho,base,n1,n2,nt,dx,dy,dt);	

    double energy;
    double previous_energy=1;
    double error;
    int iterPDHG;

    cout<<"\nXXX Starting Iterations XXX"<<endl;

    t = clock();

    for(iterPDHG=0;iterPDHG<max_iteration_PDHG;++iterPDHG){

    	memcpy(rhoTmp,rho,n1*n2*nt*sizeof(double));

		update_rho(fftps,rho,rhoTmp,mx,my,Phi,tau_PDHG,n1,n2,nt,dx,dy,dt);

    	memcpy(mxTmp,mx,n1*n2*nt*sizeof(double));
    	memcpy(myTmp,my,n1*n2*nt*sizeof(double));

		update_mx(fftps,mx,mxTmp,rho,my,Phi,tau_PDHG,n1,n2,nt,dx,dy,dt);
		update_my(fftps,my,myTmp,rho,mx,Phi,tau_PDHG,n1,n2,nt,dx,dy,dt);

    	for(int i=0;i<n1*n2*nt;++i){
    		rhoTmp[i]=2*rho[i]-rhoTmp[i];
    		mxTmp[i]=2*mx[i]-mxTmp[i];
    		myTmp[i]=2*my[i]-myTmp[i];
    	}

    	update_Phi(fftps,Phi,rhoTmp,mxTmp,myTmp,sigma_PDHG,n1,n2,nt,dx,dy,dt);

    	energy = calculateEnergy(rho, mx, my,n1, n2, nt, dx, dy, dt);
    	error = fabs((energy-previous_energy)/previous_energy);
		previous_energy = energy;

    	double dual_value = calculate_dual(Phi,rho,n1,n2,nt,dx,dy,dt);
    	double dual_gap = fabs(energy-dual_value);

    	if((iterPDHG+1)%100==0){
    		cout<<scientific;
    		cout<<"iter : "<<setw(5)<<iterPDHG+1<<" tau : "<<setw(5)<<tau_PDHG<<" sigma : "<<setw(5)<<sigma_PDHG<<"  energy : "<<scientific<<setw(12)<<energy<<"  dual gap : "<<scientific<<setw(12)<<dual_gap<<" relative error : "<<scientific<<setw(12)<<error<<endl;
    	}

    	if(error<tolerance){
    		break; 
    	}
    }

    t = clock() - t;

    printf ("CPU time for Iterations: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);

    cout<<"Iterations     : "<<iterPDHG<<endl;
    cout<<"Energy         : "<<energy<<endl;
    cout<<"Relative Error : "<<error<<endl;

    create_csv_file(rho,"./data/rho.csv",n1,n2,nt);

	free(rho);
    free(rhoTmp);
    free(mx);
    free(my);
    free(mxTmp);
    free(myTmp);
    free(Phi);
    fftps.destroy_all_fftps();
}