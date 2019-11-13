#include <iostream>
#include <fstream>
#include <cmath>
#include <fftw3.h>
#include <time.h>

using namespace std;

class poisson_solver{
public:
    fftw_plan planIn;
    fftw_plan planOut;
    double *fourierWorkspace;
    double *u;
    double *kernel;

    int n1;
    int n2;
    double dx;
    double dy;

    poisson_solver(int n1, int n2, double dx, double dy) {
    	this->n1=n1;
    	this->n2=n2;
    	this->dx=dx;
    	this->dy=dy;

        fourierWorkspace =(double*) fftw_malloc(n1*n2*sizeof(double));

		planIn = fftw_plan_r2r_2d(n2,n1, fourierWorkspace, fourierWorkspace, FFTW_REDFT10,FFTW_REDFT10, FFTW_MEASURE);
		planOut = fftw_plan_r2r_2d(n2,n1, fourierWorkspace, fourierWorkspace, FFTW_REDFT01,FFTW_REDFT01, FFTW_MEASURE);

        kernel=new double[n1*n2];
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

    void perform_inverse_laplacian(){

        fftw_execute(planIn);
        fourierWorkspace[0]=0;

        for(int i=1;i<n1*n2;++i){
            fourierWorkspace[i]/=4*(n1)*(n2)*kernel[i];
        }

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

void create_csv_file(double* rho, int size, string filename){
	ofstream outfile;
	outfile.open(filename);
	for(int i=0;i<size;++i){
		outfile<<rho[i]<<" ";
	}
	outfile.close();
}


void create_initial_densities_two_boxes_positive(double* rho0, double* rho1, int n1, int n2, double dx, double dy){
    double sum0=0;
    double sum1=0;

    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            double x=j*dx;
            double y=i*dy;

            if(abs(x-1.0/3)<0.1 && abs(y-1.0/3)<0.1){
             rho0[i*n1+j]=1; 
            }

            if(abs(x-2.0/3)<0.1 && abs(y-2.0/3)<0.1){
             rho1[i*n1+j]=1; 
            }

            sum0+=rho0[i*n1+j];
            sum1+=rho1[i*n1+j];
        }
    }

    for(int i=0;i<n1*n2;++i){
        rho0[i]/=sum0*dx*dy;
        rho1[i]/=sum1*dy*dy;
    }

}

void create_initial_densities_two_boxes_negative(double* rho0, double* rho1, int n1, int n2, double dx, double dy){
    double sum0=0;
    double sum1=0;

    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            double x=j*dx;
            double y=i*dy;

            if(abs(x-1.0/3)<0.1 && abs(y-2.0/3)<0.1){
             rho0[i*n1+j]=1; 
            }

            if(abs(x-2.0/3)<0.1 && abs(y-1.0/3)<0.1){
             rho1[i*n1+j]=1; 
            }

            sum0+=rho0[i*n1+j];
            sum1+=rho1[i*n1+j];
        }
    }

    for(int i=0;i<n1*n2;++i){
        rho0[i]/=sum0*dx*dy;
        rho1[i]/=sum1*dy*dy;
    }
}

void create_initial_densities_c4(double* rho0, double* rho1, int n1, int n2, double dx, double dy){
    double sum0=0;
    double sum1=0;

    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            double x=j*dx;
            double y=i*dy;

            if(abs(x-0.5)<0.1 && abs(y-0.5)<0.1){
             rho0[i*n1+j]=1; 
            }
            if(abs(x-0.2)<0.05 && abs(y-0.2)<0.05){
             rho1[i*n1+j]=+1; 
            }
            if(abs(x-0.2)<0.05 && abs(y-0.8)<0.05){
             rho1[i*n1+j]=+1; 
            }
            if(abs(x-0.8)<0.05 && abs(y-0.2)<0.05){
             rho1[i*n1+j]=+1; 
            }
            if(abs(x-0.8)<0.05 && abs(y-0.8)<0.05){
             rho1[i*n1+j]=+1; 
            }

            sum0+=rho0[i*n1+j];
            sum1+=rho1[i*n1+j];
        }
    }

    for(int i=0;i<n1*n2;++i){
        rho0[i]/=sum0*dx*dy;
        rho1[i]/=sum1*dy*dy;
    }
}

void shrink(double* xx, double mx, double my, double h) {
    double* y=new double[2];
    y[0]=0;y[1]=0;
    if(mx < my) {
        shrink(y,my,mx,h);
        xx[0] = y[1];
        xx[1] = y[0];
    } else {
        xx[0]=0;
        xx[1]=0;
        if(mx < 0) {
            shrink(y,-my,-mx,h);
            xx[0] = -y[1];
            xx[1] = -y[0];
        } else {
            if(my>=0){
                if(mx>=my+1){
                    xx[0]=h*(mx-1);
                    xx[1]=h*my;
                }else if(mx+my>=1){
                    xx[0]=h*(mx+my-1)/2;
                    xx[1]=h*(mx+my-1)/2;
                }
            }else{
                if(mx>=1){
                    xx[0]=h*(mx-1);
                }
                if(my<=-1){
                    xx[1]=h*(my+1);
                }
            }
        }
    }
    free(y);
}

double calculate_tropical_norm(double mx, double my){
    if(mx>=my && my>=0){
        return mx;
    }else if(my>=mx && mx>=0){
        return my;
    }else if(mx>=0 && 0 >= my){
        return mx - my;
    }else if(my>=0 && 0 >= mx){
        return my - mx;
    }else if(0>=mx && mx >= my){
        return - my;
    }else if(0>=my && my >= mx){
        return - mx;
    }

    cout <<"wrong answer"<<endl;;

    exit(1);
    return 0;
}

double calculate_primal(double* mx, double* my, int n1, int n2, double dx, double dy){
	double sum=0;
	for(int i=0;i<n2;++i){
		for(int j=0;j<n1;++j){
			double mxvalue=mx[i*n1+j];
			double myvalue=my[i*n1+j];
			sum+=calculate_tropical_norm(mxvalue,myvalue);
		}
	}
	return sum*dx*dy;
}

void update_m(double* mx, double* my, const double* mxTmp, const double* myTmp, const double* Phi, double* xx, double h, int n1, int n2, double dx, double dy){

    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            double mxvalue=mxTmp[i*n1+j]+h/dx*(Phi[i*n1+int(fmin(n1-1,j+1))]-Phi[i*n1+j]);
            double myvalue=myTmp[i*n1+j]+h/dy*(Phi[int(fmin(n2-1,i+1))*n1+j]-Phi[i*n1+j]);

            shrink(xx,mxvalue,myvalue,h);

            mx[i*n1+j]=xx[0];
            my[i*n1+j]=xx[1];
        }
    }
}

void update_Phi(poisson_solver& fftps,double* Phi,const double* mx,const double* my,const double* mxTmp,const double* myTmp,const double* rho0, const double* rho1, double tau,int n1,int n2,double dx,double dy){

    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            double mxvalue,myvalue;
            if(j>0){
                mxvalue=2.0*mx[i*n1+j]-mxTmp[i*n1+j]-2.0*mx[i*n1+j-1]+mxTmp[i*n1+j-1];
            }else {
                mxvalue=2.0*mx[i*n1+j]-mxTmp[i*n1+j];
            }

            if(i>0){
                myvalue=2.0*my[i*n1+j]-myTmp[i*n1+j]-2.0*my[(i-1)*n1+j]+myTmp[(i-1)*n1+j];
            }else {
                myvalue=2.0*my[i*n1+j]-myTmp[i*n1+j];
            }

            double div_value=mxvalue/dx+myvalue/dy;

            fftps.fourierWorkspace[i*n1+j]=div_value+rho1[i*n1+j]-rho0[i*n1+j];
        }
    }

    fftps.perform_inverse_laplacian();

    for(int i=0;i<n1*n2;++i){
        Phi[i]+=tau*fftps.fourierWorkspace[i];
    }
}

int main(int argc, char **argv)
{

    if(argc!=7){
        cout << "Need to do the following : " << endl;
        cout << "./main.exe [n1] [n2] [h] [tau] [tolerance] [max_iteration]" << endl;
        exit(1);
    }

    // Parameters for Grids
    
    int n1 = stoi(argv[1]);  // x panels
    int n2 = stoi(argv[2]);  // y panels

    // Parameters for iterations

    double h   = stod(argv[3]);
    double tau = stod(argv[4]);

    double tolerance  = stod(argv[5]);
    int max_iteration = stoi(argv[6]);

    // delta x and delta y

    double dx=1.0/n1;
    double dy=1.0/n2;

    ofstream outfile;
    outfile.open("./data/Parameters.csv");
    outfile<<n1<<" "<<n2;
    outfile.close();

    double* rho0=new double[n1*n2];
    double* rho1=new double[n1*n2];
    double* Phi=new double[n1*n2];

    double* mx=new double[n1*n2];
    double* my=new double[n1*n2];

    double* mxTmp=new double[n1*n2];
    double* myTmp=new double[n1*n2];

    double* xx=new double[2];

    // --------------- Initial Densities ------------------

    // create_initial_densities_two_boxes_positive(rho0,rho1,n1,n2,dx,dy);
    // create_initial_densities_two_boxes_negative(rho0,rho1,n1,n2,dx,dy);
    create_initial_densities_c4(rho0,rho1,n1,n2,dx,dy);

    // ----------------------------------------------------

    create_csv_file(rho0,n1*n2, "./data/rho0.csv");
    create_csv_file(rho1,n1*n2, "./data/rho1.csv");

    cout<<"n1: "<<n1<<" n2: "<<n2<<" h: "<<h<<" tau: "<<tau<<endl;
    cout<<"max iteration : "<<max_iteration<<" tolerance : "<<tolerance<<endl;

    double energy=1;
    double previous_energy=1;
    double error=10;
    int iter;

    clock_t t;
    t=clock();
    poisson_solver fftps(n1,n2,dx,dy);
    t=clock()-t;

    printf ("\nCPU time for setting up FFT : %f seconds\n",((float)t)/CLOCKS_PER_SEC);

    cout<<"\nXXX Starting G-Prox PDHG XXX\n"<<endl;
    t=clock();

    for(iter=0;iter<max_iteration;++iter){
    	memcpy(mxTmp, mx, n1*n2*sizeof(double));
    	memcpy(myTmp, my, n1*n2*sizeof(double));

        update_m(mx,my,mxTmp,myTmp,Phi,xx,h,n1,n2,dx,dy);
        update_Phi(fftps,Phi,mx,my,mxTmp,myTmp,rho0,rho1,tau,n1,n2,dx,dy);
    
    	energy=calculate_primal(mx,my,n1,n2,dx,dy);

    	error=fabs((energy-previous_energy)/previous_energy);

    	previous_energy=energy; 

        if(iter % 100 == 99){
            cout<<"iter : "<<iter<<" error : "<<error<<endl;
        }

        if(error < tolerance){
            break;
        }	
    }

    t=clock()-t;

    cout<<"Total Iterations :" << iter<<endl;
    cout<<"Relative Error   :" << error<<endl;
    printf ("CPU time for iterations : %f seconds\n",((float)t)/CLOCKS_PER_SEC);

    create_csv_file(mx,n1*n2, "./data/mx.csv");
    create_csv_file(my,n1*n2, "./data/my.csv");

    free(rho0);
    free(rho1);
    free(Phi);
    free(mx);
    free(my);

    fftps.destroy_all_fftps();
}
