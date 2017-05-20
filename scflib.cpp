#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <time.h>
#include <iomanip>
extern"C"
{
#include "cint.h"
#include "mkl.h"
}

#define _USE_MATH_DEFINES

using namespace std;

int Cint_P_Forming()
{
        cout<<"__________Density Matrix Formation Subroutine__________"<<endl;
	ifstream f_para;
	ofstream f_P;
	f_P.open("data/P_density.txt");
	f_para.open("data/para_c.txt");
	int N_con,N_e,i,m,n;
	string buf;
	cout<<"Reading combination parameters."<<endl;
	f_para>>N_e>>N_con;
	cout<<"Number of electrons:"<<N_e<<endl;
	cout<<"Number of contracted gaussians:"<<N_con<<endl;
        double *P=new double[N_con*N_con];
	double *C=new double[N_con*N_con];
        while(!f_para.eof()&&f_para.peek()!=EOF)
	{
		/*In the cycles we use in.peek() to read the next char to the current file pointer,thus we prevent taking in EOF note 0xff as the last line!(WHen we finished reading the last line, we actually hadn't reach the EOF so we will go into the next cycle, and take the single EOF note as the last line!)*/
		getline(f_para,buf);
	}
	istringstream ss_para;
	ss_para.str(buf);
	for (i=0;i<N_con*N_con;i++)
	{
		ss_para>>C[i];
		cout<<setprecision(12)<<C[i]<<' ';
	}
	cout<<endl;
	cout<<"Parameters read.";
	for (m=0;m<N_con;m++)
	{
		for(n=0;n<N_con;n++)
		{
		       P[m*N_con+n]=0;
               	       for (i=0;i<N_e/2;i++)
			       P[m*N_con+n]+=2*C[m*N_con+i]*C[n*N_con+i];
	         	/*Only valid for close shell structures.*/
		}
	}
	/*The initial guess of MO parameters shall be sorted by eigenvalues in an increasing order.*/
	cout<<"Density matrix elements:"<<endl;
	for (i=0;i<N_con*N_con;i++)
	{
			cout<<setprecision(12)<<P[i]<<' ';
			f_P<<setprecision(12)<<P[i]<<' ';
	}
	cout<<endl;
	f_P<<endl;
	cout<<"Density matrix written."<<endl;
	f_para.close();
	f_P.close();
	delete []P;
	delete []C;
	return 0;
}
/*This subprogram only works for restricted closed shell.*/


int Cint_F_Forming()
{
	cout<<"__________Fock Matrix Forming Subroutine__________"<<endl;
	ifstream f_P;
	ifstream f_T;
	ifstream f_V;
	ifstream f_2e;
	ifstream f_para;
	ofstream f_F;
	f_P.open("data/P_density.txt");
	f_T.open("data/T_integral.txt");
	f_V.open("data/V_integral.txt");
	f_para.open("data/para_c.txt");
	f_2e.open("data/2e_integral.txt");
	f_F.open("data/Fock.txt");
	int N_e,N_con,i,k,l,m,n;
	f_para>>N_e>>N_con;
	cout<<"Number of electrons:"<<N_e<<endl;
	cout<<"Number of contracted gaussians:"<<N_con<<endl;
	double *P=new double[N_con*N_con];
	double *T=new double[N_con*N_con];
	double *V=new double[N_con*N_con];
	double *e2=new double[N_con*N_con*N_con*N_con];
	double *F=new double[N_con*N_con];
	cout<<"Reading density matrix"<<endl;
	for (i=0;i<N_con*N_con;i++)
	{
		f_P>>P[i];
		cout<<setprecision(12)<<P[i]<<' ';
	}
	cout<<endl;
	cout<<"Fock matrix elements:"<<endl;
        for (k=0;k<N_con;k++)
	{
		for (l=0;l<N_con;l++)
		{
			f_T>>T[k*N_con+l];
			f_V>>V[k*N_con+l];
			F[k*N_con+l]=T[k*N_con+l]+V[k*N_con+l];
			for (m=0;m<N_con;m++)
				for (n=0;n<N_con;n++)
					f_2e>>e2[k*N_con*N_con*N_con+l*N_con*N_con+m*N_con+n];
		}
	}
       for (k=0;k<N_con;k++)
	{	
		for (l=0;l<N_con;l++)
		{
	                for (m=0;m<N_con;m++)
			        for (n=0;n<N_con;n++)
                                        F[k*N_con+l]+=P[m*N_con+n]*(e2[k*N_con*N_con*N_con+l*N_con*N_con+m*N_con+n]-0.5*e2[k*N_con*N_con*N_con+n*N_con*N_con+m*N_con+l]);
		}
	}
        for (i=0;i<N_con*N_con;i++)
        {
                cout<<setprecision(12)<<F[i]<<' ';
                f_F<<setprecision(12)<<F[i]<<' ';
        }
	cout<<endl;
	f_F<<endl;
	cout<<"Finshed reading integrals and Fock calculation!"<<endl;
	f_P.close();
	f_T.close();
	f_V.close();
	f_para.close();
	f_2e.close();
	f_F.close();
	delete []P;
	delete []T;
	delete []V;
	delete []e2;
	delete []F;
	return 0;
}
/*This subprogram forms the Fock matrix from the calculated and remapped 1e and 2e integrals matrix and the density matrix under close-shell restricted approximation.*/

int Cint_F_Rescaling()
{
        cout<<"_________Fock Rescaling Subroutine__________"<<endl;
	ifstream f_Fock;
	ifstream f_para;
	ifstream f_X;
	ofstream f_Fock_pri;
	f_Fock.open("data/Fock.txt");
	f_para.open("data/para_c.txt");
	f_X.open("data/X_transmat.txt");
	f_Fock_pri.open("data/Fock_pri.txt");
	int i;
	int N_e,N_con;
	f_para>>N_e>>N_con;
	cout<<"Number of electrons:"<<N_e<<endl;
	cout<<"Number of contracted gaussians:"<<N_con<<endl;
	double *F=new double[N_con*N_con];
        double *F_pri=new double[N_con*N_con];
        double *X=new double[N_con*N_con];
	double *C=new double[N_con*N_con];
	cout<<"Reading original Fock Matrix and the rescaling matrix X,and initializing buffers"<<endl;
	for (i=0;i<N_con*N_con;i++)
	{
		f_Fock>>F[i];
		F_pri[i]=0;
		f_X>>X[i];
		C[i]=0;
	}
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,N_con,N_con,N_con,1.0,F,N_con,X,N_con,0,C,N_con);
        cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,N_con,N_con,N_con,1.0,X,N_con,C,N_con,0,F_pri,N_con);
        cout<<"Fock rescalation successful!Fock prime:"<<endl;
	for (i=0;i<N_con*N_con;i++)
	{
		cout<<setprecision(12)<<F_pri[i]<<' ';
		f_Fock_pri<<setprecision(12)<<F_pri[i]<<' ';
	}
	cout<<endl;
	f_Fock_pri<<endl;
	f_Fock.close();
	f_para.close();
	f_X.close();
	f_Fock_pri.close();
	delete []F;
        delete []F_pri;
	delete []X;
	delete []C;
	return 0;
}

bool Cint_F_Diag_And_C_New(double eigenval_conv)
{
	cout<<"_________Fock_prime Diagonization and C_prime Rescaling Subroutine________"<<endl;
	ifstream f_Fock_pri;
	ifstream f_para;
	ifstream f_X;
	ifstream f_eigenval_former;
	ofstream f_new_C;
	ofstream f_eigenval;
	f_Fock_pri.open("data/Fock_pri.txt");
	f_para.open("data/para_c.txt");
	f_X.open("data/X_transmat.txt");
	f_eigenval_former.open("data/eigenval.txt");
	int N_e,N_con,i,j;
	bool EIG_CONVERGE=true;
	f_para>>N_e>>N_con;
	cout<<"Number of electrons:"<<N_e<<endl;
	cout<<"Number of contracted gaussians"<<N_con<<endl;
	f_para.close();
	f_new_C.open("data/para_c.txt",ios::app);//Use appenditive output.
	double *F_pri=new double[N_con*N_con];
	double *C_new=new double[N_con*N_con];
	double *eigs_former=new double[N_con];
        double *eigs=new double[N_con];
	double *e=new double[N_con-1];
	double *tau=new double[N_con-1];
	double *X=new double[N_con*N_con];
	int *ipiv=new int[N_con];
	for (i=0;i<N_con*N_con;i++)
	{
		f_Fock_pri>>F_pri[i];
		f_X>>X[i];
		C_new[i]=0;
	}
	cout<<"F_pri and X matrices read."<<endl;
	cout<<"Reading former eigenvalues"<<endl;
	string buf;
        while(!f_eigenval_former.eof()&&f_eigenval_former.peek()!=EOF)
	{
		/*In the cycles we use in.peek() to read the next char to the current file pointer,thus we prevent taking in EOF note 0xff as the last line!(WHen we finished reading the last line, we actually hadn't reach the EOF so we will go into the next cycle, and take the single EOF note as the last line!)*/
		getline(f_eigenval_former,buf);
	}
	istringstream ss_eigenval_former;
	ss_eigenval_former.str(buf);
	for (i=0;i<N_con;i++)
		eigs_former[i]=0;
	for (i=0;i<N_con;i++)
	{
		ss_eigenval_former>>eigs_former[i];
		cout<<setprecision(12)<<eigs_former[i]<<' ';
	}
	cout<<endl;
	cout<<"Former eigenvalues read."<<endl;
	f_eigenval_former.close();
	if (LAPACKE_dsytrd(LAPACK_ROW_MAJOR,'U',N_con,F_pri,N_con,eigs,e,tau)==0)
		cout<<"Tridiagonal trasformation successful!"<<endl;
	if (LAPACKE_dorgtr(LAPACK_ROW_MAJOR,'U',N_con,F_pri,N_con,tau)==0)
		cout<<"Tridiagonal transformation matrix Q written to matrix F_pri."<<endl;
	if (LAPACKE_dsteqr(LAPACK_ROW_MAJOR,'V',N_con,eigs,e,F_pri,N_con)==0)
		cout<<"Diagonal trasformation of F_pri finished!Final C_pri written to F_pri by column!Eigenvalues written to eigs!C_pri:"<<endl;
	for (i=0;i<N_con*N_con;i++)
		cout<<setprecision(12)<<F_pri[i]<<' ';
	cout<<endl;
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,N_con,N_con,N_con,1.0,X,N_con,F_pri,N_con,0,C_new,N_con);
	cout<<"New C generated!"<<endl;
	for (i=0;i<N_con*N_con;i++)
	{
                cout<<setprecision(12)<<C_new[i]<<' ';
	        f_new_C<<setprecision(12)<<C_new[i]<<' ';
	}
	cout<<endl;
	f_new_C<<endl;
	cout<<"New C written!"<<endl;
	f_eigenval.open("data/eigenval.txt",ios::app);
	for (i=0;i<N_con;i++)
	{
		cout<<setprecision(12)<<eigs[i]<<' ';
		f_eigenval<<setprecision(12)<<eigs[i]<<' ';
	}
	cout<<endl;
	f_eigenval<<endl;
	cout<<"Eigenvalues written!"<<endl;
	f_Fock_pri.close();
	f_X.close();
        f_new_C.close();
	f_eigenval.close();
	for (i=0;i<N_con;i++)
		EIG_CONVERGE=(EIG_CONVERGE and (abs(eigs[i]-eigs_former[i])<eigenval_conv));
	cout<<"Eigenvals converged?"<<EIG_CONVERGE<<endl;
	delete []F_pri;
	delete []C_new;
	delete []eigs;
	delete []eigs_former;
	delete []e;
	delete []tau;
	delete []X;
	delete []ipiv;
	return EIG_CONVERGE;
}

bool Cint_Elec_Energy(double ene_conv)
{
	cout<<"__________Total Energy Subroutine__________"<<endl;
	ifstream f_para;
	ifstream f_P;
	ifstream f_T;
	ifstream f_V;
	ifstream f_F;
	ifstream f_elecene_former;
	ofstream f_elecene;
	f_para.open("data/para_c.txt");
	f_P.open("data/P_density.txt");
	f_T.open("data/T_integral.txt");
	f_V.open("data/V_integral.txt");
	f_F.open("data/Fock.txt");
	f_elecene_former.open("data/elecene.txt");
	int N_e,N_con,i,j;
	bool ENE_CONVERGE=true;
	f_para>>N_e>>N_con;
	cout<<"Number of electrons:"<<N_e<<endl;
	cout<<"Number of contracted gaussians:"<<N_con<<endl;
	double *P=new double[N_con*N_con];
	double *T=new double[N_con*N_con];
	double *V=new double[N_con*N_con];
	double *F=new double[N_con*N_con];
	double elecene_former=0;
	string buf_E;
	while(!f_elecene_former.eof()&&f_elecene_former.peek()!=EOF)
		getline(f_elecene_former,buf_E);
	istringstream ss_elecene_former;
	ss_elecene_former.str(buf_E);
	cout<<"Reading last energy:"<<endl;
	ss_elecene_former>>elecene_former;
	cout<<elecene_former<<endl;
	cout<<"Last energy read!"<<endl;
	f_elecene_former.close();
	cout<<"Reading density matrix:"<<endl;
	for (i=0;i<N_con*N_con;i++)
	{
		f_P>>P[i];
		cout<<setprecision(12)<<P[i]<<' ';
	}
	cout<<endl;
	cout<<"Density matrix read."<<endl;
	cout<<"Reading kinetic integrals:"<<endl;
	for (i=0;i<N_con*N_con;i++)
	{
		f_T>>T[i];
		cout<<setprecision(12)<<T[i]<<' ';
	}
	cout<<endl;
	cout<<"Kinetic integrals read!"<<endl;
	cout<<"Reading Vnuc integrals:"<<endl;
	for (i=0;i<N_con*N_con;i++)
	{
		f_V>>V[i];
		cout<<setprecision(12)<<V[i]<<' ';
	}
	cout<<endl;
	cout<<"Vnuc integrals read!"<<endl;
	cout<<"Reading Fock matrix:"<<endl;
	for (i=0;i<N_con*N_con;i++)
	{
		f_F>>F[i];
		cout<<setprecision(12)<<F[i]<<' ';
	}
	cout<<endl;
	cout<<"Fock matrix read!"<<endl;
	f_elecene.open("data/elecene.txt",ios::app);
	double elecene_new=0;
	for (i=0;i<N_con;i++)
		for (j=0;j<N_con;j++)
			elecene_new+=0.5*P[i*N_con+j]*(T[j*N_con+i]+V[j*N_con+i]+F[j*N_con+i]);
	cout<<"New electronic energy:"<<setprecision(12)<<elecene_new<<endl;
	f_elecene<<setprecision(12)<<elecene_new<<endl;
	ENE_CONVERGE=ENE_CONVERGE and (abs(elecene_new-elecene_former)<ene_conv);
	f_para.close();
	f_P.close();
	f_T.close();
	f_V.close();
	f_F.close();
	f_elecene.close();
	delete []P;
	delete []T;
	delete []V;
	delete []F;
	cout<<"Electronic energy converged?"<<ENE_CONVERGE<<endl;
	return ENE_CONVERGE;
}

int Cint_Summary(int n)
{
	cout<<"_________OOParts Summary Subroutine_________"<<endl;
	ifstream f_elecene;
	ifstream f_eigenval;
	ifstream f_eigenvec;
	ifstream f_nucene;
	f_elecene.open("data/elecene.txt");
	f_eigenval.open("data/eigenval.txt");
	f_eigenvec.open("data/para_c.txt");
	f_nucene.open("data/nucene.txt");
	int N_e,N_con,i,j;
	f_eigenvec>>N_e>>N_con;
	cout<<"Number of electrons:"<<N_e<<endl;
	cout<<"Number of contracted gaussians(=N_mo):"<<N_con<<endl;
	string buf_E,buf_es,buf_ev;
	while(!f_elecene.eof()&&f_elecene.peek()!=EOF)
		getline(f_elecene,buf_E);
	while(!f_eigenval.eof()&&f_eigenval.peek()!=EOF)
		getline(f_eigenval,buf_es);
	while(!f_eigenvec.eof()&&f_eigenvec.peek()!=EOF)
		getline(f_eigenvec,buf_ev);
	istringstream ene,val,vec;
	ene.str(buf_E);
	val.str(buf_es);
	vec.str(buf_ev);
	double elecene;
	double *eigenval=new double[N_con];
	double *eigenvec=new double[N_con*N_con];
	ene>>elecene;
	cout<<"Convergence after "<<n<<" SCF cycles."<<endl;
	cout<<"Electronic Energy Ee(RHF):"<<setprecision(12)<<elecene<<endl;
	cout<<"MO energies(ascending order):"<<endl;
	for (i=0;i<N_con;i++)
	{
		val>>eigenval[i];
		cout<<setprecision(12)<<eigenval[i]<<' ';
	}
	cout<<endl;
	cout<<"LCAO-MO parameters(in correspondence with MO's above):"<<endl;
	for (i=0;i<N_con*N_con;i++)
		vec>>eigenvec[i];
	for (i=0;i<N_con;i++)
	{
		for (j=0;j<N_con;j++)
			cout<<setprecision(12)<<eigenvec[j*N_con+i]<<' ';
		cout<<endl;
	}
	double nucene;
	f_nucene>>nucene;
	cout<<"Nuclear repulsion energy:"<<setprecision(12)<<nucene<<endl;
	cout<<"Total energy:"<<setprecision(12)<<(nucene+elecene)<<endl;
	f_elecene.close();
	f_eigenval.close();
	f_eigenvec.close();
	f_nucene.close();
	delete []eigenval;
	delete []eigenvec;
	return 0;
}
