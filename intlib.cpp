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


int ind_remap1e(int i,int j,int k,int l,int N_con,int *bas)
{
	int sum1=0;
	int sum2=0;
	int x;
	for (x=0;x<i;x++)
		sum1+=bas[x*8+3];
	for (x=0;x<k;x++)
		sum2+=bas[x*8+3];
	return (sum1+j)*N_con+sum2+l;
}

int ind_remap2e(int i,int j,int k,int l,int p,int q,int r,int s,int N_con,int *bas)
{
	int sum1=0;
	int sum2=0;
	int sum3=0;
	int sum4=0;
	int x;
	for (x=0;x<i;x++)
		sum1+=bas[x*8+3];
	for (x=0;x<k;x++)
		sum2+=bas[x*8+3];
	for (x=0;x<p;x++)
		sum3+=bas[x*8+3];
	for (x=0;x<r;x++)
		sum4+=bas[x*8+3];
	return (sum1+j)*N_con*N_con*N_con+(sum2+l)*N_con*N_con+(sum3+q)*N_con+sum4+s;
}


int Cint_All1e()
{
	cout<<"_________1e Integration Subroutine__________"<<endl;
	ofstream f_S;
	ofstream f_V;
	ofstream f_T;
	ifstream f_env;
	ifstream f_atm;
	ifstream f_bas;
	ifstream f_para;
	f_S.open("data/S_integral.txt",ios::out);
        f_V.open("data/V_integral.txt",ios::out);
        f_T.open("data/T_integral.txt",ios::out);
	f_env.open("data/env.txt");
	f_atm.open("data/atom.txt");
	f_bas.open("data/bas.txt");
	f_para.open("data/para_c.txt");
	/*env array of the called cint function has been written according to cint manual in the info.txt file*/
        int i,j,k,l,natm,nbas,N_con,N_e;
	f_atm>>natm;
	cout<<"Number of atoms:"<<natm<<endl;
	f_bas>>nbas;
        cout<<"Number of shells:"<<nbas<<endl;
	f_para>>N_e>>N_con;
	cout<<"Number of electrons:"<<N_e<<endl;
	cout<<"Number of contracted gaussians:"<<N_con<<endl;
	int *atm=new int[natm*6];
	int *bas=new int[nbas*8];
	cout<<"Reading atomic information.atm:"<<endl;
        for (i=0;i<natm*6;i++)
	{
		f_atm>>atm[i];
		cout<<atm[i]<<' ';
	}
	cout<<endl;
	cout<<"Reading basis information,bas:"<<endl;
	for (i=0;i<nbas*8;i++)
	{
		f_bas>>bas[i];
		cout<<bas[i]<<' ';
	}
	cout<<endl;
	int size_env = natm*3;
	for (i=0;i<nbas;i++)
		size_env+=(bas[i*8+2]+bas[i*8+2]*bas[i*8+3]);
	double *env = new double[size_env];
	f_env>>natm>>nbas;
	for (i=0;i<size_env;i++)
		f_env>>env[i];
        int m,n;
	double *all_s_remap = new double[N_con*N_con];
	double *all_v_remap = new double[N_con*N_con];
	double *all_t_remap = new double[N_con*N_con];
	for (m=0;m<nbas;m++)
	{
		for (n=0;n<nbas;n++)
		{
	               int shls[2]={n,m};
	               double *buf_s =new double[bas[m*8+3]*bas[n*8+3]];
		       /*We need to consider different polarization directions in the same shell.Revise after consulting libcint author!*/
		       /*No revision required since bas[i*8+1] already has included polarization direction, ex. bas[i*8+1]=-1,0 or 1 stands for three different p func!*/
	               double *buf_v =new double[bas[m*8+3]*bas[n*8+3]];
                       double *buf_t =new double[bas[m*8+3]*bas[n*8+3]];
		       cint1e_ovlp_cart(buf_s,shls,atm,natm,bas,nbas,env);
                       cout<<"Overlap "<<n<<','<<m<<"successful!:"<<endl;
		       for (i=0;i<bas[n*8+3]*bas[m*8+3];i++)
		       {
		               f_S<<setprecision(12)<<buf_s[i]<<' ';
			       cout<<setprecision(12)<<buf_s[i]<<' ';
		       }
		       cout<<endl;
		       delete []buf_s;
                       cint1e_nuc_cart(buf_v,shls,atm,natm,bas,nbas,env);
                       cout<<"NuclearV "<<n<<','<<m<<"successful!:"<<endl;
		       for (i=0;i<bas[n*8+3]*bas[m*8+3];i++)
		       {
			       f_V<<setprecision(12)<<buf_v[i]<<' ';
			       cout<<setprecision(12)<<buf_v[i]<<' ';
		       }
		       cout<<endl;
		       delete []buf_v;
                       cint1e_kin_cart(buf_t,shls,atm,natm,bas,nbas,env);
		       cout<<"Kinetic "<<n<<','<<m<<"successful!:"<<endl;
		       for (i=0;i<bas[n*8+3]*bas[m*8+3];i++)
		       {
			       f_T<<setprecision(12)<<buf_t[i]<<' ';
			       cout<<setprecision(12)<<buf_t[i]<<' ';
		       }
		       cout<<endl;
		       delete []buf_t;
		}
	}
	cout<<"Setup files (atomic,basis,environmental) closed!"<<endl;
	f_S<<endl;
	f_S.close();
	f_V<<endl;
	f_V.close();
	f_T<<endl;
	f_T.close();
	f_atm.close();
	f_bas.close();
        f_env.close();
	f_para.close();
	ifstream f2_S;
	ifstream f2_V;
	ifstream f2_T;
	f2_S.open("data/S_integral.txt");
        f2_V.open("data/V_integral.txt");
        f2_T.open("data/T_integral.txt");
	for (k=0;k<nbas;k++)
	{
		for (i=0;i<nbas;i++)
		{
			for (l=0;l<bas[k*8+3];l++)
			{
				for (j=0;j<bas[i*8+3];j++)
				{
					int ind_re=ind_remap1e(i,j,k,l,N_con,bas);
					f2_S>>all_s_remap[ind_re];
					f2_V>>all_v_remap[ind_re];
					f2_T>>all_t_remap[ind_re];
				}
			}
		}
	}
	f2_S.close();
	f2_V.close();
	f2_T.close();
	ofstream f3_S,f3_V,f3_T;
	f3_S.open("data/S_integral.txt");
        f3_V.open("data/V_integral.txt");
        f3_T.open("data/T_integral.txt");
	for (i=0;i<N_con*N_con;i++)
	{
		f3_S<<setprecision(12)<<all_s_remap[i]<<' ';
		f3_V<<setprecision(12)<<all_v_remap[i]<<' ';
		f3_T<<setprecision(12)<<all_t_remap[i]<<' ';	
	}
	f3_S<<endl;
	f3_V<<endl;
	f3_T<<endl;
	f3_S.close();
	f3_V.close();
	f3_T.close();
	delete []atm;
	delete []bas;
	delete []env;
	delete []all_s_remap;
	delete []all_v_remap;
	delete []all_t_remap;
	return 0;
}
/*This version is transplatable if the bas.txt and the env.txt are set correctly!*/
/*In generated 1e and 2e files, shell orders are removed and contracted orbitals are sorted without shell. ex. (01) represents the 1st contracted orbital under shell 0, which we now index as 0.Max index=N_con-1.Thus the parameter files shall be written and intepreted in the same sequence that we have formulated.*/

int Cint_All2e()
{
	cout<<"_________2e Integration Subroutine_________"<<endl;
	ofstream f_2e;
	ifstream f_env;
	ifstream f_atm;
	ifstream f_bas;
	ifstream f_para;
	f_2e.open("data/2e_integral.txt");
	f_env.open("data/env.txt");
	f_atm.open("data/atom.txt");
	f_bas.open("data/bas.txt");
	f_para.open("data/para_c.txt");
	/*env array of the called cint function has been written according to cint manual in the info.txt file*/
        int i,j,k,l,p,q,r,s,natm,nbas,N_con,N_e;
	f_atm>>natm;
	cout<<"Number of atoms:"<<natm<<endl;
	f_bas>>nbas;
        cout<<"Number of shells:"<<nbas<<endl;
	f_para>>N_e>>N_con;
	cout<<"Number of electrons:"<<N_e<<endl;
	cout<<"Number of contracted gaussians:"<<N_con<<endl;
	int *atm=new int[natm*6];
	int *bas=new int[nbas*8];
	cout<<"Reading atomic information.atm:"<<endl;
        for (i=0;i<natm*6;i++)
	{
		f_atm>>atm[i];
		cout<<atm[i]<<' ';
	}
	cout<<endl;
	cout<<"Reading basis information,bas:"<<endl;
	for (i=0;i<nbas*8;i++)
	{
		f_bas>>bas[i];
		cout<<bas[i]<<' ';
	}
	cout<<endl;
	int size_env = natm*3;
	for (i=0;i<nbas;i++)
		size_env+=(bas[i*8+2]+bas[i*8+2]*bas[i*8+3]);
	double *env = new double[size_env];
	f_env>>natm>>nbas;
	for (i=0;i<size_env;i++)
		f_env>>env[i];
        int m,n;
	double *all_2e_remap= new double[N_con*N_con*N_con*N_con];
	for (m=0;m<nbas;m++)
	{
		for (n=0;n<nbas;n++)
		{
			for (p=0;p<nbas;p++)
			{
				for (q=0;q<nbas;q++)
				{
                                        int shls[4]={q,p,n,m};
					CINTOpt *opt=NULL;
	                                double *buf_2e =new double[bas[q*8+3]*bas[p*8+3]*bas[n*8+3]*bas[m*8+3]];
					cint2e_cart_optimizer(&opt,atm,natm,bas,nbas,env);
		                        cint2e_cart(buf_2e,shls,atm,natm,bas,nbas,env,opt);
                                        cout<<"2e integral "<<q<<','<<p<<','<<n<<','<<m<<"successful!:"<<endl;
		                        for (i=0;i<bas[q*8+3]*bas[p*8+3]*bas[n*8+3]*bas[m*8+3];i++)
		                        {
					        f_2e<<setprecision(12)<<buf_2e[i]<<' ';
			                        cout<<setprecision(12)<<buf_2e[i]<<' ';
		                        }
		                        cout<<endl;
		                        delete []buf_2e;
				}
			}
		}
	}
	cout<<"Setup files (atomic,basis,environmental) closed!"<<endl;
	f_2e<<endl;
	f_2e.close();
	f_atm.close();
	f_bas.close();
        f_env.close();
	f_para.close();
	ifstream f2_2e;
	f2_2e.open("data/2e_integral.txt");
	for (r=0;r<nbas;r++)
	{
		for (p=0;p<nbas;p++)
		{
			for (k=0;k<nbas;k++)
			{
				for (i=0;i<nbas;i++)
				{
					 for (s=0;s<bas[r*8+3];s++)
					 {
						for (q=0;q<bas[p*8+3];q++)
						 {
							 for (l=0;l<bas[k*8+3];l++)
							 {
								 for (j=0;j<bas[i*8+3];j++)
								 {
	                                 			        int ind_re=ind_remap2e(i,j,k,l,p,q,r,s,N_con,bas);
					                                f2_2e>>all_2e_remap[ind_re];
								 }  
							 }
						 }
					 }
				}
			}  
		}
        }
	f2_2e.close();
	ofstream f3_2e;
	f3_2e.open("data/2e_integral.txt");
	for (i=0;i<N_con*N_con*N_con*N_con;i++)
	{
		f3_2e<<setprecision(12)<<all_2e_remap[i]<<' ';
	}
	f3_2e<<endl;
	f3_2e.close();
	delete []atm;
	delete []bas;
	delete []env;
	delete []all_2e_remap;
	return 0;
}	

int Cint_S_Diag()
{
	cout<<"__________Basis Orthonormalization Subroutine__________"<<endl;
	ifstream f_S;
	ifstream f_para;
	ofstream f_U;
	ofstream f_X;
	f_S.open("data/S_integral.txt");
	f_para.open("data/para_c.txt");
	f_U.open("data/U_of_S.txt");
	f_X.open("data/X_transmat.txt");
	int i,j,N_e,N_con;
	f_para>>N_e>>N_con;
	cout<<"Number of electrons:"<<N_e<<endl;
	cout<<"Number of contracted gaussians:"<<N_con<<endl;
	double *S=new double[N_con*N_con];
	double *S_buf=new double[N_con*N_con];
	double *U=new double[N_con*N_con];
	double *d=new double[N_con];
	double *e=new double[N_con-1];
	double *tau=new double[N_con-1];
	double *X=new double[N_con*N_con];
	cout<<"Reading overlap integrals:"<<endl;
	for (i=0;i<N_con*N_con;i++)
	{
		f_S>>S[i];
		S_buf[i]=S[i];
	}
	if (LAPACKE_dsytrd(LAPACK_ROW_MAJOR,'U',N_con,S,N_con,d,e,tau)==0)
		cout<<"Tridiagonal trasformation successful!"<<endl;
	if (LAPACKE_dorgtr(LAPACK_ROW_MAJOR,'U',N_con,S,N_con,tau)==0)
		cout<<"Tridiagonal transformation matrix Q written to matrix S, original S buffered!"<<endl;
	if (LAPACKE_dsteqr(LAPACK_ROW_MAJOR,'V',N_con,d,e,S,N_con)==0)
		cout<<"Diagonal trasformation of S finished!Final U written to S!Eigenvalues written to d!"<<endl;
	for (i=0;i<N_con*N_con;i++)
		f_U<<setprecision(12)<<S[i]<<' ';
	f_U<<endl;
	cout<<"U written to file."<<endl;
	for (i=0;i<N_con;i++)
		for (j=0;j<N_con;j++)
			X[i*N_con+j]=S[i*N_con+j]/sqrt(d[j]);
	for (i=0;i<N_con*N_con;i++)
		f_X<<setprecision(12)<<X[i]<<' ';
	f_X<<endl;
	cout<<"X computed and successfully written!"<<endl;
	f_S.close();
	f_para.close();
	f_U.close();
	f_X.close();
        delete []S;
	delete []S_buf;
	delete []U;
	delete []d;
	delete []e;
	delete []tau;
	delete []X;
	return 0;
	/*Here, the U matrix is the same as in Szabo's book,p143!*/
}

int Cint_Nuc_Energy()
{
        cout<<"_________Nuclear Repulsion Subroutine__________"<<endl;
	ifstream f_atm;
	ifstream f_env;
	ifstream f_bas;
	ofstream f_nucene;
	f_atm.open("data/atom.txt");
	f_env.open("data/env.txt");
	f_bas.open("data/bas.txt");
	f_nucene.open("data/nucene.txt");
	int natm,nbas,i,j;
	f_atm>>natm;
	cout<<"Number of atoms:"<<natm<<endl;
	f_bas>>nbas;
	cout<<"Number of base shells:"<<nbas<<endl;
	int *atm=new int[natm*6];
	int *bas=new int[nbas*8];
	cout<<"Atom information:"<<endl;
        for (i=0;i<natm*6;i++)
	{
		f_atm>>atm[i];
		cout<<atm[i]<<' ';
	}
	cout<<endl;
	cout<<"Basis information:"<<endl;
	for (i=0;i<nbas*8;i++)
	{
		f_bas>>bas[i];
		cout<<bas[i]<<' ';
	}
	cout<<endl;
	f_env>>natm>>nbas;
	int size_env=natm*3;
	for (i=0;i<nbas;i++)
		size_env+=(bas[i*8+2]+bas[i*8+2]*bas[i*8+3]);
	double *env = new double[size_env];
	cout<<"Env information:"<<endl;
	for (i=0;i<size_env;i++)
	{
		f_env>>env[i];
		cout<<setprecision(12)<<env[i]<<' ';
	}
	cout<<endl;
	double nucene=0;
	for (i=0;i<natm;i++)
		for (j=i+1;j<natm;j++)
			nucene+=0.5*(atm[i*6+0]*atm[j*6+0])/sqrt(pow(env[atm[i*6+1]]-env[atm[j*6+1]],2)+pow(env[atm[i*6+1]+1]-env[atm[j*6+1]+1],2)+pow(env[atm[i*6+1]+2]-env[atm[j*6+1]+2],2));
	cout<<"Nuclear repulsion energy:"<<setprecision(12)<<nucene<<endl;
	f_nucene<<setprecision(12)<<nucene<<endl;
	f_atm.close();
	f_bas.close();
	f_env.close();
	f_nucene.close();
	delete []atm;
	delete []bas;
	delete []env;
	return 0;
}
/*End of OOParts integration library.*/
