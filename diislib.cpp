/*OOParts DIIS library*/
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
extern"C"
{
#include "mkl.h"
}
#define min(a,b) (a<b)?a:b

using namespace std;

bool Cint_DIIS(double diis_conv)
{
	cout<<"__________OOParts DIIS subroutine__________"<<endl;
	bool DIIS_CONVERGE=true;
	ifstream f_para_former;
	ofstream f_para;
	f_para_former.open("data/para_c.txt");
	int N_e,N_con,i,j,k;
	f_para_former>>N_e>>N_con;
	cout<<"Number of electrons:"<<N_e<<endl;
	cout<<"Number of contracted gaussians:"<<N_con<<endl;
	vector <double> P_s;
	string line_buf;
	cout<<"Reading former C parameters!"<<endl;
	f_para_former.ignore();
	/*Use ignore here to prevent taking \n as an independent line!*/
	while(!f_para_former.eof()&&f_para_former.peek()!=EOF)
	{
		getline(f_para_former,line_buf);
	        istringstream line;
		line.str(line_buf);
	        double p;
        	for (i=0;i<N_con*N_con;i++)
	       	{
	        	line>>p;
	        	cout<<setprecision(12)<<p<<' ';
	                P_s.push_back(p);
        	}
	        cout<<endl;
	}
	cout<<"Former parameters read."<<endl;
	f_para_former.close();
	int n_l=P_s.size()/(N_con*N_con);
	cout<<"Current number of lines:"<<n_l<<endl;
	int m=min(n_l-1,5);
	double *B_pri=new double[(m+1)*(m+1)];
	cout<<"Constructing B inner product matrix of delta_p's:"<<endl;
	for (i=0;i<m;i++)
	{
		for (j=0;j<m;j++)
		{
			B_pri[i*(m+1)+j]=0;
			for (k=0;k<N_con*N_con;k++)
				B_pri[i*(m+1)+j]+=(P_s[(n_l-i-1)*N_con*N_con+k]-P_s[(n_l-i-2)*N_con*N_con+k])*(P_s[(n_l-j-1)*N_con*N_con+k]-P_s[(n_l-j-2)*N_con*N_con+k]);
			cout<<setprecision(12)<<B_pri[i*(m+1)+j]<<' ';
			/*The lines of P_s were read in a reversive manner.*/
		}
		cout<<endl;
	}
	cout<<"Setting redundant elements of B_pri."<<endl;
	for (i=0;i<m;i++)
	{
		B_pri[m*(m+1)+i]=-1;
		B_pri[i*(m+1)+m]=-1;
	}
	B_pri[m*(m+1)+m]=0;
	cout<<"B_pri successfully written in row major!"<<endl;
	double *X=new double[m+1];
	for(i=0;i<m;i++)
		X[i]=0;
	X[m]=-1;
	int *ipiv=new int[m+1];
        LAPACKE_dgetrf(LAPACK_ROW_MAJOR,m+1,m+1,B_pri,m+1,ipiv);
	cout<<"LU factorization of B_pri successful."<<endl;
	LAPACKE_dgetrs(LAPACK_ROW_MAJOR,'N',m+1,1,B_pri,m+1,ipiv,X,1);
	cout<<"DIIS linear functions solved.C's written to vector X:"<<endl;
	for (i=0;i<m;i++)
		cout<<setprecision(12)<<X[i]<<' ';
	cout<<endl;
	double *P_new=new double[N_con*N_con];
	cout<<"New parameters:"<<endl;
	f_para.open("data/para_c.txt");
	f_para<<N_e<<endl;
	f_para<<N_con<<endl;
	for (i=0;i<n_l-1;i++)
	{
		for (j=0;j<N_con*N_con;j++)
			f_para<<setprecision(12)<<P_s[i*N_con*N_con+j]<<' ';
		f_para<<endl;
	}
	for (i=0;i<N_con*N_con;i++)
	{
                P_new[i]=0;
		for (j=0;j<m;j++)
		{
			P_new[i]+=X[j]*P_s[(n_l-j-1)*N_con*N_con+i];
		}
		cout<<setprecision(12)<<P_new[i]<<' ';
		f_para<<setprecision(12)<<P_new[i]<<' ';
	}
	cout<<endl;
	f_para<<endl;
	for (i=0;i<N_con*N_con;i++)
		DIIS_CONVERGE = DIIS_CONVERGE && (abs(P_new[i]-P_s[(n_l-1)*N_con*N_con+i])<diis_conv);
	/* B_pri*C_pri=X */
	cout<<"DIIS converged?"<<DIIS_CONVERGE<<endl;
	f_para.close();
	delete []B_pri;
	delete []X;
	delete []ipiv;
	delete []P_new;
	return DIIS_CONVERGE;
}
/*This DIIS subroutine uses 5 previous paremeters vectors to generate a new DIIS guess. See algorithum on http://vergil.chemistry.gatech.edu/notes/diis/node2.html*/
