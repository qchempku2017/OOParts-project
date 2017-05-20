#include <time.h>
#include <iostream>
#include "scflib.h"
#include "intlib.h"
#include "diislib.h"

using namespace std;

int main()
{
        int i=0;
	double EIG_CONV,ENE_CONV,DIIS_CONV;
	cout<<"Welcome to OOParts-RHF, please specify eigenvalue convergence criteria(>=10E-10):"<<endl;
	cin>>EIG_CONV;
	cout<<"Specify energy convergence criteria(>=N_mo*10E-10):"<<endl;
	cin>>ENE_CONV;
	cout<<"Specify diis parameters convergence criteria(>10E-7)"<<endl;
	cin>>DIIS_CONV;
	time_t START,END;
	START = time(NULL);
	Cint_All1e();
	Cint_All2e();
        Cint_S_Diag();
	Cint_Nuc_Energy();
	cout<<"__________NOTICE: ENTERING SCF PROCESS!__________"<<endl;
	while(true)
	{
                cout<<"Iteration number:"<<i<<endl;
		cout<<endl;
		Cint_P_Forming();
		Cint_F_Forming();
         	Cint_F_Rescaling();
		bool eig_conv=Cint_F_Diag_And_C_New(EIG_CONV);
		bool ene_conv=Cint_Elec_Energy(ENE_CONV);
		bool diis_conv=Cint_DIIS(DIIS_CONV);
		if (eig_conv && ene_conv && diis_conv)
			/*You shouldn't write like if ((Cint_F_Diag_And_C_New())&&Cint_Elec_Energy()), because if you write this, when Cint_F_Diag...returns flase, Cint_Elec_Energy won't be called!*/
		{
			cout<<"All convergence requirements met!"<<endl;
			break;
		}
		else 
		{
			if(i>128)
			{
				cout<<"Failed converging to required condition(>_<)...Try another inital guess? I hope you know what you are doing."<<endl;
				break;
			}
		}
		i++;
	}
	Cint_Summary(i+1);
	END = time(NULL);
	cout<<"Total CPU time(s):"<<(END-START)<<endl;
	return 0;
}
