#include "payoff.h"
#include "option.h"
#include "type.h"
#include "fdm.h"
#include <iostream>
using namespace std;

int main(int argc, char **argv) {
	int optiontype, exestyle, barrier;
	double S0,K,r,v,T,d,down,up;
	int J,N,exercise_time;
	char choice = 'y';

	while (1) {
		cout<<"Choose option type (Call/ Put)\n";
		cin>>optiontype;
		cout<<"Choose exercise style (European/ American/ Bermuda)\n";
		cin>>exestyle;
		if (exestyle==3) {
			cout<<"Choose number of exercise time per year:\n";
			cin>>exercise_time;
		}

		cout<<"Input Stock price, K, risk-free rate, volatility, maturity, dividend:\n";
		cin>>S0>>K>>r>>v>>T>>d;
		cout<<"Choose Barrier (Up&Out / Down&Out / Double&Out / NA)\n";
		cin>>barrier;
		up = 5*S0;
		down = S0/5;
		if (barrier == 1) { cout<<"Input up barrier:\n"; cin>>up; }
		if (barrier == 2) { cout<<"Input down barrier:\n"; cin>>down; }
		if (barrier == 3) { cout<<"Input up&down barrier:\n"; cin>>up>>down; }

		cout<<"Input Spatial step N, time step J:\n";
		cin>>N>>J;
		
		PayOff* pay_off = &PayOffPut(K);
		VanillaOption* vanilla_option = &VanillaOption(S0, K, r, d, T, v, pay_off, down, up);
		Option_Type* option = &Eur_Call(N, J, vanilla_option);
		FDMBase* fdm_explicit;
		FDMBase* fdm_implicit;
		FDMBase* fdm_cn;
		if (optiontype==1) {
			pay_off = new PayOffCall(K);
			vanilla_option = new VanillaOption(S0, K, r, d, T, v, pay_off, down, up);
			if (exestyle == 1) option = new Eur_Call(N, J, vanilla_option);
			if (exestyle == 2) option = new Ame_Call(N, J, vanilla_option);
			if (exestyle == 3) option = new Ber_Call(N, J, vanilla_option);
		}
		else {
			pay_off = new PayOffPut(K);
			vanilla_option = new VanillaOption(S0, K, r, d, T, v, pay_off, down, up);
			if (exestyle == 1) option = new Eur_Put(N, J, vanilla_option);
			if (exestyle == 2) option = new Ame_Put(N, J, vanilla_option);
			if (exestyle == 3) option = new Ber_Put(N, J, vanilla_option);
		}

		// Run the FDM solver
		fdm_explicit = new FDMEulerExplicit(option);
		cout<<"\nExplicit\n";
		fdm_explicit->step_march();
		fdm_explicit->greeks();
		fdm_implicit = new FDMEulerImplicit(option);
		cout<<"\nImplicit\n";
		fdm_implicit->step_march();
		fdm_implicit->greeks();
		fdm_cn = new FDMCN(option);
		cout<<"\nCrank-Nicolson\n";
		fdm_cn->step_march();
		fdm_cn->greeks();

		// Delete the PDE, PayOff and Option objects
		delete fdm_cn;
		delete fdm_implicit;
		delete fdm_explicit;
		delete option;
		delete vanilla_option;
		delete pay_off;

		cout<<"Continue? (y/n)\n";
		cin>>choice;
		if(choice == 'n') break; 
	}
	
	cin.get();
	return 0;
}

