#ifndef __VANILLA_OPTION_H
#define __VANILLA_OPTION_H

#include "payoff.h"

class VanillaOption {
public:
	PayOff* pay_off;

	double S0;
	double K;
	double r;
	double d;
	double T;
	double sigma;
	double down_barrier;
	double up_barrier;

	VanillaOption();
	VanillaOption(double _S0, double _K, double _r, double _d, double _T,
		double _sigma, PayOff* _pay_off, double _down_barrier, double _up_barrier);
};



#endif
