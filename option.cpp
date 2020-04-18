#ifndef __VANILLA_OPTION_CPP
#define __VANILLA_OPTION_CPP

#include "option.h"

VanillaOption::VanillaOption() {}

VanillaOption::VanillaOption(double _S0, double _K, double _r, double _d, double _T,
	double _sigma, PayOff* _pay_off, double _down_barrier, double _up_barrier) :
	S0(_S0), K(_K), r(_r), d(_d), T(_T), sigma(_sigma), pay_off(_pay_off) ,
	down_barrier(_down_barrier), up_barrier(_up_barrier){}

#endif