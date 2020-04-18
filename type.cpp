#ifndef __TYPE_CPP
#define __TYPE_CPP

#include "type.h"
#include <math.h>
#include <iostream>


Option_Type::Option_Type(int _N, int _J, VanillaOption* _option) : option(_option), N(_N), J(_J) {
	mu = option->r - option->d - 0.5*option->sigma*option->sigma;
	sigma = option->sigma;
	rf = option->r;
	d = option->d;
	T = option->T;
	K = option->K;
	down_barrier = option->down_barrier;
	up_barrier = option->up_barrier;
}

// Initial condition (vanilla call option)
double Option_Type::init_cond(double x) const {
	return option->pay_off->operator()(x);
}

// ==========
// Call, Put boundary
// ==========
Call::Call(int _N, int _J, VanillaOption* _option): Option_Type(_N, _J, _option) {};
double Call::boundary_Smin(double t, double s) const { return 0.0; }
double Call::boundary_Smax(double t, double s) const { return s*exp(-d*t) - K*exp(-rf*t); }

Put::Put(int _N, int _J, VanillaOption* _option): Option_Type(_N, _J, _option) {};
double Put::boundary_Smin(double t, double s) const { return K*exp(-rf*t); }
double Put::boundary_Smax(double t, double s) const { return 0.0; }

// ==========
// Eur Call, Put boundary
// ==========
Eur_Call::Eur_Call(int _N, int _J, VanillaOption* _option): Call(_N, _J, _option) {};
Eur_Put::Eur_Put(int _N, int _J, VanillaOption* _option): Put(_N, _J, _option) {};

// ==========
// Ame Call, Put boundary
// ==========
Ame_Call::Ame_Call(int _N, int _J, VanillaOption* _option): Call(_N, _J, _option) {};
void Ame_Call::boundary_EarlyExercise(VectorXd& new_vec, const VectorXd& vec0) const {
	new_vec = new_vec.binaryExpr(vec0, [](const double &a,const double &b) {return a>b? a:b;});
}
Ame_Put::Ame_Put(int _N, int _J, VanillaOption* _option): Put(_N, _J, _option) {};
void Ame_Put::boundary_EarlyExercise(VectorXd& new_vec, const VectorXd& vec0) const {
	new_vec = new_vec.binaryExpr(vec0, [](const double &a, const double &b) {return a>b? a:b;});
}

// ==========
// Ber Call, Put boundary
// ==========
Ber_Call::Ber_Call(int _N, int exercise_time_, VanillaOption* _option): Call(_N, exercise_time_*roundf(_option->T*20), _option) {
	exercise_time=exercise_time_;
	j=0;
};
void Ber_Call::boundary_Bermuda(VectorXd& new_vec, const VectorXd& vec0) {
	if (j % (int)roundf(T*20.0) == 0) new_vec = new_vec.binaryExpr(vec0, [](const double &a, const double &b) {return a>b? a:b;});
	++j;
}
Ber_Put::Ber_Put(int _N, int exercise_time_, VanillaOption* _option): Put(_N, exercise_time_*roundf(_option->T*20), _option) {
	j=0;
	exercise_time = exercise_time_;
};
void Ber_Put::boundary_Bermuda(VectorXd& new_vec, const VectorXd& vec0) {
	if (j % (int)roundf(T*20.0) == 0) new_vec = new_vec.binaryExpr(vec0, [](const double &a, const double &b) {return a>b? a:b;});
	++j;
}
#endif