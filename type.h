#ifndef __TYPE_H
#define __TYPE_H

#include "option.h"
#include <Eigen/Dense>

using namespace Eigen;

// Option type: Call/Put
// Exercise type: European/American/Bermuda
class Option_Type {
public:
	VanillaOption* option;
	double mu;
	double sigma;
	double rf;
	double d;
	double T;
	double K;
	double down_barrier;
	double up_barrier;
	int N;
	int J;

	Option_Type(int _N, int _J, VanillaOption* _option);
	virtual ~Option_Type() = 0 {};
	double init_cond(double x) const;
	virtual double boundary_Smin(double t, double s) const = 0;
	virtual double boundary_Smax(double t, double s) const = 0;
	virtual void boundary_EarlyExercise(VectorXd& new_vec, const VectorXd& vec0) const = 0;
	virtual void boundary_Bermuda(VectorXd& new_vec, const VectorXd& vec0) = 0;
};

//Call, Put
class Call:public Option_Type {
public:
	Call(int _N, int _J, VanillaOption* _option);
	virtual ~Call() {};
	double boundary_Smin(double t, double s) const;
	double boundary_Smax(double t, double s) const;
};
class Put:public Option_Type {
public:
	Put(int _N, int _J, VanillaOption* _option);
	virtual ~Put() {};
	double boundary_Smin(double t, double s) const;
	double boundary_Smax(double t, double s) const;
};

//Eur Call, Put
class Eur_Call:public Call {
public:
	Eur_Call(int N, int J, VanillaOption* _option);
	virtual ~Eur_Call() {};
	void boundary_EarlyExercise(VectorXd& new_vec, const VectorXd& vec0) const {};
	void boundary_Bermuda(VectorXd& new_vec, const VectorXd& vec0) {};
};
class Eur_Put:public Put {
public:
	Eur_Put(int N, int J, VanillaOption* _option);
	virtual ~Eur_Put() {};
	void boundary_EarlyExercise(VectorXd& new_vec, const VectorXd& vec0) const {};
	void boundary_Bermuda(VectorXd& new_vec, const VectorXd& vec0) {};
};

//Ame Call, Put
class Ame_Call:public Call {
public:
	Ame_Call(int N, int J, VanillaOption* _option);
	virtual ~Ame_Call() {};
	void boundary_EarlyExercise(VectorXd& new_vec, const VectorXd& vec0) const;
	void boundary_Bermuda(VectorXd& new_vec, const VectorXd& vec0) {};
};
class Ame_Put:public Put {
public:
	Ame_Put(int N, int J, VanillaOption* _option);
	virtual ~Ame_Put() {};
	void boundary_EarlyExercise(VectorXd& new_vec, const VectorXd& vec0) const;
	void boundary_Bermuda(VectorXd& new_vec, const VectorXd& vec0) {};
};

//Bermuda Call, Put
class Ber_Call:public Call {
	int exercise_time;
	int j;
public:
	Ber_Call(int _N, int exercise_time, VanillaOption* _option);
	virtual ~Ber_Call() {};
	void boundary_Bermuda(VectorXd& new_vec, const VectorXd& vec0);
	void boundary_EarlyExercise(VectorXd& new_vec, const VectorXd& vec0) const {};
};
class Ber_Put:public Put {
	int exercise_time;
	int j;
public:
	Ber_Put(int _N, int exercise_time, VanillaOption* _option);
	virtual ~Ber_Put() {};
	void boundary_Bermuda(VectorXd& new_vec, const VectorXd& vec0);
	void boundary_EarlyExercise(VectorXd& new_vec, const VectorXd& vec0) const {};
};


#endif
