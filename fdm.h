#ifndef __FDM_H
#define __FDM_H

#include "type.h"
#include <vector>
#include <Eigen/Dense>

using namespace Eigen;

// Finite Difference Method - Abstract Base Class
class FDMBase {
protected:
	Option_Type* option;

	// Space discretisation
	double x_max;      // Spatial extent [0.0, x_max]
	double x_min;
	double S0, S0_down, S0_up;
	int S0_index;
	int N;			   // Number of spatial differencing points
	double dx;         // Spatial step size
	std::vector<double> x_values;  // Stores the log(S/S0) values
	std::vector<double> s_values;
	VectorXd indicator;
	double price;

	// Time discretisation
	double T;      // Temporal extent [0.0, t_dom]
	int J;		   // Number of time differencing points
	double dt;     // Step size 

	// Time-marching
	double tau;   // Time to maturity

	// Differencing coefficients
	double d, c, u;
	MatrixXd D;

	// Storage
	VectorXd new_result;   // New solution (becomes N+1)
	VectorXd old_result;   // Old solution (becomes N)
	VectorXd result_0;	   // Value at tau=0

	// Greeks
	double delta,gamma,theta;

	// Override these virtual methods in derived classes for 
	// specific FDM techniques, such as explicit Euler, Crank-Nicolson, etc.
	void calculate_step_sizes();
	void set_initial_conditions();
	void calculate_boundary_conditions();
	virtual void calculate_matrix_coef() = 0;
	void calculate_inner_domain();

public:
	// Constructor
	FDMBase(Option_Type* _option);
	virtual ~FDMBase() {};

	// Carry out the actual time-stepping
	void step_march();
	void greeks();
};

//****Explict****//
class FDMEulerExplicit : public FDMBase {
protected:
	void calculate_matrix_coef();

public:
	FDMEulerExplicit(Option_Type* _option);
	virtual ~FDMEulerExplicit() {};
};

//****Implict****//
class FDMEulerImplicit : public FDMBase {
protected:
	void calculate_matrix_coef();

public:
	FDMEulerImplicit(Option_Type* _option);
	virtual ~FDMEulerImplicit() {};
};

//****Crank-Nicolson****//
class FDMCN: public FDMBase {
protected:
	void calculate_matrix_coef();

public:
	FDMCN(Option_Type* _option);
	virtual ~FDMCN() {};
};
#endif 
