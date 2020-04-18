#ifndef __FDM_CPP
#define __FDM_CPP

#include <iostream>
#include <algorithm>
#include <vector>
#include "fdm.h"
using namespace std;

FDMBase::FDMBase(Option_Type* _option)
	: option(_option) {
	N = option->N;
	N /= 6;
	N *= 6;
	J = option->J;
	T = option->T;
	S0 = option->option->S0;
	if (option->up_barrier >= 5 * S0) x_max = 5 * option->sigma*pow(T, 0.5);
	else x_max = 1.2*log(option->up_barrier / S0);
	if (option->down_barrier <= S0 / 5) x_min = -5 * option->sigma*pow(T, 0.5);
	else x_min = 1.2*log(option->down_barrier / S0);
	
	//cout << S0 * exp(x_max) << endl;
	//cout << option->down_barrier << ' ' << option->up_barrier << endl;
	calculate_step_sizes();
	set_initial_conditions();
}

void FDMBase::calculate_step_sizes() {
	dx = (x_max - x_min) / N;
	dt = T / J;
	// find S0 index on grid in [S0_index, S0_index+1)
	S0_index = int((0 - x_min) / dx);
	
}

void FDMBase::set_initial_conditions() {
	// Spatial settings
	double cur_x = 0.0;
	double cur_s = 0.0;

	old_result.resize(N+1);
	new_result.resize(N+1);
	x_values.resize(N+1, 0.0);
	s_values.resize(N+1, 0.0);
	indicator.resize(N + 1);

	for (int n = 0; n <= N; ++n) {
		cur_x = x_min + n * dx;
		cur_s = S0*exp(cur_x);
		new_result(n) = option->init_cond(cur_s);
		result_0 = new_result;
		x_values[n] = cur_x;
		s_values[n] = cur_s;
		if (cur_s >= option->up_barrier || cur_s <= option->down_barrier)
			indicator(n) = 0.0;
		else indicator(n) = 1.0;
	}

	// Temporal settings
	tau = 0.0;
}

void FDMBase::calculate_boundary_conditions() {
	new_result(0) = option->boundary_Smin(tau, s_values[0]);
	new_result(N) = option->boundary_Smax(tau, s_values[N]);
	option->boundary_EarlyExercise(new_result,result_0);
	option->boundary_Bermuda(new_result, result_0);
}

void FDMBase::calculate_inner_domain() {
	new_result = D*new_result;
}

void FDMBase::step_march() {
	for (int j = 0; j<J; ++j) {
		calculate_inner_domain();
		calculate_boundary_conditions();
		new_result = new_result.cwiseProduct(indicator);
		tau += dt;
		if (j == J-2) old_result = new_result;
	}
	
	S0_down = S0 * exp(x_values[S0_index]);
	S0_up = S0 * exp(x_values[S0_index + 1]);
	price = (new_result(S0_index + 1) - new_result(S0_index)) / (S0_up - S0_down)*(S0 - S0_down) + new_result(S0_index);
	//price = (new_result(S0_index + 1) - new_result(S0_index)) / dx * (0 - x_values[S0_index]) + new_result(S0_index);
	cout << "S0: is between " << S0_down << " and " << S0_up << endl
		<< "price: " << price << endl;

}
void FDMBase::greeks() {
	delta = (new_result(S0_index + 1) - new_result(S0_index)) / dx / S0;
	gamma = (new_result(S0_index + 1) - 2 * new_result(S0_index) + new_result(S0_index - 1))
		/ pow(dx, 2) / pow(S0, 2) - delta / S0;
	double old_price = (old_result(S0_index + 1)- old_result(S0_index)) / (S0_up - S0_down)*(S0 - S0_down) 
		+ old_result(S0_index);
	theta = (old_price - price) / dt /365;

	cout<<"\Delta, Gamma, Theta:\n";
	std::vector<double> greek{delta,gamma,theta};
	for_each(greek.begin(), greek.end(), [](double &i) {std::cout<<i<<endl;});
}

//****Explicit****//
FDMEulerExplicit::FDMEulerExplicit(Option_Type* _option) : FDMBase(_option) {
	calculate_matrix_coef();
	//set_initial_conditions();
}

void FDMEulerExplicit::calculate_matrix_coef() {
	double h = dt / dx;
	double k = dt / pow(dx, 2);
	double sigma = option->sigma;
	double mu = option->mu;
	double r = option->rf;
	d = 0.5*(k*sigma*sigma - h * mu);
	c = 1 - (k*sigma*sigma + r * dt);
	u = 0.5*(k*sigma*sigma + h * mu);
	D = MatrixXd::Zero(N + 1, N + 1);
	for (int i = 0; i <= N; i++)
	{
		D(i, i) = c;
		if (i > 0) D(i - 1, i) = u;
		if (i < N) D(i + 1, i) = d;
	}
	
}
//****Implicit****//
FDMEulerImplicit::FDMEulerImplicit(Option_Type* _option): FDMBase(_option) {
	calculate_matrix_coef();
	//set_initial_conditions();
}

void FDMEulerImplicit::calculate_matrix_coef() {
	double h = dt / dx;
	double k = dt / pow(dx, 2);
	double sigma = option->sigma;
	double mu = option->mu;
	double r = option->rf;
	d = -0.5*(k*sigma*sigma - h * mu);
	c = 1 + (k*sigma*sigma + r * dt);
	u = -0.5*(k*sigma*sigma + h * mu);
	D = MatrixXd::Zero(N+1, N+1);
	for (int i = 0; i <= N; i++)
	{
		D(i, i) = c;
		if (i > 0) D(i - 1, i) = u;
		if (i < N) D(i + 1, i) = d;
	}
	D = D.inverse();
}

//****Crank-Nicolson****//
FDMCN::FDMCN(Option_Type* _option): FDMBase(_option) {
	calculate_matrix_coef();
	//set_initial_conditions();
}

void FDMCN::calculate_matrix_coef() {
	double h = dt / dx;
	double k = dt / pow(dx, 2);
	double sigma = option->sigma;
	double mu = option->mu;
	double r = option->rf;
	double c_;
	d = -0.25*(k*sigma*sigma - h*mu);
	c = 1.0 + k*sigma*sigma/2;
	u = -0.25*(k*sigma*sigma + h*mu);
	c_ = c-k*sigma*sigma-r*dt;
	D = MatrixXd::Zero(N+1, N+1);
	MatrixXd D_ = MatrixXd::Zero(N+1, N+1);
	for (int i = 0; i <= N; i++)
	{
		D(i, i) = c;
		D_(i, i) = c_;
		if (i > 0) {
			D(i - 1, i) = u;
			D_(i - 1, i) = -u;
		}
		if (i < N) {
			D(i + 1, i) = d;
			D_(i + 1, i) = -d;
		}
	}
	D = D.inverse()*D_;
}

#endif