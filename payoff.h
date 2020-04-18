#ifndef __PAY_OFF_HPP
#define __PAY_OFF_HPP

#include <algorithm> 

class PayOff {
public:
	PayOff();			  // Default constructor
	virtual ~PayOff() {}; // Virtual destructor

	// Overloaded () operator, turns the PayOff into an abstract function object
	virtual double operator() (const double& S) const = 0;
};

class PayOffCall : public PayOff {
private:
	double K; // Strike price

public:
	PayOffCall(const double& K_);
	virtual ~PayOffCall() {};

	// Virtual function is now over-ridden
	virtual double operator() (const double& S) const;
};

class PayOffPut : public PayOff {
private:
	double K; // Strike

public:
	PayOffPut(const double& K_);
	virtual ~PayOffPut() {};
	virtual double operator() (const double& S) const;
};

#endif

