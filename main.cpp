#include <iostream>

/*
	NUMERIC INTEGRATION
	******BEGIN
*/

template<typename Method, typename F,typename Float>
double numericIntegration(F function, Float a, Float b, int steps, Method m)
{
	double s = 0;
	double h = (b - a) / steps;
	for (int i = 0; i < steps; ++i){
		s += m(function, a + h * i, h);
	}
	return h * s;
}

template <typename Method, typename F, typename Float>
double numericDerivative(F function, Float a, Float h, Method m)
{
	return m(function, a, h);
}

class rect {

public:
	enum position_type {left, middle, right};
	rect(position_type pos) : position(pos) {}
	template<typename F, typename Float>
	double operator()(F f, Float x, Float h) const
	{
		switch (position)
		{
		case left:
			return f(x);
		case middle:
			return f(x + h / 2);
		case right:
			return f(x + h);
		}
		return 0;
	}
private:
	const position_type position;
};

class simpson
{
public:
	template<typename F, typename Float>
	double operator()(F f, Float x, Float h) const
	{
		return (f(x) + 4 * f(x + h / 2) + f(x + h)) / 6;
	}
};

class trapezium 
{
public:
	template<typename F, typename Float>
	double operator()(F f, Float x, Float h) const
	{
		return (f(x) + f(x + h)) / 2;
	}
};
/*
//simple usage

double f(double x) { return x * x * x;} // provided function

double rl = numericIntegration(f, 0.0, 1.0, 10, rect(rect::left));
double rm = numericIntegration(f, 0.0, 1.0, 10, rect(rect::middle));
double rr = numericIntegration(f, 0.0, 1.0, 10, rect(rect::right));
double t  = numericIntegration(f, 0.0, 1.0, 10, trapezium());
double s  = numericIntegration(f, 0.0, 1.0, 10, simpson());

	*****END
*/

/*
	NUMERIC Differentiation
	*****BEGIN
*/

class derivative 
{
public:
	enum derivative_type {backward, center, forward};
	derivative(derivative_type t) : d_type(t) {}
	template<typename F, typename Float>
	double operator()(F f, Float x, Float h) const
	{
		switch (d_type){
		case backward:
			return (f(x + h) - f(x)) / h;
		case center:
			return (f(x + 0.5 * h) - f(x - 0.5*h)) / h;
		case forward:
			return (f(x) - f(x - h)) / h;
		}
		return 0;
	}
private:
	const derivative_type d_type;
};

/*
this same pattern of usage
numericDerivative(f, 4.0, 0.001, derivative(derivative::center));
numericDerivative(f, 4.0, 0.001, derivative(derivative::backward));
numericDerivative(f, 4.0, 0.001, derivative(derivative::forward));

	*****END
*/
/*
	DISCRETE PID CONTROLLER
	*****BEGIN
*/

class DiscretePID { // backward Euler implementation of discrete PID controller
private:
	struct parameters {
		double Kp; // proportional gain
		double Ki; // integral gain
		double Kd; // derivative gain
		double Ts; //sampling time
		int N; // filter coefficients
		int UMAX; // max values
		int UMIN; // min values
	} properties;
	double error;
	double previous_error1;
	double previous_error2;
	double u0;
	double previous_u1;
	double previous_u2;
	double setpoint;
	double a0, a1, a2, b0, b1, b2, ku1, ku2, ke0, ke1, ke2; //some variables used in PID computation
public:
	DiscretePID(double Kp, double Ki, double Kd, double Ts, int N = 100, int UMAX = 2048, int UMIN = 0)
		: properties{ Kp, Ki, Kd, Ts, N, UMAX, UMIN },
		error(0),
		previous_error1(0),
		previous_error2(0),
		u0(0),
		previous_u1(0),
		previous_u2(0),
		setpoint(0),
		a0(0),
		a1(0),
		a2(0),
		b0(0),
		b1(0),
		b2(0),
		ku1(0),
		ku2(0),
		ke0(0),
		ke1(0),
		ke2(0)
	{
	}
	double setPoint(double point){
		setpoint = point;
	}

	double calculateOutput(double input){
		double N = static_cast<double>(properties.N);
		double Ts = properties.Ts;
		double Kp = properties.Kp;
		double Ki = properties.Ki;
		double Kd = properties.Kd;
		a0 = (1 + N * Ts);
		a1 = -(2 + N * Ts);
		a2 = 1;
		b0 = Kp * a0 + Ki * Ts * a0 + Kd * N;
		b1 = -(Kp * (-a1) + Ki * Ts + 2 * Kd * N);
		b2 = Kp + Kd * N;
		ku1 = a0 / a1;
		ku2 = a2 / a0;
		ke0 = b0 / a0;
		ke1 = b1 / a0;
		ke2 = b2 / a0;

		previous_error2 = previous_error1;
		previous_error1 = error;
		previous_u2 = previous_u1;
		previous_u1 = u0;
		error = setpoint - input;
		u0 = -ku1*previous_u1 - ku2 * previous_u2 + ke0 * error + ke1 * previous_error1 + ke2 * previous_error2;
		if (u0 > properties.UMAX)
			u0 = properties.UMAX;
		if (u0 < properties.UMIN)
			u0 = properties.UMIN;
		return u0;

	}


};


double f(double x) { return pow(x, 2) + x; };
int main()
{	
	/*
	DiscretePID PID(1, 1, 1, 0.02);
	PID.setPoint(4);
	PID.calculateOutput(50); //timer needed for cycle interruption
	*/
	std::cout << numericIntegration(f, 1.0, 10.5, 100, rect(rect::middle)) << std::endl;
	std::cout << numericDerivative(f, 3.0, 0.001, derivative(derivative::forward)) << std::endl;
	std::cout << "Press ENTER to continue...";
	std::cin.ignore(std::numeric_limits <std::streamsize> ::max(), '\n');
	return 0;
}