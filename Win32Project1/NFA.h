#pragma once

#define TABSIZE	100000
#define MAX_LUT_SIZE	1024

#ifndef TRUE
#define TRUE	1
#endif

#ifndef M_LN10
#define M_LN10	2.30258
#endif

#ifndef M_PI
#define M_PI	3.1415926
#endif

#define RELATIVE_ERROR_FACTOR	100.0

class NFALUT
{
public:
	NFALUT(int size, double _prob, double _logNT);
	~NFALUT();

	int *LUT;
	int LUTSize;
	double prob;
	double logNT;

	bool checkValidationByNFA(int n, int k);
	static double myAtan2(double yy, double xx);

private:
	double nfa(int n, int k);
	static double log_gamma_lanczos(double x);
	static double log_gamma_windschitl(double x);
	static double log_gamma(double x);
	static int double_equal(double a, double b);
};