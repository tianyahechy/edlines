#include "NFA.h"
#include <math.h>
#include <float.h>

NFALUT::NFALUT(int size, double _prob, double _logNT)
{
	LUTSize = size;
	LUT = new int[LUTSize];
	prob = _prob;
	logNT = _logNT;
	LUT[0] = 1;
	int j = 1;
	for (int i = 0; i < LUTSize; i++)
	{
		LUT[i] = LUTSize + 1;
		double ret = nfa(i, j);
		if (ret < 0)
		{
			while (j < i)
			{
				j++;
				ret = nfa(i, j);
				if ( ret >= 0 )
				{
					break;
				}
			}
			if ( ret < 0)
			{
				continue;
			}
		}
		LUT[i] = j;
	}
}

NFALUT::~NFALUT()
{
	delete[] LUT;
}

bool NFALUT::checkValidationByNFA(int n, int k)
{
	if ( n>= LUTSize)
	{
		return nfa(n, k) >= 0;
	}
	else
	{
		return k >= LUT[n];
	}
}

double NFALUT::myAtan2(double yy, double xx)
{
	static double LUT[MAX_LUT_SIZE + 1];
	static bool tableInited = false;
	if ( !tableInited )
	{
		for (int i = 0; i < MAX_LUT_SIZE; i++)
		{
			LUT[i] = atan((double)i / MAX_LUT_SIZE);
		}
		tableInited = true;

	}
	double y = fabs(yy);
	double x = fabs(xx);
	bool invert = false;
	if ( y > x )
	{
		double t = x;
		x = y;
		y = t;
		invert = true;
	}
	if ( x == 0)
	{
		x = 0.000001;
	}
	double ratio = y / x;
	double angle = LUT[(int)(ratio * MAX_LUT_SIZE)];
	if ( xx >= 0 )
	{
		if (yy >= 0)
		{
			if (invert)
			{
				angle = M_PI / 2 - angle;

			}
		}
		else
		{
			if (invert == false)
			{
				angle = M_PI - angle;
			}
			else
			{
				angle = M_PI / 2 + angle;
			}
		}
	}
	else
	{
		if (yy >= 0)
		{
			if (invert == false)
			{
				angle = M_PI - angle;
			}
			else
			{
				angle = M_PI / 2 + angle;
			}

		}
		else
		{
			if (invert)
			{
				angle = M_PI / 2 - angle;

			}
		}
	}
	return angle;
}

double NFALUT::nfa(int n, int k)
{
	static double inv[TABSIZE];
	double tolerance = 0.1;
	if (n < 0 || k < 0 || k > n || prob <= 0.0 || prob >= 1.0)
	{
		return -1.0;
	}
	if (n == 0 || k == 0)
	{
		return -logNT;
	}
	if ( n == k )
	{
		return -logNT - (double)n * log10(prob);
	}
	double p_term = prob / (1.0 - prob);
	double logItem = log_gamma((double)n + 1.0) - log_gamma((double)k + 1.0) - log_gamma((double)(n - k) + 1.0) + (double)k * log(prob) + double(n - k) * log(1.0 - prob);
	double term = exp(logItem);

	if (double_equal(term, 0.0 ))
	{
		if ((double)k > (double)n * prob)
		{
			return -logItem / M_LN10 - logNT;
		}
		else
		{
			return -logNT;
		}
	}

	double bin_tail = term;
	for (int i = k+1; i <= n; i++)
	{
		double bin_term = (double)(n - i + 1) * (i < TABSIZE ? (inv[i] != 0.0 ? inv[i] : (inv[i] = 1.0 / (double)i)) : 1.0 / (double)i);
		double mult_term = bin_term * p_term;
		term *= mult_term;
		bin_tail += term;
	}

}