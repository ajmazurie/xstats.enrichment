#include <stdio.h>
#include <math.h>
#include <string.h>
#include <Python.h>

static PyObject *hypergeometric_distribution_ (PyObject *self, PyObject *args);
double hypergeometric_distribution (int i, int n, int B, int N);

static PyObject *mHG_pvalue_ (PyObject *self, PyObject *args);
long double mHG_pvalue (int B, int N, int max_size, double min_hgt);

////////////////////////////////////////////////////////////////////////////////

static PyObject *hypergeometric_distribution_ (PyObject *self, PyObject *args)
{
	int i, n, B, N;
	if (!PyArg_ParseTuple(args, "iiii", &i, &n, &B, &N))
		return NULL;
	
	double result = hypergeometric_distribution(i, n, B, N);

	return Py_BuildValue("d", result);
}

// Logarithm of the number of combinations of 'n' objects taken 'k' at a time
double ln_n_choose_k (int n, int k)
{
	return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1);
}

// Compute the hypergeometric distribution, or probability that a list of
// n objects should contain i ones with a particular property when the
// list has been selected randomly without replacement from a set of N
// objects in which B exhibit the same property
double hypergeometric_distribution (int i, int n, int B, int N)
{
	return exp(ln_n_choose_k(B, i) + ln_n_choose_k(N - B, n - i) - ln_n_choose_k(N, n));
}

/* For education purpose: original Python implementation
 
def __hypergeometric_probability (i, n, B, N):
	return exp(
	  __lncombination(B, i) +
	  __lncombination(N - B, n - i) -
	  __lncombination(N, n)
	)
 
def __lncombination (n, p):
	return \
	  __lnfactorial(n) - \
	  __lnfactorial(p) - \
	  __lnfactorial(n - p)
 
# Logarithm of n! with algorithmic approximation
# Reference:
#   Lanczos, C. 'A precision approximation of the gamma function',
#   J. SIAM Numer. Anal., B, 1, 86-96, 1964."
#   http://www.matforsk.no/ola/fisher.htm 
def __lnfactorial (n):
	if (n <= 1):
		return 0
	else:
		return __lngamma(n + 1)
 
def __lngamma (z):
	x = 0
	x += 0.1659470187408462e-06 / (z + 7)
	x += 0.9934937113930748e-05 / (z + 6)
	x -= 0.1385710331296526 / (z + 5)
	x += 12.50734324009056 / (z + 4)
	x -= 176.6150291498386 / (z + 3)
	x += 771.3234287757674 / (z + 2)
	x -= 1259.139216722289 / (z + 1)
	x += 676.5203681218835 / (z)
	x += 0.9999999999995183
 
	return log(x) - 5.58106146679532777 - z + (z - 0.5) * log(z + 6.5)
*/

////////////////////////////////////////////////////////////////////////////////

static PyObject *mHG_pvalue_ (PyObject *self, PyObject *args)
{
	int B, N, max_size;
	double min_hgt;
	if (!PyArg_ParseTuple(args, "iiid", &B, &N, &max_size, &min_hgt))
		return NULL;

	long double result = mHG_pvalue(B, N, max_size, min_hgt);

	char buffer[256];
	memset(buffer, 0, sizeof(char) * 256);
	sprintf(buffer, "%Lg", result);
	return Py_BuildValue("s", buffer);
}

// Taken from Calculate_pValue() from MotifTable.ccp
long double mHG_pvalue (int B, int N, int max_size, double min_hgt)
{
	int N_ = (N < max_size) ? N : max_size;
	int B_ = (B < max_size) ? B : max_size;

	// we allocate m on the heap rather than on the stack, to
	// prevent any memory issue in case B and/or N are large
	long double **m;
	m = (long double **)malloc(sizeof(long double *) * (B_ + 1));
	int i;
	for (i = 0; i < (B_ + 1); i++)
	{
		m[i] = (long double *)malloc(sizeof(long double) * (N_ + 1));
		memset(m[i], 0, sizeof(long double) * (N_ + 1));
	}
	m[0][0] = 1;
	
	long double base_hg = 1;
	
	int n;
	for (n = 1; n <= N_; n++)
	{
		int min_nB;
		if (B >= n)
		{
			min_nB = n;
			base_hg = base_hg * (B - n + 1) / (N - n + 1);
		}
		else
		{
			min_nB = B;
			base_hg = base_hg * n / (n - B);
		}
		
		long double tail_hg = base_hg;
		long double curr_hg = base_hg;

//		printf("%d	%g\n", n, base_hg);

		// first loop - sum up the tail, until the sum is bigger than min_hgt
		int b;
		for (b = min_nB; (tail_hg <= min_hgt) && (b > 0); b--)
		{
			m[b][n] = 0;

			curr_hg = curr_hg * (b * (N - B - n + b)) / ((n - b + 1) * (B - b + 1));
			tail_hg += curr_hg;
		}

		// second loop, starts when b is the maximal for which
		// HGT(N,B,n,b) > min_hgt
		for (; b > 0; b--)
		{
			m[b][n] = 0;

			// calculate current cell value by two optional
			// cells from which it can be reached

			// 1. last element in vector is 0
			if (m[b][n - 1] <= 1)
				m[b][n] += m[b][n - 1] * (N - B - n + b + 1) / (N - n + 1);

			// 2. last element in vector is 1
			if (m[b - 1][n - 1] <= 1)
				m[b][n] += m[b - 1][n - 1] * (B - b + 1) / (N - n + 1);
		}

		m[b][n] = 0;
		m[0][n] += m[0][n - 1] * (N - B - n + 1) / (N - n + 1);
	}
	
	long double r = 0;
	for (i = 0; i < B_ + 1; i++)
		r += m[i][N_];

	for (i = 0; i < (B_ + 1); i++)
	{
		free(m[i]);
		m[i] = NULL;
	}
	free(m);
	m = NULL;
	
	return (1.0 - r);
}

/* For education purpose: original Python implementation

def mHGT_ (s, max_size = 1000):
	B, N = len(filter(lambda x: x, s)), len(s)
	B_, N_ = min(B, max_size), min(N, max_size)

	# calculate the minimum HGT for the vector 's' of
	# successes, considering all possible partitions
	min_hgt = 1
	p = 0
	b = 0
	for n in range(1, N_ - 1):
		if (s[n-1]):
			b += 1

		pvalue = HGT(b, n, B, N)[1]
		if (pvalue <= min_hgt):
			min_hgt = pvalue
			p = n

	# calculate the corresponding exact p-value
	m = [[0 for j in range(N_ + 1)] for i in range(B_ + 1)]

	m[0][0] = 1.0
	base_hg = 1.0

	n = 1
	while (n <= N_):
		if (B >= n):
			min_nb = n
			base_hg *= (B - n + 1) / (N - n + 1.0)
		else:
			min_nb = B
			base_hg *= n / (n - B + 0.0)

		tail_hg = base_hg
		curr_hg = base_hg
		b = min_nb

		while (tail_hg <= min_hgt) and (b > 0):
			m[b][n] = 0.0

			curr_hg *= (b * (N - B - n + b)) / ((n - b + 1) * (B - b + 1.0))
			tail_hg += curr_hg

			b -= 1

		while (b > 0):
			m[b][n] = 0.0

			if (m[b][n-1] <= 1):
				m[b][n] += m[b][n-1] * (N - B - n + b + 1) / (N - n + 1.0)

			if (m[b-1][n-1] <= 1):
				m[b][n] += m[b-1][n-1] * (B - b + 1) / (N - n + 1.0)

			b -= 1

		m[0][n] = m[0][n-1] * (N - B - n + 1) / (N - n + 1.0)
		n += 1

		r = 0.0
		for i in range(B_ + 1):
			r += m[i][N_]

	return p, 1 - r
*/

////////////////////////////////////////////////////////////////////////////////

static PyMethodDef methods[] =
{
	{"hypergeometric_distribution", hypergeometric_distribution_, METH_VARARGS, "hypergeometric_distribution"},
	{"mHG_pvalue", mHG_pvalue_, METH_VARARGS, "mHG_pvalue"},
	{NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC initenrichment_ (void)
{
	(void)Py_InitModule("enrichment_", methods);
}

/*
main()
{
	mHG_pvalue(219, 1052, 1000, 1.29526774746e-07); // 5.70201200378e-06
	mHG_pvalue(50, 184, 1000, 2.47088282429e-07); // 3.83217156386e-06
	mHG_pvalue(32, 142, 1000, 9.71238952567e-08); // 9.4672248252e-07
	mHG_pvalue(293, 1271, 1000, 8.4386601215e-08); // 3.4188064657e-06
	mHG_pvalue(403, 1974, 1000, 1.1497812128e-07); // 2.26210630483e-06 !!
	mHG_pvalue(399, 1880, 1000, 1.69507675318e-07); // 3.35349159275e-06 !!
	mHG_pvalue(263, 1106, 1000, 5.0641563001e-07); // 2.08113077316e-05
	mHG_pvalue(17, 50, 1000, 1.89325091067e-06); // 4.20926993394e-06
	mHG_pvalue(57, 175, 1000, 7.55872628014e-06); // 0.000110935865582
	mHG_pvalue(2457, 4105, 1000, 0.000452071375635);
}
*/
