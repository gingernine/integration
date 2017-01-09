#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;

double factorial(int l) {
	if (l==0) {
		return 1.0;
	} else if (l >= 1) {
		double prod = 1.0;
		for (int i=1; i <= l; i++) {
			prod *= 1.0 * i;
		}
		return prod;
	}
}

double legendre(double x, int n) {
	double sum = 0.0;
	double numer = 0.0;
	double denom = 0.0;
	for (int k=0; k <= n/2; k++) {
		numer = factorial(2*(n-k)) * pow(-1, k) * pow(x, n-2*k);
		denom = factorial(n-2*k) * factorial(k) * factorial(n-k);
		sum += numer/denom;
	}
	return sum / pow(2, n);
}

double legendre_deriv(double x, int n) {
	return n * (x * legendre(x, n) - legendre(x, n-1)) / (x*x - 1);
}

vector<double> newton_method(int n) {
	vector<double> root;
	double x, tmp;
	for (int i=0; i<n; i++) {
		x = cos((i+0.75) * M_PI / (n+0.5));
		tmp = 1;
		while (abs(x-tmp) >= 1e-10) {
			tmp = x;
			x = x - legendre(x, n) / legendre_deriv(x, n);
		}
		root.push_back(x);
	}
	return root;
}

double gauss_legendre_weight(double y, int n) {
    double denom = n*legendre(y, n-1);
    return 2*(1-y*y) / (denom*denom);
}

int main() {
    ofstream ofs("C:\\Users\\User\\Documents\\Cpp\\integration\\legendre_roots.csv");
    vector<double> roots;
    for (int n=1; n<=32; n++) {
        roots = newton_method(n);
    	for (int i=0; i<n; ++i) {
    		ofs << n << "," << roots[i] << "," << gauss_legendre_weight(roots[i], n) << endl;
    	}
    }
	return 0;
}
