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

double combination(int n, int i) {
	if (i==0) {
		return 1.0;
	} else {
		double prod = 1.0;
		for(int k=1; k<=i; ++k) {
			prod *= 1.0*(n-k+1)/k;
		}
		return prod;
	}
}

double laguerre(double x, int n) {
	double sum = 0.0;
	double comb;
	for (int i=0; i<=n; i++) {
		comb = combination(n, i);
		sum += pow(-x, i) * comb * comb * factorial(n-i);
	}
	return sum;
}

double laguerre_deriv(double x, int n) {
	return n * (laguerre(x, n) - n*laguerre(x, n-1)) / x;
}

vector<double> newton_method(int n) {
	vector<double> roots;
	double x, tmp;
	for (int i=0; i<n; i++) {
		x = M_PI*M_PI*(i+0.75)*(i+0.75) / (4*n);
		tmp = 1;
		while (abs(x-tmp) >= 1e-10) {
			tmp = x;
			x = x - laguerre(x, n) / laguerre_deriv(x, n);
		}
		roots.push_back(x);
	}
	return roots;
}

double gauss_laguerre_weight(double y, int n) {
    double l = factorial(n) / laguerre(y, n+1);
    return l*l*y;
}

int main() {
    ofstream ofs("C:\\Users\\User\\Documents\\Cpp\\integration\\laguerre_roots.csv");
    vector<double> roots;
    for (int n=1; n<=15; n++) {
        roots = newton_method(n);
    	for (int i=0; i<n; ++i) {
			cout << n << "," << roots[i] << "," << gauss_laguerre_weight(roots[i], n) << endl;
    		ofs << n << "," << roots[i] << "," << gauss_laguerre_weight(roots[i], n) << endl;
    	}
    }
	return 0;
}
