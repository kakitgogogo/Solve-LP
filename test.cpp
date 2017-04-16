#include <iostream>
#include "lp.h"
using namespace std;

#define M 1000

void test_LP()
{
	cout << "==================== test-LP: ====================" << endl;
	vector<double> param{3, -1, -1, 0, 0, -M, -M};
	vector<vector<double> > restrain{{1, -2, 1, 1, 0, 0, 0, 11}, {-4, 1, 2, 0, -1, 1, 0, 3}, {-2, 0, 1, 0, 0, 0, 1, 1}};
	vector<int> base{3, 5, 6};
	LP lp(param, restrain, base);

	vector<double> res;
	int state = lp.solve(res);

	if(state == LP::INFINITE_OPTIMAL_SOLUTION)
	{
		cout << "INFINITE_OPTIMAL_SOLUTION" << endl;
		return;
	}
	if(state == LP::NO_OPTIMAL_SOLUTION)
	{
		cout << "NO_OPTIMAL_SOLUTION" << endl;
		return;
	}

	cout << "optimal solution: (";
	for(int i = 0; i < res.size() - 1; ++i)
	{
		cout << res[i] << ", ";
	}
	cout << ")" << endl;
	cout << "best result: " << res.back() << endl;
}

void test_make_lp()
{
	cout << "==================== test-make_lp: ====================" << endl;
	auto lp = make_lp("program.txt");

	vector<double> res;
	int state = lp->solve(res);

	if(state == LP::INFINITE_OPTIMAL_SOLUTION)
	{
		cout << "INFINITE_OPTIMAL_SOLUTION" << endl;
		return;
	}
	if(state == LP::NO_OPTIMAL_SOLUTION)
	{
		cout << "NO_OPTIMAL_SOLUTION" << endl;
		return;
	}

	cout << "optimal solution: (";
	for(int i = 0; i < res.size() - 1; ++i)
	{
		cout << res[i] << ", ";
	}
	cout << ")" << endl;
	cout << "best result: " << res.back() << endl;
}

int main()
{
	test_LP();
	test_make_lp();
}