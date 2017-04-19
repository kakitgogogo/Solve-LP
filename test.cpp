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
	auto lp = make_lp("program/program0.txt");

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

void test_mix()
{
	cout << "==================== test-mix: ====================" << endl;
	cout << sizeof(LP) << endl;

	auto lp = make_lp("program/program0.txt");

	vector<double> res;
	int n = 10, state;
	while(n--) state = lp->solve(res);
	cout << "best result: " << res.back() << endl;
}

void test_update_after_solve()
{
	cout << "==================== test-update_after_solve: ====================" << endl;
	auto lp = make_lp("SensitivityAnalysis/program1.txt");

	vector<double> res;
	int state = lp->solve(res);
	cout << "best result: " << res.back() << endl;

	vector<double> c = {4, 2, 1, 0, 0, 0};
	state = lp->update(LP::CHANGE_C, c);
	cout << "state after update: " << state << endl;

	state = lp->solve(res);
	cout << "best result: " << res.back() << endl;

	vector<double> c1 = {-1, 2, 1, 0, 0, 0};
	state = lp->update(LP::CHANGE_C, c1);
	cout << "state after update: " << state << endl;

	state = lp->solve(res);
	cout << "best result: " << res.back() << endl;
}

void test_update_before_solve()
{
	cout << "==================== test-update_before_solve: ====================" << endl;
	auto lp = make_lp("SensitivityAnalysis/program1.txt");

	vector<double> res;
	vector<double> c = {4, 2, 1, 0, 0, 0};
	int state = lp->update(LP::CHANGE_C, c);
	cout << "state after update: " << state << endl;

	state = lp->solve(res);
	cout << "best result: " << res.back() << endl;
}

void test_change_b()
{
	cout << "==================== test-change_b: ====================" << endl;
	auto lp = make_lp("SensitivityAnalysis/program2.txt");

	vector<double> res;
	vector<double> b = {3, 2, 3};
	int state = lp->update(LP::CHANGE_B, b);

	state = lp->solve(res);
	cout << "best result: " << res.back() << endl;
}

void test_add_xn()
{
	cout << "==================== test-add_xn: ====================" << endl;
	auto lp = make_lp("SensitivityAnalysis/program3.txt");

	vector<double> res;
	int state = lp->solve(res);
	cout << "best result: " << res.back() << endl;

	vector<double> pnc = {3, 1, -3, 3};
	state = lp->update(LP::ADD_XN, pnc);

	state = lp->solve(res);

	cout << "state: " << state << endl;
	cout << "best result: " << res.back() << endl;
}

void test_add_restrain()
{
	cout << "==================== test-add_restrain: ====================" << endl;
	auto lp = make_lp("SensitivityAnalysis/program4.txt");

	vector<double> res;
	int state = lp->solve(res);
	cout << "best result: " << res.back() << endl;

	vector<double> an = {-3, 1, 6, 0, 0, 0, 1, 17};
	state = lp->update(LP::ADD_RESTRAIN, an);

	lp->showtab();

	state = lp->solve(res);

	cout << state << endl;
	cout << "best result: " << res.back() << endl;
}

int main()
{
	//test_LP();
	//test_make_lp();
	//test_mix();
	//test_update_after_solve();
	//test_update_before_solve();
	//test_change_b();
	test_add_xn();
	test_add_restrain();
}