#include <iostream>
#include "lp.h"
using namespace std;

int main(int argc, char* argv[])
{
	if(argc < 2) 
	{
		fprintf(stderr, "Usage: %s programfile\n", argv[0]);
		return 1;
	}

	auto lp = make_lp(argv[1]);
	if(lp == NULL)
	{
		fprintf(stderr, "make_lp failed\n");
		return 1;
	}

	vector<double> res;
	int state = lp->solve(res);

	if(state == LP::INFINITE_OPTIMAL_SOLUTION)
	{
		cout << "INFINITE_OPTIMAL_SOLUTION" << endl;
		return 1;
	}
	if(state == LP::NO_OPTIMAL_SOLUTION)
	{
		cout << "NO_OPTIMAL_SOLUTION" << endl;
		return 1;
	}

	cout << "optimal solution: (";
	for(int i = 0; i < res.size() - 1; ++i)
	{
		cout << res[i] << ", ";
	}
	cout << ")" << endl;
	cout << "best result: " << res.back() << endl;
}