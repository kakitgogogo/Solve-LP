#ifndef LP_H
#define LP_H

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <limits>
using namespace std;

#define DEBUG_MODE 0

class LP
{
public:
	enum 
	{
		MORE_ITERATION,
		ONE_OPTIMAL_SOLUTION,
		INFINITE_OPTIMAL_SOLUTION,
		NO_OPTIMAL_SOLUTION,
	};

	LP(vector<double> param, vector<vector<double> > restrain, vector<int>base): C(param), Ab(restrain), XB(base)
	{
		m = restrain.size(); 
		n = restrain[0].size();
		C.push_back(0);
		Sigma.resize(n);
		Theta.resize(m);
	}

	int solve(vector<double>& res);

private:

	int m, n;
	vector<double> C;
	vector<vector<double> > Ab;
	vector<int> XB;
	vector<double> Sigma;
	vector<double> Theta;

	int cal_sigma();
	int cal_theta(int col);
	int step();
};

int LP::cal_sigma()
{
	for(int i = 0; i < n; ++i)
	{
		Sigma[i] = 0;
		for(int j = 0; j < m; ++j)
		{
			Sigma[i] += C[XB[j]] * Ab[j][i];
		}
		Sigma[i] = C[i] - Sigma[i];
	}

#if DEBUG_MODE
	cout << "Sigma: ";
	for(int i = 0; i < n; ++i)
	{
		cout << Sigma[i] << " ";
	}
	cout << endl;
#endif

	int nzero = 0, i;

	for(i = 0; i < n - 1; ++i)
	{
		if(Sigma[i] > 0) return MORE_ITERATION;
		if(Sigma[i] == 0) nzero++;
	}
	
	return nzero > m ? INFINITE_OPTIMAL_SOLUTION : ONE_OPTIMAL_SOLUTION;
}

int LP::cal_theta(int col)
{
	bool has_optimal_solution = false;
	for(int i = 0; i < m; ++i)
	{
		if(Ab[i][col] <= 0) 
		{
			Theta[i] = (numeric_limits<double>::max)();
			continue;
		}
		Theta[i] = Ab[i].back() / Ab[i][col];
		has_optimal_solution = true;
	}

#if DEBUG_MODE
	cout << "Theta: ";
	for(int i = 0; i < m; ++i)
	{
		cout << Theta[i] << " ";
	}
	cout << endl;
#endif

	return has_optimal_solution ? MORE_ITERATION : NO_OPTIMAL_SOLUTION;
}

int LP::step()
{
	int state, row = 0, col = 0;

	state = cal_sigma();
	if(state != 0) return state;

	for(int i = 0; i < n-1; ++i)
	{
		if(Sigma[col] < Sigma[i])
		{
			col = i;
		}
	}

	state = cal_theta(col);
	if(state != 0) return state;

	for(int i = 0; i < m; ++i)
	{
		if(Theta[row] > Theta[i])
		{
			row = i;
		}
	}

	XB[row] = col;

#if DEBUG_MODE
	cout << "XB: ";
	for(int i = 0; i < m; ++i)
	{
		cout << XB[i] << " ";
	}
	cout << endl;
#endif

	int pivot = Ab[row][col];
	for(int i = 0; i < n; ++i)
	{
		Ab[row][i] /= pivot;
	}

	for(int i = 0; i < m; ++i)
	{
		if(i == row) continue;
		int multi = Ab[i][col];
		for(int j = 0; j < n; ++j)
		{
			Ab[i][j] -= multi * Ab[row][j];
		}
	}

	return MORE_ITERATION;
}

int LP::solve(vector<double>& res)
{
	int state;
	while((state = step()) == MORE_ITERATION);
	if(state != ONE_OPTIMAL_SOLUTION) return state;

	res.resize(n);
	res[n-1] = -Sigma.back();
	for(int i = 0; i < XB.size(); ++i)
	{
		res[XB[i]] = Ab[i].back();
	}
	return state;
} 

LP make_lp(string filename)
{
	
}

#endif