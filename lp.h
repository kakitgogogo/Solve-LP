#ifndef LP_H
#define LP_H

#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <memory>
#include <unordered_map>
#include <assert.h>
using namespace std;

#define DEBUG_MODE 1

class LP
{
public:
	enum 
	{
		MORE_ITERATION,
		MORE_ITERATION_DUAL,
		ONE_OPTIMAL_SOLUTION,
		INFINITE_OPTIMAL_SOLUTION,
		NO_OPTIMAL_SOLUTION,
	};

	enum
	{
		CHANGE_C,
		CHANGE_B,
		CHANGE_A,
		ADD_XN,
		ADD_RESTRAIN,
	};

	LP() {}

	LP(vector<double> param, vector<vector<double> > restrain, vector<int>base): C(param), Ab(restrain), XB(base), Ab_Origin(restrain), XB_Origin(base)
	{
		m = restrain.size(); 
		n = restrain[0].size();
		C.push_back(0);
		Sigma.resize(n);
		Theta.resize(m);

		nstep = 0;
		state = MORE_ITERATION;
	}

	int solve(vector<double>& res);

	int update(int change_mode, const vector<double>& vec);
	int update(int change_mode, const vector<vector<double> >& Anew);

private:
	int m, n;
	int nstep;
	int state;
	vector<double> C;
	vector<vector<double> > Ab;
	vector<int> XB;
	vector<double> Sigma;
	vector<double> Theta;

	vector<vector<double> > Ab_Origin;
	vector<int> XB_Origin;

	int cal_sigma();
	int cal_theta(int col);
	int iteration();
	int iteration_dual();

	int change_c(const vector<double>& c);
	int change_b(const vector<double>& b);
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

int LP::iteration()
{
	int row = 0, col = 0;

	++nstep;

#if DEBUG_MODE
	cout << "----------------- Intaration-" << nstep << " -----------------" << endl;
#endif

	state = cal_sigma();
	if(state != MORE_ITERATION) 
	{
#if DEBUG_MODE
	cout << "--------------- Intaration-Finish ---------------" << endl;
#endif		
		return state;
	}

	for(int i = 0; i < n-1; ++i)
	{
		if(Sigma[col] < Sigma[i])
		{
			col = i;
		}
	}

	state = cal_theta(col);
	if(state != MORE_ITERATION) 
	{
#if DEBUG_MODE
	cout << "--------------- Intaration-Finish ---------------" << endl;
#endif	
		return state;
	}

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

	double pivot = Ab[row][col];
	for(int i = 0; i < n; ++i)
	{
		Ab[row][i] /= pivot;
	}

	for(int i = 0; i < m; ++i)
	{
		if(i == row) continue;
		double multi = Ab[i][col];
		for(int j = 0; j < n; ++j)
		{
			Ab[i][j] -= multi * Ab[row][j];
		}
	}

#if DEBUG_MODE
	cout << "Ab: " << endl;
	for(int i = 0; i < m; ++i)
	{
		for(int j = 0; j < n; ++j)
		{
			cout << Ab[i][j] << " ";
		}
		cout << endl;
	}
#endif

	return MORE_ITERATION;
}

int LP::iteration_dual()
{
	int row = 0, col = 0;

	++nstep;

#if DEBUG_MODE
	cout << "----------------- Intaration-" << nstep << " -----------------" << endl;
#endif

	for(int i = 0; i < m; ++i)
	{
		if(Ab[i].back() < 0)
		{
			state = MORE_ITERATION;
			if(Ab[i].back() < Ab[row].back())
			{
				row = i;
			}
		}
	}
	if(state != MORE_ITERATION) 
	{
#if DEBUG_MODE
	cout << "--------------- Intaration-Finish ---------------" << endl;
#endif	
		return state;
	}

	double min_ratio = 0;
	for(int i = 0; i < n-1; ++i)
	{
		if(Ab[row][i] < 0 && Sigma[row]/Ab[row][i] < min_ratio)
		{
			col = i;
			min_ratio = Sigma[row]/Ab[row][i];
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

	double pivot = Ab[row][col];
	for(int i = 0; i < n; ++i)
	{
		Ab[row][i] /= pivot;
	}

	for(int i = 0; i < m; ++i)
	{
		if(i == row) continue;
		double multi = Ab[i][col];
		for(int j = 0; j < n; ++j)
		{
			Ab[i][j] -= multi * Ab[row][j];
		}
	}

	cal_sigma();

#if DEBUG_MODE
	cout << "Ab: " << endl;
	for(int i = 0; i < m; ++i)
	{
		for(int j = 0; j < n; ++j)
		{
			cout << Ab[i][j] << " ";
		}
		cout << endl;
	}
#endif

	return MORE_ITERATION;
}

int LP::solve(vector<double>& res)
{
	if(state == MORE_ITERATION) 
	{
		while((state = iteration()) == MORE_ITERATION);
	}
	else if(state == MORE_ITERATION_DUAL)
	{
		while((state = iteration_dual()) == MORE_ITERATION);
	}

	if(state == NO_OPTIMAL_SOLUTION) return state;

	res.resize(n);
	res[n-1] = -Sigma.back();
	for(int i = 0; i < XB.size(); ++i)
	{
		res[XB[i]] = Ab[i].back();
	}
	return state;
} 

/* -------------------- Sensitivity Analysis -------------------- */

int LP::change_c(const vector<double>& c)
{
	assert(c.size() == n);

	C = c;

	for(int i = 0; i < n; ++i)
	{
		double d = 0;
		for(int j = 0; j < m; ++j)
		{
			d += C[XB[j]] * Ab[j][i];
		}
		Sigma[i] = C[i] - d;
		if(Sigma[i] > 0) state = MORE_ITERATION;
	}

	return state;
}

int LP::change_b(const vector<double>& b)
{
	assert(b.size() == m);

	for(int i = 0; i < m; ++i)
	{
		Ab_Origin[i].back() = b[i];
	}

	for(int i = 0; i < m; ++i)
	{
		double d = 0;
		for(int j = 0; j < XB_Origin.size(); ++j)
		{
			d += Ab[i][XB_Origin[j]] * Ab_Origin[j].back();
		}
		Ab[i].back() = d;
		if(Ab[i].back() < 0) state = MORE_ITERATION_DUAL;
	}

	return state;
}

int LP::update(int change_mode, const vector<double>& vec)
{
	switch(change_mode)
	{
		case CHANGE_C:
			state = change_c(vec);
			break;
		case CHANGE_B:
			state = change_b(vec);
			break;
	}
	return state;
}


/* --------------------- Read Program File --------------------- */

static vector<string> get_line(FILE* fin)
{
	char str[500];
	vector<string> ret;

	if(fgets(str, 500, fin) == NULL)
	{
		return ret;
	}
	if(strcmp(str, "\n") == 0)
	{
		return ret;
	}

	if(str[0] == '#') return get_line(fin);

	char delimiter[] = " ;.,\n";
	char* item;
	ret.push_back(strtok(str, delimiter));
	while(item = strtok(NULL, delimiter))
		ret.push_back(item);

	return ret;
}


shared_ptr<LP> make_lp(string filename)
{
	FILE* fin = fopen(filename.c_str(), "r");
	if(fin == NULL)
	{
		perror("Open File Failed");
		return NULL;
	}

	int max_index = 0, new_var_idx = 0;
	double M, BigNum = (numeric_limits<double>::max)();
	vector<string> target; 
	vector<string> line;
	vector<double> param;
	vector<vector<double> > A;
	vector<int> b;
	vector<int> base;

	target= get_line(fin);
	for(int i = 3; i < target.size(); ++i)
	{
		size_t pos = 0;
		double c;

		if(isalpha(target[i][0])) c = 1;
		else if(target[i][0] == '+' && isalpha(target[i][1]))
		{
			c = 1;
			pos = 1;
		}
		else if(target[i][0] == '-' && isalpha(target[i][1]))
		{
			c = -1;
			pos = 1;
		}
		else c = stod(target[i], &pos);

		M = max(M, fabs(c));

		int idx = stod(target[i].substr(pos+1))-1;
		if(idx >= param.size())
		{
			param.resize(idx+1);
		}
		new_var_idx = max(new_var_idx, idx+1);
		param[idx] = c;
	}

	if(target[0] == "min")
	{
		for(int i = 0; i < param.size(); ++i)
		{
			param[i] = -param[i];
		}
	}

	unordered_map<string, int> RelationMap;
	RelationMap["="] = 0;
	RelationMap["<="] = RelationMap["<"] = 1;
	RelationMap[">="] = RelationMap[">"] = -1;

	while((line = get_line(fin)).empty());
	while(!line.empty())
	{
		int rel = 2;
		vector<double> Arow(new_var_idx, 0);
		for(int i = 0; i < line.size(); ++i)
		{
			switch(line[i][0])
			{
				case '=':
				case '<':
				case '>':
					if(rel != 2) return NULL;
					rel = RelationMap[line[i]];
					break;
				default:
					size_t pos = 0;
					double c;

					if(isalpha(line[i][0])) c = 1;
					else if(line[i][0] == '+' && isalpha(line[i][1]))
					{
						c = 1;
						pos = 1;
					}
					else if(line[i][0] == '-' && isalpha(line[i][1]))
					{
						c = -1;
						pos = 1;
					}
					else c = stod(line[i], &pos);

					M = max(M, fabs(c));

					if(pos == line[i].size())
					{
						if(i != line.size() - 1) return NULL;
						b.push_back(c);
						break;
					}
					int idx = stod(line[i].substr(pos+1))-1;
					if(idx >= Arow.size())
					{
						Arow.resize(idx+1);
					}
					new_var_idx = max(new_var_idx, idx+1);
					Arow[idx] = c;
					break;
			}
		}
		if(b.back() < 0)
		{
			b.back() = -b.back();
			rel = -rel;
			for(int i = 0; i < Arow.size(); ++i)
			{
				Arow[i] = -Arow[i];
			}
		}
		if(rel == 1)
		{
			Arow.resize(new_var_idx+1);
			Arow[new_var_idx] = 1;

			base.push_back(new_var_idx);

			param.resize(new_var_idx+1);
			param[new_var_idx++] = 0;
		}
		else if(rel == 0)
		{
			Arow.resize(new_var_idx+1);
			Arow[new_var_idx] = 1;

			base.push_back(new_var_idx);

			param.resize(new_var_idx+1);
			param[new_var_idx++] = BigNum;
		}
		else if(rel == -1)
		{
			Arow.resize(new_var_idx+2);
			Arow[new_var_idx] = -1;
			Arow[++new_var_idx] = 1;

			base.push_back(new_var_idx);

			param.resize(new_var_idx+1);
			param[new_var_idx++] = BigNum;
		}

		A.push_back(Arow);

		line = get_line(fin);
	}

	fclose(fin);

	M *= 100;

	for(int i = 0; i < param.size(); ++i)
	{
		if(param[i] == BigNum) 
			param[i] = -M;
	}

	for(int i = 0; i < A.size(); ++i)
	{
		for(int j = 0; j < A[i].size(); ++j)
		{
			A[i].resize(new_var_idx);
		}
		A[i].push_back(b[i]);
	}

#if DEBUG_MODE
	cout << "----------------- Init-Table -----------------" << endl;
	cout << "M: " << M << endl;
	cout << "target: ";
	for(int i = 0; i < param.size(); ++i)
	{
		cout << param[i] << " ";
	}
	cout << endl;
	cout << "A:" << endl;
	for(int i = 0; i < A.size(); ++i)
	{
		for(int j = 0; j < A[i].size(); ++j)
		{
			cout << A[i][j] << " ";
		}
		cout << endl;
	}
	cout << "base: " << endl;
	for(int i = 0; i < base.size(); ++i) 
	{
		cout << base[i] << " ";
	}
	cout << endl;
#endif

	return shared_ptr<LP>(new LP(param, A, base));
}

#endif