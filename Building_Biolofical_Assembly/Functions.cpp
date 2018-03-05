#include"Functions.h"
void lgg_crystal(vector<double>& cell,
	vector<double>& cells,
	vector<vector<double>>& deor,
	vector<vector<double>>& orth,
	vector<vector<double>>& deors,
	vector<vector<double>>& orths)
{
	double const PI = 3.14159265;
	vector<double> B1(3), B2(3), IUF(3);
	double AP = cell[3] * PI / 180.0;
	double BA = cell[4]* PI / 180.0;
	double GA = cell[5] * PI / 180.0;
	double coastar = (cos(BA)* cos(GA) - cos(AP)) / (sin(BA)*sin(GA));
	double siastar = sqrt(1. - coastar * coastar);
	for (int i = 0; i < 3; i++)
	{
		vector<double> arr;
		for (int j = 0; j < 3; j++)
		{
			arr.push_back(0);
		}
		orth.push_back(arr);
	}
	orth[0][0] = cell[0];
	orth[0][1] = cell[1] * cos(GA);
	orth[1][1] = cell[1] * sin(GA);
	orth[0][2] = cell[2] * cos(BA);
	orth[1][2] = -1 * cell[2] * sin(BA)*coastar;
	orth[2][2] = cell[2] * sin(BA) * siastar;

	deor = orth; // copy
	auto DE = determinanet(deor, 0.0000001);
	deors = transpose(orth);
	orths = transpose(deor);

	for (int i = 0; i < 3; i++)
		cells[i] = distancE(orths[i], 3)* 180 /PI;

	cells[3] = angle(orths[1], orths[2]) * 180.0 / PI;
	cells[4] = angle(orths[2], orths[0]) * 180.0 / PI;
	cells[5] = angle(orths[0], orths[1]) * 180.0 / PI;
}

//output is DE in return which is the determenaint.
// and A is the inverse of A.
double determinanet(vector<vector<double>>  &A, double EP)
{
	double DE = 1.0;
	int N = A.size();
	int k = 0 , j;
	vector<int> ME, C(N),B(N);
	//Initilize ME with rows numbers
	for (int i = 0; i < N; i++)
		ME.push_back(i);
	
	for (int i = 0; i < N; i++)
	{
		double y = 0.0;
		for (j = i; j < N; j++)
		{
			if (abs(A[i][j]) < abs(y)) continue;
				k = j;
				y = A[i][j];
		}
		DE = DE*y;
		if (abs(y) < EP)
		{
			cout << "Illconditioned Matrix:" << endl;
			for (int p = 0; p < N; p++)
			{
				for (int P = 0; P < N; P++)
				{
					cout << A[p][P] << " ";
				}cout << endl;
			}
			cout << "The pivot of the matrix is " << y << " N = " << N << endl;
			// Stop 4444 ??return;
		}
		y = 1.0 / y;
		for ( j = 0; j < N; j++)
		{
			C[j] = A[j][k];
			A[j][k] = A[j][i];
			A[j][i] = -1 * C[j] * y;
			B[j] = A[i][j] * y;
			A[i][j] = A[i][j] * y;
		}
		A[i][i] = y;
		j = ME[i];
		ME[i] = ME[k];
		ME[k] = j;
		for (k = 0; k < N; k++)
		{
			if (k == i) continue;
			for ( j = 0; j < N; j++)
			{
				if (j == i) continue;
				A[k][j] = A[k][j] - B[j] * C[k];
			}
		}

	}
	
	for (int i = 0; i < N; i++)
	{
		for ( k = 0; k < N; k++)
		{
			if (ME[k] == i) break;
		}
		if (k == i) continue;
		for ( j = 0; j < N; j++)
		{ //Swap k and i
			auto w = A[i][j];
			A[i][j] = A[k][j];
			A[k][j] = w;
		}
		auto IW = ME[i];
		ME[i] = ME[k];
		ME[k] = IW;
		DE = -1 * DE;
	}
	return DE;
}

vector<vector<double>> transpose(vector<vector<double>>& A)
{
	vector < vector<double>> trans(A[0].size(),
				vector<double>(A.size()));
	for (int i = 0; i < A.size(); i++)
	{
		for (int j = 0; j < A[i].size(); j++)
		{
			trans[j][i] = A[i][j];
		}

	}
	return trans;
}

double distancE(const vector<double>& v, int size)
{
	double c = 0;
	for (int i = 0; i < size; i++)
		c = c + v[i] * v[i];

	return sqrt(c);
}

double angle(vector<double>& v1, vector<double>& v2)
{
	auto a = distancE(v1, 3);
	auto b = distancE(v2, 3);
	auto c = dotproduct(v1, v2);
	auto arc = c / ( a* b);
	if (abs(arc) > 1.0)
	{
		if (abs(arc) - 1 > 0.00005) cout << "Warning arccod > 1 \n";
		if (arc > 1) arc = 1;
		if (arc < -1) arc = -1;

	}
	return acos(arc);
}

double dotproduct(vector<double>& v1, vector<double>& v2)
{
	auto res = 0;
	if (v1.size() != v2.size())
	{
		cout << "Two vectors do not have the same dimention \n";
		return res;
	}

	for (int i = 0; i < v1.size(); i++)
		res = v1[i] * v2[i] + res;

	return res;
}

vector<string> extract_HAs(istream& in, string ha_file)
{
	vector<string> HAs;
	ofstream of(ha_file);
	if (!of.is_open()) {
		cout << "error in opeining file for reading \n";
		return HAs;
	}

	string line;
	int count = 0;
	while (!in.eof())
	{
		++count;
		getline(in, line);
		if (line.size() == 0)
			continue;

		string x=  line.substr(76, 2), y = line.substr(0, 6);

		/*bool b = y.compare("ATOM")==0,
			c = is_HA(x),
			d = x.compare("S ") == 0;*/
		if (y.compare("CRYST1") == 0) {
			of << "cell " << line.substr(6, 49) << endl;
			of << "spacegroup " << line.substr(55, 10) << endl;
			continue;
		}

		if (y.compare("ATOM  " )==0 && is_HA(x))
		{
			of << line << endl;
			HAs.push_back(line);
		}
	}

	of.close();

	return HAs;
}

bool is_HA(string str)
{

	if (str.compare(" S")== 0 ||
		str.compare("NA")== 0 ||
		str.compare("MG")== 0 ||
		str.compare(" K")== 0 ||
		str.compare("CA")== 0 ||
		str.compare("MN")== 0 ||
		str.compare("FE")== 0 ||
		str.compare("CO")== 0 ||
		str.compare("NI")== 0 ||
		str.compare("CU")== 0 ||
		str.compare("ZN")== 0)
		return true;

	return false;
}

