#pragma once
#include<vector>
#include<iostream>
#include <string>
#include <fstream>
#include <iomanip>
using namespace std;

void lgg_crystal(vector<double>& cell,
	vector<double>& cells,
	vector<vector<double>>& deor,
	vector<vector<double>>& orth,
	vector<vector<double>>& deors,
	vector<vector<double>>& orths);

double determinanet(vector<vector<double>>  &A, double EP);
vector<vector<double>> transpose(vector<vector<double>>& A);
double distancE(const vector<double>& v, int size);
double angle(vector<double>& v1, vector<double>& v2);
double dotproduct(vector<double>& v1, vector<double>& v2);
//This function checks if this atom is HA or not according to this paper
//Metals in protein structures: a review of their principal features by Harding 2010
//http://dx.doi.org/10.1080/0889311X.2010.485616 
bool is_HA(string str);
//This function takes file handle and extracts all HAs.
// returns list of strings of HAs and writes them on a separete file
vector<string> extract_HAs(istream& in, string ha_file);
//Searching for the min and max values of X, Y and Z coordinates 
vector<double> range_XYZ(istream& in, string out_f);
//Matrix multiplication 
vector<double> vec_multi_matrix(vector<vector<double>> m, vector<double> vec);
vector<vector<double>> identity_vec(int size);
double** identity_i(int size);
vector<vector<double>> identity_i_vec(int size);
