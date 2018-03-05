#pragma once
#include<vector>
#include<iostream>
#include <string>
#include <fstream>
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
