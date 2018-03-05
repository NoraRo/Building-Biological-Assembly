#include<iostream>
#include<fstream>
#include<string>
#include <iomanip>
#include <vector>
#include <sstream>
#include<algorithm>
#include<map>
#include <filesystem>
#include "Functions.h"
using namespace std;
namespace fs = std::experimental::filesystem;

void generate_files(string virus)
{
	string part1 = "", part2 = "", biomits = "";
	string path = "E:/New Matrices/icosahedron_2btv/";
	string p1_file = "part1_" + virus + ".txt", p2_file = "part2_" + virus + ".txt", biomits_file = "full BIOMIT.txt", generated_file = "p";
	//reading part 1
	ifstream p(path + p1_file);
	if (!p.is_open()) {
		cout << "Error in file 1\n";  return;
	}
	part1.assign((std::istreambuf_iterator<char>(p)), (std::istreambuf_iterator<char>()));
	p.close();
	//reading part 2
	p.open(path + p2_file);
	if (!p.is_open()) {
		cout << "Error in file 2\n";  return;
	}
	part2.assign((std::istreambuf_iterator<char>(p)), (std::istreambuf_iterator<char>()));
	p.close();

	//reading biomits and creating new files
	p.open(path + biomits_file);
	if (!p.is_open()) {
		cout << "Error in file biomits file\n";  return;
	}
	string temp = "";
	int n = 0;
	ofstream of;
	while (getline(p, temp))
	{
		biomits = biomits + temp + "\n";//BIOMIT1
		getline(p, temp);
		biomits = biomits + temp + "\n";//BIOMIT2
		getline(p, temp);
		biomits = biomits + temp + "\n";//BIOMIT3

		//create new file and add its content
		generated_file = "p" + std::to_string(++n) + ".pdb";
		of.open(path + generated_file);

		if (!of.is_open()) {
			cout << "Error in file " + generated_file << endl;  return;
		}
		of << part1 << "\n" << biomits << part2;

		of.close();
	}
	p.close();
	cout << "End.\n";
}

void edit_file()
{
	string name_f = "E:/New Matrices/icosahedron_2btv/full BIOMIT.txt";
	ifstream rf(name_f);
	if (!rf.is_open()) {
		cout << "error in opeining file for reading \n";
		return;
	}

	ofstream wf("E:/New Matrices/icosahedron_2btv/full BIOMIT_new.txt");
	if (!wf.is_open()) {
		cout << "error in opeining file for wrting \n";
		return;
	}

	string tempS;
	double tempD;
	for (int i = 0, j = 1; i < 180; i++) {
		rf >> tempS;
		wf << tempS << " ";

		rf >> tempS;
		wf << tempS << "   ";

		rf >> tempS;
		wf << tempS << "   ";


		wf << std::setw(2) << j << "   ";

		for (size_t k = 0; k < 4; k++)
		{
			rf >> tempD;
			wf << std::setw(10) << fixed << setprecision(5) << tempD << "\t";
		}
		wf << endl;

		if ((i + 1) % 3 == 0) j++;
	}
	rf.close();
	wf.close();

}

vector<double**> read_smtry_records(ifstream & rf)
{
	vector<double**> list;
	const string REMARK = "REMARK";
	const string RENum = "290";
	const string SMTRY1 = "SMTRY1";
	const string SMTRY2 = "SMTRY2";
	const string SMTRY3 = "SMTRY3";
	string line = "", temp1, temp2;
	double ** matrix = nullptr;
	bool matrix_constructed = false, found_remark = false;
	int row = 0;
	double t = 0;
	while (!rf.eof()) {
		getline(rf, line);

		istringstream iss(line);
		iss >> temp1; // REMARK
		iss >> temp2; //290
		if (temp1.compare(REMARK) == 0 && temp2.compare(RENum) == 0)
		{

			iss >> temp1; //SMTRY1/SMTRY2/SMTRY3
			if (temp1.compare(SMTRY1) == 0)
			{
				found_remark = true;
				row = 0;
				matrix = new double*[4];
				for (int i = 0; i < 4; i++)
				{
					matrix[i] = new double[4];
				}
			}
			else if (temp1.compare(SMTRY3) == 0)
			{
				matrix_constructed = true;
				row = 2;
			}
			else if (temp1.compare(SMTRY2) == 0)
			{
				row = 1;
			}
			else if (found_remark)
				return list;
			else
				continue;


			iss >> temp1; // matrix number
			for (int i = 0; i < 4; i++)
			{
				iss >> t;
				matrix[row][i] = t;
			}

			if (matrix_constructed)
			{
				matrix[3][0] = matrix[3][1] = matrix[3][2] = 0;
				matrix[3][3] = 1;
				list.push_back(matrix);
				matrix_constructed = false;
			}
		}
	}


	return list;
}

vector<double**> read_MTRIX_records(ifstream & rf)
{
	vector<double**> list;
	const string MTRIX1 = "MTRIX1";
	const string MTRIX2 = "MTRIX2";
	const string MTRIX3 = "MTRIX3";
	string line = "", temp;
	double ** matrix = nullptr;
	bool matrix_constructed = false, found_remark = false;
	int row = 0;
	double t = 0;
	while (!rf.eof()) {
		getline(rf, line);

		istringstream iss(line);
		iss >> temp; // MTRIX1/MTRIX2/MTRIX3
		if (temp.compare(MTRIX1) == 0)
		{
			found_remark = true;
			row = 0;
			matrix = new double*[4];
			for (int i = 0; i < 4; i++)
			{
				matrix[i] = new double[4];
			}
		}
		else if (temp.compare(MTRIX3) == 0)
		{
			matrix_constructed = true;
			row = 2;
		}
		else if (temp.compare(MTRIX2) == 0)
		{
			row = 1;
		}
		else if (found_remark)
			return list;
		else
			continue;

		iss >> temp; // matrix number
		for (int i = 0; i < 4; i++)
		{
			iss >> t;
			matrix[row][i] = t;
		}

		if (matrix_constructed)
		{
			matrix[3][0] = matrix[3][1] = matrix[3][2] = 0;
			matrix[3][3] = 1;
			list.push_back(matrix);
			matrix_constructed = false;
		}
	}
	return list;
}

double ** multi(double** m, double ** l)
{
	double** r = new double*[4];
	for (int i = 0; i < 4; i++)
	{
		r[i] = new double[4];
	}

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			r[i][j] = 0;
		}
		for (int j = 0; j < 4; j++)
		{

			for (int k = 0; k < 4; k++)
			{
				r[i][j] += m[i][k] * l[k][j];
			}

		}
	}

	return r;
}

vector<double**> total_multiplication(const vector<double**> &slist, const vector<double**> &mlist, int size_s = 999, int size_m = 999)
{
	vector<double**> result;
	for (int i = 0 ; i<slist.size() && i< size_s; i++)
	{
		for (int j = 0; j<mlist.size() && j< size_m; j++)
		{
			auto r = multi(slist[i], mlist[j]);
			result.push_back(r);
		}

	}
	return result;
}

vector<double**> partial_multiplication(string file_name, const vector<double**> &slist, const vector<double**> &mlist)
{
	vector<double**> result;
	ifstream rf(file_name);
	if (!rf.is_open())
	{
		cout << "Partial multiplication is falied to open file \n";
		return result;
	}
	
	while (!rf.eof())
	{
		int sym = 0, mTRIX = 0;
		rf >> sym;
		for (auto m : mlist)
		{
			auto r = multi(slist[sym - 1], m);
			result.push_back(r);
		}
		
	}
	return result;
}

void rewrite_biomits(string path, string virus, string edited_name, const vector<double**> &new_biomits , int size =0)
{

	ifstream rf(path + virus);
	if (!rf.is_open()) {
		cout << "Error in original file\n";  return;
	}

	ofstream wf(path + edited_name + virus);
	if (!wf.is_open()) {
		cout << "Error in edited file\n";  return;
	}

	string line, temp1, temp2;
	const string REMARK = "REMARK";
	const string RENum = "350";
	int found_remark = 0;

	while (!rf.eof()) {
		getline(rf, line);

		istringstream iss(line);
		iss >> temp1; // REMARK
		iss >> temp2; //350
		if (temp1.compare(REMARK) == 0 && temp2.compare(RENum) == 0)
			found_remark++;
		else
			wf << line << endl;

		
		if (found_remark == 1)
		{
			for (int i = 0; i < new_biomits.size() && (size == 0 ? true : i < size); i++)
			{
				for (int j = 0; j < 3; j++)
				{
					wf << "REMARK 350   BIOMT"<<j+1<<"  ";
					wf << std::setw(2) << i+1 << "   ";
					for (int k = 0; k < 4; k++)
					{
						wf << std::setw(10) << fixed << setprecision(6) << new_biomits[i][j][k] << "\t";
					}
					wf << endl;
				}

			}
		}
	}
	wf.close();
}

void rewrite_symmtry(string path, string virus, string edited_name, const vector<double**> &new_sym)
{

	ifstream rf(path + virus);
	if (!rf.is_open()) {
		cout << "Error in original file\n";  return;
	}

	ofstream wf(path + edited_name + virus);
	if (!wf.is_open()) {
		cout << "Error in edited file\n";  return;
	}

	string line, temp1, temp2;
	const string REMARK = "REMARK";
	const string RENum = "290";
	int found_remark = 0;

	while (!rf.eof()) {
		getline(rf, line);

		istringstream iss(line);
		iss >> temp1; // REMARK
		iss >> temp2; //290
		if (temp1.compare(REMARK) == 0 && temp2.compare(RENum) == 0)
			found_remark++;
		else
			wf << line << endl;


		if (found_remark == 1)
		{
			for (int i = 0; i < new_sym.size(); i++)
			{
				for (int j = 0; j < 3; j++)
				{
					wf << "REMARK 290   SMTRY" << j + 1 << "  ";
					wf << std::setw(2) << i + 1 << "   ";
					for (int k = 0; k < 4; k++)
					{
						wf << std::setw(10) << fixed << setprecision(5) << new_sym[i][j][k] << "\t";
					}
					wf << endl;
				}

			}
		}
	}
	wf.close();
}

void sort_space_group()
{
	string name_f = "E:/New Matrices/report_spaceGroup.txt";
	ifstream rf(name_f);
	if (!rf.is_open()) {
		cout << "error in opeining file for reading \n";
		return;
	}

	ofstream wf("E:/New Matrices/report_Sorted_space_group.txt");
	if (!wf.is_open()) {
		cout << "error in opeining file for wrting \n";
		return;
	}
	vector<string> list;
	string temp;
	while (!rf.eof())
	{
		getline(rf, temp);
		list.push_back(temp);
	}
	sort(list.begin(), list.end());
	for (int i = 0; i < list.size(); i++)
	{
		wf << list[i] << endl;
	}
	rf.close();
	wf.close();
}

void count_space_group()
{
	string name_f = "E:/New Matrices/report_spaceGroup.txt";
	ifstream rf(name_f);
	if (!rf.is_open()) {
		cout << "error in opeining file for reading \n";
		return;
	}

	map<string, int> list;

	string temp;
	while (!rf.eof())
	{
		getline(rf, temp);
		if (list.find(temp) == list.end())
			list[temp] = 1;
		else list[temp]++;
	}
	//using it= map<string, int>::iterator ;
	//sort(list.begin(), list.end(), [](auto a, auto b) {return a->second > b->second; });
	for (auto x : list)
	{
		cout << x.first << "\t" << x.second << endl;
	}
	rf.close();
	
}

double ** generate_n_fold_Matrix(int n, char axis)
{
	double ** matrix = new double*[4];
	for (int i = 0; i < 4; i++)
	{
		matrix[i] = new double[4];
	}
	matrix[3][0] = 0; matrix[3][1] = 0;	  matrix[3][2] = 0;	matrix[3][3] = 1;

	double const PI = 3.14159265;
	double radius = (360.0 / n )* PI / 180.0;
	double COSINE = cos(radius), SINE = sin(radius);
	switch (axis)
	{
	case 'x':
		matrix[0][0] = 1;		matrix[0][1] = 0;	   matrix[0][2] = 0;		 matrix[0][3] = 0;
		matrix[1][0] = 0;		matrix[1][1] = COSINE; matrix[1][2] = -1 * SINE; matrix[1][3] = 0;
		matrix[2][0] = 0;		matrix[2][1] = SINE;   matrix[2][2] = COSINE;	 matrix[2][3] = 0;
		break;
	case 'y':
		matrix[0][0] = COSINE;  matrix[0][1] = 0;	   matrix[0][2] = SINE;		 matrix[0][3] = 0;
		matrix[1][0] = 0;		matrix[1][1] = 1;	   matrix[1][2] = 0;		 matrix[1][3] = 0;
		matrix[2][0] = -1 * SINE; matrix[2][1] = 0;	   matrix[2][2] = COSINE;	 matrix[2][3] = 0;
		break;
	case 'z':
		matrix[0][0] = COSINE;		matrix[0][1] = -1 * SINE; matrix[0][2] = 0;		 matrix[0][3] = 0;
		matrix[1][0] = SINE;		matrix[1][1] = COSINE; matrix[1][2] = 0;		 matrix[1][3] = 0;
		matrix[2][0] = 0;		matrix[2][1] = 0;   matrix[2][2] = 1;		 matrix[2][3] = 0;
		break;

	default:
		break;
	}

	return matrix;
}

vector<double**> generate_symtry(string file_name)
{
	vector<double**> list;

	ifstream rf(file_name);
	if (!rf.is_open()) {
		cout << "error in opeining file for reading \n";
		return list;
	}

	int n = 0, count = 0;
	char c = ' ';
	rf >> count;
	for (int i = 0; i < count; i++)
	{
		rf >> n >> c;
		if (n < 10) // one digit
			list.push_back(generate_n_fold_Matrix(n, c));
		//else if(n<99)//two digit
		//	list.push_back(generate_screw_Matrix(n, c));
		else
			cout << "error can't generate this fold " << n << " and " << c << endl;
	}

	return list;
}

void body()
{
	string name_f = "E:/New Matrices/", virus = "pdb1v9u.pdb";
	ifstream rf(name_f + virus);
	if (!rf.is_open()) {
		cout << "error in opeining file for reading \n";
		return ;
	}
	//vector<double**> list_SMTRY = generate_symtry("E:/New Matrices/fold_matrices.txt");
	vector<double**> list_SMTRY = read_smtry_records(rf);
	vector<double**> list_MTRIX = read_MTRIX_records(rf);
	vector<double**> list_answer = total_multiplication(list_SMTRY, list_MTRIX );
	//vector<double**> list_answer = partial_multiplication("E:/New Matrices/list_matrices.txt",list_SMTRY, list_MTRIX);
	rf.close();
	//rewrite_biomits(name_f,virus, "edited_sym_MTRIX", list_SMTRY);
	rewrite_biomits(name_f, virus, "edited_biomits", list_MTRIX, 60);
	//rewrite_symmtry(name_f, virus, "edited_REwrite_sym", list_SMTRY);
}

void body_1()
{
	string name_f = "E:/New Matrices/", virus = "pdb1auy.pdb";
	ifstream rf(name_f + virus);
	if (!rf.is_open()) {
		cout << "error in opeining file for reading \n";
		return;
	}
	vector<double**> list_SMTRY = generate_symtry("E:/New Matrices/fold_matrices.txt");
	//vector<double**> list_SMTRY = read_smtry_records(rf);
	vector<double**> list_MTRIX = read_MTRIX_records(rf);
	vector<double**> list_answer = total_multiplication(list_SMTRY, list_MTRIX);
	//vector<double**> list_answer = partial_multiplication("E:/New Matrices/list_matrices.txt",list_SMTRY, list_MTRIX);
	rf.close();
	//rewrite_biomits(name_f,virus, "edited_sym_MTRIX", list_SMTRY);
	rewrite_biomits(name_f, virus, "edited_biomits", list_SMTRY);
	//rewrite_symmtry(name_f, virus, "edited_REwrite_sym", list_SMTRY);
}


void body3()
{
	vector<double> cell = { 100.37,  151.33,  115.35,   90.00,  114.76,   90.00 };
	vector<double> cells(6);
	vector<vector<double>> deor , orth , deors, orths;

	lgg_crystal(cell, cells, deor, orth, deors, orths);
}

void preprocessing_input()
{
	
	string path ="E:/New Matrices/", mid_f = "Archive pdb/",
		ha_f = "HAs/", virus = "";

	
	for (auto& file : fs::directory_iterator(path + mid_f))
	{
		cout << file << endl;
		if (fs::path(file).extension() != fs::path(".pdb").extension())
			continue;


		ifstream rf(file);
		if (!rf.is_open()) {
			cout << file <<" error in opening file for reading \n";
			return;
		}

		string ha_path = path + ha_f + fs::path(file).stem().u8string() + ".ha";
		try {
			extract_HAs(rf, ha_path);
		}
		catch (exception e)
		{
			cout << e.what();
		}

		rf.close();

	}
}

int main()
{
	//edit_file();
	//generate_files("2btv");
	//body();
	//sort_space_group();
	//count_space_group();
	//body3();

	preprocessing_input();
	return 0;
}