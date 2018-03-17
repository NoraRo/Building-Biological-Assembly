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

vector<vector<vector<double>>> read_smtry_records_vec(ifstream & rf)
{
	vector<vector<vector<double>>> list;
	const string REMARK = "REMARK";
	const string RENum = "290";
	const string SMTRY1 = "SMTRY1";
	const string SMTRY2 = "SMTRY2";
	const string SMTRY3 = "SMTRY3";
	string line = "", temp1, temp2;
	vector<vector<double>> matrix (4);
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
				for (int i = 0; i < 4; i++)
				{
					vector<double> mm(4);
					matrix[i] = mm;
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

vector<vector<vector<double>>> read_BIOMT_records_vec(ifstream & rf)
{
	vector<vector<vector<double>>> list;
	const string REMARK = "REMARK";
	const string RENum = "350";
	const string SMTRY1 = "BIOMT1";
	const string SMTRY2 = "BIOMT2";
	const string SMTRY3 = "BIOMT3";
	string line = "", temp1, temp2;
	vector<vector<double>> matrix(4);
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
				for (int i = 0; i < 4; i++)
				{
					vector<double> mm(4);
					matrix[i] = mm;
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

vector<vector<vector<double>>> read_MTRIX_records_vec(ifstream & rf)
{
	vector<vector<vector<double>>> list;
	const string MTRIX1 = "MTRIX1";
	const string MTRIX2 = "MTRIX2";
	const string MTRIX3 = "MTRIX3";
	string line = "", temp;
	vector<vector<double>> matrix(4);
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
			
			for (int i = 0; i < 4; i++)
			{
				vector<double> mm(4);
				matrix[i] = mm;
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


vector<vector<double>> read_ATOM_records(ifstream & rf, vector<vector<string>> & labels)
{
	const string ATOM = "ATOM";
	string line = "";
	vector<vector<double>> matrix;
	double t;

	int row = 0;
	while (!rf.eof()) {
		getline(rf, line);
		if (line.size() == 0)
			continue;
		istringstream iss(line);
		ostringstream oss;
		string temp = "", str = "";
		iss >> temp; //ATOM
		
		if (temp.compare(ATOM) == 0)
		{
			vector<double> record;
			//passing to XYZ parts
			vector<string> l(2);
			//oss <<  setw(6) <<left<< ATOM ;			
			//iss >> str;
			//oss << setw(4) <<right<< str;
			//iss >> str;
			//oss << setw(1) << " ";
			//oss << setw(3) <<left << str;
			//iss >> str;
			//oss << setw(1) << str;
			//iss >> str;
			//oss << setw(2) <<right<< str;
			//iss >> str;
			//oss << setw(2) <<" ";
			//oss << setw(3) <<right<< str;
			char* w = new char[27];
			iss.read(w, 26);
			w[26] = '\0';
			l[0] = ATOM + w;
			delete[]w;
			//for (int j = 0; j < 5; j++)
			//{
			//	iss >> str;
			//	oss << setw(6) << str;
			//	
			//}
			//l[0] = oss.str();
			/*iss >> t;
			iss >> str;
			iss >> str;
			iss >> t;
			iss >> t;*/
			//reading XYZ
			for (int j = 0; j < 3; j++)
			{
				iss >> t;
				record.push_back(t);
			}
			//reading/saving the rest of the line
			/*for (int j = 0; j < 2; j++)
			{
				iss >> str;
				l[1] = l[1] + str + " ";
			}*/
			w = new char[26];
			iss.read(w, 25);
			w[25] = '\0';
			l[1] = w;
			delete[]w;
			labels.push_back(l);
			//convert it to 4x1 vector
			record.push_back(1);
			matrix.push_back(record);
			row++;
		}
		
	}
	return matrix;
}

void writing_ATOMS_only(ostream&of, vector<vector<double>> new_ATOMs, 
			vector<vector<string>> lables)
{
	
	for (int i = 0; i < new_ATOMs.size(); i++)
	{
		of << lables[i][0];
		for (int j = 0; j < new_ATOMs[0].size() - 1; j++)
		{
			of << std::setw(8) << fixed << setprecision(3)
				<< new_ATOMs[i][j];
		}
		of << lables[i][1] << endl;
	}
}
//This function extract rotation from ncs.ofm file
//That is the output of FINDNCS program
//one based
vector<vector<double>> read_ROTATION_records(ifstream & rf, int ncs_number)
{
	const string label = ".LSQ_RT_NCS" + to_string(ncs_number);
	string line = "", temp;
	vector<vector<double>> matrix;
	double t;
	
	while (!rf.eof()) {
		getline(rf, line);

		istringstream iss(line);
		iss >> temp; //ATOM

		if (temp.compare(label) == 0 )
		{
			//Initializing matrix
			//Last record is 0 0 0 1
			for (int j = 0; j < 4; j++)
			{
				vector<double> record(4);
				matrix.push_back(record);
			}

			for (int j = 0; j < 4 && !rf.eof(); j++)
			{
				//reading the next line
				getline(rf, line);
				istringstream iss2(line);
				//inserting the last colomn ==> T
				if (j == 3)
				{
					for (int k = 0; k < 3; k++)
					{
						iss2 >> t;
						matrix[k][3] = t;
					}
					break;
				}
				//reading XYZ ==> R
				for (int k = 0; k < 3; k++)
				{
					iss2 >> t;
					matrix[j][k] = t;
				}
			}
			//Last row filling it to 4x4 matrix
			matrix[3][0] = 0;
			matrix[3][1] = 0;
			matrix[3][2] = 0;
			matrix[3][3] = 1;
			break;
		}

	}
	return matrix;
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

vector<vector<double>> multi_vec(vector<vector<double>> m, vector<vector<double>> l)
{
	vector<vector<double>> r ;
	for (int i = 0; i < m.size(); i++)
	{
		vector<double> t(m[0].size());
		r.push_back(t);
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

vector<double**> total_multiplication( vector<double**> &slist,  vector<double**> &mlist, int size_s = 999, int size_m = 999)
{
	//Handling zero matrix
	if (slist.size() == 0)
	{
		auto x = identity_i(4);
		slist.push_back(x);
	}
	if (mlist.size() == 0)
	{
		auto x = identity_i(4);
		mlist.push_back(x);
	}

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

//If you need the new files in the same place then pass the newPath as originalPath
void rewrite_biomits(string originalPath, string newPath, string virus, string edited_name, const vector<double**> &new_biomits , int size =0)
{

	ifstream rf(originalPath + virus);
	if (!rf.is_open()) {
		cout << "Error in original file\n";  return;
	}

	ofstream wf(newPath + edited_name + virus);
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
	string name_f = "E:/New Matrices/", virus = "pdb1a34.pdb";
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
	//rewrite_biomits(name_f, name_f,virus, "edited_sym_MTRIX", list_SMTRY);
	rewrite_biomits(name_f, name_f, virus, "edited_biomits", list_MTRIX);
	//rewrite_symmtry(name_f, name_f, virus, "edited_REwrite_sym", list_SMTRY);
}

//Regenerate BIOMITs from multiplying MTRIX with SYM
void body_1()
{
	string name_f = "E:/New Matrices/", mid_f = "Archive pdb/", virus = "pdb1bbt.pdb";
	string r_path = name_f + mid_f;
	string w_path = name_f + "Try/";
	ifstream rf(r_path + virus);
	if (!rf.is_open()) {
		cout << "error in opeining file for reading \n";
		return;
	}
	//vector<double**> list_SMTRY = generate_symtry("E:/New Matrices/fold_matrices.txt");
    vector<double**> list_SMTRY = read_smtry_records(rf);
	vector<double**> list_MTRIX = read_MTRIX_records(rf);
	vector<double**> list_answer = total_multiplication(list_SMTRY, list_MTRIX);
	//vector<double**> list_answer = partial_multiplication("E:/New Matrices/list_matrices.txt",list_SMTRY, list_MTRIX);
	rf.close();
	//rewrite_biomits(name_f, name_f,virus, "edited_sym_MTRIX", list_SMTRY);
	rewrite_biomits(r_path, w_path, virus, "edited_biomits", list_answer);
	//rewrite_symmtry(name_f, name_f, virus, "edited_REwrite_sym", list_SMTRY);
}

//Reads part of Atoms and RT in NCS and generates the new Atoms
void body3()
{
	string name_f = "E:/New Matrices/Fortran/out/", virus = "part_1bev.pdb";
	ifstream rf(name_f + virus);
	if (!rf.is_open()) {
		cout << "error in opeining file for reading \n";
		return;
	}
	
	vector<vector<string>> lables;
	vector<vector<double>> list_ATOMs = read_ATOM_records(rf, lables);
	rf.close();

	int ncsMax = 10;
	for (int k = 0; k < ncsMax; k++)
	{


		string file = "ncs.ofm";
		ifstream rf2(name_f + file);
		if (!rf2.is_open()) {
			cout << "error in opeining file for reading \n";
			return;
		}

		//Choose which NCS
		int NCS_num = k+1;
		vector<vector<double>> rotation = read_ROTATION_records(rf2, NCS_num);
		rf2.close();

		//Start logic
		vector<vector<double>> new_ATOMs;
		for (int i = 0; i < list_ATOMs.size(); i++)
		{
			new_ATOMs.push_back(vec_multi_matrix(rotation, list_ATOMs[i]));
		}

		//Write them to file
		ofstream of(name_f + "1bev-FINDNCS/new_"+ to_string(NCS_num)+"_" + virus);
		if (!of.is_open()) {
			cout << "error in opeining file for reading \n";
			return;
		}
		//should be a function
		for (int i = 0; i < new_ATOMs.size(); i++)
		{
			of << lables[i][0];
			for (int j = 0; j < new_ATOMs[0].size(); j++)
			{
				of << std::setw(8) << fixed << setprecision(3)
					<< new_ATOMs[i][j] << "  ";
			}
			of << lables[i][1]<< endl;
		}
		of.close();
	}
}

void body4()
{
	vector<double> cell = { 100.37,  151.33,  115.35,   90.00,  114.76,   90.00 };
	vector<double> cells(6);
	vector<vector<double>> deor , orth , deors, orths;

	lgg_crystal(cell, cells, deor, orth, deors, orths);

	vector<vector<double>> sxyz = { { 193.831,  12.734, 136.923 } 
		,{ 221.393,   8.258, 125.614 }
		,{ 196.435,  12.976, 147.965 }
		,{ 190.474, -5.181, 136.633	 }
		,{ 201.054, -26.972, 150.251 }
		,{ 197.978, -21.072, 123.994 }	
	};
	vector<vector<double>> sfxyz;
	
	for (int i = 0; i < sxyz.size(); i++)
	{
		sfxyz.push_back(vec_multi_matrix(deor, sxyz[i]));
	}
	
	string name_f = "E:/New Matrices/", virus = "receprocal_1bev.pdb";
	ofstream of(name_f + virus);
	if (!of.is_open()) {
		cout << "error in opeining file for reading \n";
		return;
	}
	
	for (int i = 0; i < sfxyz.size(); i++)
	{
		for (int j = 0; j < sfxyz[i].size(); j++)
		{
			of << sfxyz[i][j] << "\t";
		}
		of << endl;
	}

	of.close();
}

//Generate all possible BIOMITS 
void body_11()
{
	string name_f = "E:/New Matrices/", mid_f = "Archive pdb/", virus = "";
	string originalPath = name_f + mid_f;
	string newPath = name_f + "Report all BIOMITs/";
	
	for (auto& file : fs::directory_iterator(originalPath))
	{
		cout << file << endl;
		if (fs::path(file).extension() != fs::path(".pdb").extension())
			continue;


		ifstream rf(file);
		if (!rf.is_open()) {
			cout << file << " error in opening file for reading \n";
			return;
		}
		virus = fs::path(file).stem().u8string() + ".pdb";


		//vector<double**> list_SMTRY = generate_symtry("E:/New Matrices/fold_matrices.txt");
		vector<double**> list_SMTRY = read_smtry_records(rf);
		vector<double**> list_MTRIX = read_MTRIX_records(rf);
		vector<double**> list_answer = total_multiplication(list_SMTRY, list_MTRIX);
		//vector<double**> list_answer = partial_multiplication("E:/New Matrices/list_matrices.txt",list_SMTRY, list_MTRIX);
		rf.close();
		//rewrite_biomits(name_f, name_f,virus, "edited_sym_MTRIX", list_SMTRY);
		rewrite_biomits(originalPath, newPath, virus, "edited_biomits", list_answer);
		//rewrite_symmtry(name_f, name_f, virus, "edited_REwrite_sym", list_SMTRY);
	}
}
//Generating report
void body0()
{
	string path = "E:/New Matrices/", mid_f = "Archive pdb/",
		new_f = "Report MTRIX as BIOMITs/", virus = "";

	string originalPath = path + mid_f;
	string newPath = path + new_f;
	vector<vector<int>> repo;
	vector<string> viruses;
	for (auto& file : fs::directory_iterator(originalPath))
	{
		vector<int> r;
		cout << file << endl;
		if (fs::path(file).extension() != fs::path(".pdb").extension())
			continue;


		ifstream rf(file);
		if (!rf.is_open()) {
			cout << file << " error in opening file for reading \n";
			return;
		}
		virus = fs::path(file).stem().u8string() + ".pdb";
		vector<double**> list_SYM = read_smtry_records(rf);
		vector<double**> list_MTRIX = read_MTRIX_records(rf);
		rf.close();

		int total = 0;
		if (list_MTRIX.size() == 0)
			total = list_SYM.size();
		else if (list_SYM.size() == 0)
			total = list_MTRIX.size();
		else total = list_MTRIX.size() * list_SYM.size();

		viruses.push_back(virus);
		r.push_back(list_MTRIX.size());
		r.push_back(list_SYM.size());
		r.push_back(total);
		repo.push_back(r);
	}

	string file = "report MTRIXandSYM.csv";
	ofstream of(newPath + file);
	if (!of.is_open()) {
		cout << file << " error in opening file for reading \n";
		return;
	}
	of << setw(10) << "VIRUS" << "," << setw(6) << "MTRIX" << ","
		<< setw(4) << "SYM" << "," << setw(5) << "TOTAL" << endl;
	for (int i = 0; i < viruses.size(); i++)
	{
		of << setw(10) << viruses[i] << "," << setw(6) << repo[i][0]
			<< "," << setw(4) << repo[i][1] << ","
			<< setw(5) << repo[i][2] << endl;
	}
	of.close();
}

//writing MTRIX as BIOMITs
void body5()
{
	string path = "E:/New Matrices/", mid_f = "Archive pdb/",
		new_f = "Report MTRIX as BIOMITs/", virus = "";

	string originalPath = path + mid_f;
	string newPath = path + new_f;

	for (auto& file : fs::directory_iterator(originalPath))
	{
		cout << file << endl;
		if (fs::path(file).extension() != fs::path(".pdb").extension())
			continue;


		ifstream rf(file);
		if (!rf.is_open()) {
			cout << file << " error in opening file for reading \n";
			return;
		}
		virus = fs::path(file).stem().u8string() + ".pdb";
		
		vector<double**> list_MTRIX = read_MTRIX_records(rf);
		rf.close();


		rewrite_biomits(originalPath, newPath, virus, "edited_biomits", list_MTRIX, 60);
		

	}
}
//Reads the whole atoms in pdb file and generate new atoms after applying NCS
void body6() {

		string name_f = "E:/New Matrices/Archive pdb/", virus = "pdb1bev.pdb";
		ifstream rf(name_f + virus);
		if (!rf.is_open()) {
			cout << "error in opeining file for reading \n";
			return;
		}

		vector<vector<string>> lables;
		vector<vector<double>> list_ATOMs = read_ATOM_records(rf, lables);
		rf.close();

		int ncsMax = 10;
		for (int k = 0; k < ncsMax; k++)
		{


			string file = "ncs.ofm";
			name_f = "E:/New Matrices/Fortran/out/";
			ifstream rf2(name_f + file);
			if (!rf2.is_open()) {
				cout << "error in opeining file for reading \n";
				return;
			}

			//Choose which NCS
			int NCS_num = k + 1;
			vector<vector<double>> rotation = read_ROTATION_records(rf2, NCS_num);
			rf2.close();

			//Start logic
			vector<vector<double>> new_ATOMs;
			for (int i = 0; i < list_ATOMs.size(); i++)
			{
				new_ATOMs.push_back(vec_multi_matrix(rotation, list_ATOMs[i]));
			}

			//Write them to file
			name_f = "E:/New Matrices/NCS vs MTRIX/";
			string ncsORmtrix = "ncs";
			string writing_f = name_f + "new_" + ncsORmtrix + to_string(NCS_num) + "_" + virus;
			ofstream of(writing_f);
			if (!of.is_open()) {
				cout << "error in opeining file for reading \n";
				return;
			}
			//should be a function
			cout << writing_f << endl;
			writing_ATOMS_only(of, new_ATOMs, lables);
			of.close();
		}
	}

//Reads the whole atoms in pdb file and generate new atoms after applying MTRIX
void body7() {

	string name_f = "E:/New Matrices/Archive pdb/", virus = "pdb1bev.pdb";
	ifstream rf(name_f + virus);
	if (!rf.is_open()) {
		cout << "error in opeining file for reading \n";
		return;
	}

	vector<vector<string>> lables;
	vector<vector<double>> list_ATOMs = read_ATOM_records(rf, lables);
	rf.close();

	ifstream rf1(name_f + virus);
	if (!rf1.is_open()) {
		cout << "error in opeining file for reading \n";
		return;
	}
	vector<vector<vector<double>>> rotation = read_MTRIX_records_vec(rf1);
	rf1.close();

	
	for (int k = 0; k < rotation.size(); k++)
	{

		//Start logic
		vector<vector<double>> new_ATOMs;
		for (int i = 0; i < list_ATOMs.size(); i++)
		{
			new_ATOMs.push_back(vec_multi_matrix(rotation[k], list_ATOMs[i]));
		}
		
		//Write them to file
		name_f = "E:/New Matrices/NCS vs MTRIX/";
		string ncsORmtrix = "MTRIX";
		string writing_f = name_f + "new_" + ncsORmtrix + to_string(k+1) + "_" + virus;
		ofstream of(writing_f);
		if (!of.is_open()) {
			cout << "error in opeining file for reading \n";
			return;
		}
		//should be a function THE CORRECT 1
		cout << writing_f << endl;
		writing_ATOMS_only(of, new_ATOMs, lables);
		of.close();
	}
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

void extract_ranges()
{

	string name_f = "E:/New Matrices/", virus = "pdb1bev.pdb";
	ifstream in(name_f + virus);
	if (!in.is_open()) {
		cout << "error in opeining file for reading \n";
		return ;
	}

	vector<double> XYZ = range_XYZ(in, name_f+"range_xyx_1bev.txt");

	in.close();
}
//Reads the whole atoms in pdb file and generate new atoms after applying SYM
void body8() {
	string name_f = "E:/New Matrices/Archive pdb/", virus = "pdb1bev.pdb";
	//string name_f = "E:/New Matrices/NCS vs MTRIX/", virus = "new_MTRIX1_pdb1bev.pdb";
	ifstream rf(name_f + virus);
	if (!rf.is_open()) {
		cout << "error in opeining file for reading \n";
		return;
	}

	vector<vector<string>> lables;
	vector<vector<double>> list_ATOMs = read_ATOM_records(rf, lables);
	rf.close();

	//name_f = "E:/New Matrices/Archive pdb/"; virus = "pdb1bev.pdb";
	ifstream rf1(name_f + virus);
	if (!rf1.is_open()) {
		cout << "error in opeining file for reading \n";
		return;
	}
	vector<vector<vector<double>>> rotation = read_smtry_records_vec(rf1);
	rf1.close();


	for (int k = 0; k < rotation.size(); k++)
	{

		//Start logic
		vector<vector<double>> new_ATOMs;
		for (int i = 0; i < list_ATOMs.size(); i++)
		{
			new_ATOMs.push_back(vec_multi_matrix(rotation[k], list_ATOMs[i]));
		}

		//Write them to file
		name_f = "E:/New Matrices/NCS vs MTRIX/";
		string ncsORmtrix = "SYM";
		string writing_f = name_f + "new_" + ncsORmtrix + to_string(k + 1) + "_" + virus;
		ofstream of(writing_f);
		if (!of.is_open()) {
			cout << "error in opeining file for reading \n";
			return;
		}
		//should be a function
		cout << writing_f << endl;
		writing_ATOMS_only(of, new_ATOMs, lables);
		of.close();
	}

}

//Reads the whole atoms in pdb file and generate new atoms after applying 
//All possible BIOMITs
void body10() {
	string name_f = "E:/New Matrices/Report all BIOMITs/", virus = "edited_biomitspdb1bev.pdb";
	//string name_f = "E:/New Matrices/NCS vs MTRIX/", virus = "new_MTRIX1_pdb1bev.pdb";
	ifstream rf(name_f + virus);
	if (!rf.is_open()) {
		cout << "error in opeining file for reading \n";
		return;
	}

	vector<vector<string>> lables;
	vector<vector<double>> list_ATOMs = read_ATOM_records(rf, lables);
	rf.close();

	//name_f = "E:/New Matrices/Archive pdb/"; virus = "pdb1bev.pdb";
	ifstream rf1(name_f + virus);
	if (!rf1.is_open()) {
		cout << "error in opeining file for reading \n";
		return;
	}
	vector<vector<vector<double>>> rotation = read_BIOMT_records_vec(rf1);
	rf1.close();


	for (int k = 0; k < rotation.size(); k++)
	{

		//Start logic
		vector<vector<double>> new_ATOMs;
		for (int i = 0; i < list_ATOMs.size(); i++)
		{
			new_ATOMs.push_back(vec_multi_matrix(rotation[k], list_ATOMs[i]));
		}

		//Write them to file
		name_f = "E:/New Matrices/NCS vs MTRIX/";
		string ncsORmtrix = "BIOMT";
		string writing_f = name_f + "new_" + ncsORmtrix + to_string(k + 1) + "_" + virus;
		ofstream of(writing_f);
		if (!of.is_open()) {
			cout << "error in opeining file for reading \n";
			return;
		}
		//should be a function
		cout << writing_f << endl;
		writing_ATOMS_only(of, new_ATOMs, lables);
		of.close();
	}

}

//Reads the whole atoms in pdb file and generate new atoms after applying MTRIX
//and Fixed 2nd SYM matrix for VIRUS 1bev
void body9() {

	string name_f = "E:/New Matrices/Archive pdb/", virus = "pdb1bev.pdb";
	ifstream rf(name_f + virus);
	if (!rf.is_open()) {
		cout << "error in opeining file for reading \n";
		return;
	}

	vector<vector<string>> lables;
	vector<vector<double>> list_ATOMs = read_ATOM_records(rf, lables);
	rf.close();

	ifstream rf1(name_f + virus);
	if (!rf1.is_open()) {
		cout << "error in opeining file for reading \n";
		return;
	}
	vector<vector<vector<double>>> MTRIX = read_MTRIX_records_vec(rf1);
	rf1.close();
	//Initialize the fixed 2nd rotation SYM
	vector<vector<double>> SYM;
	for (int l = 0; l < 4; l++)
	{
		vector<double> mSYM(4);
		SYM.push_back(mSYM);
	}
	SYM[0][0] = -1.000000;  SYM[0][1] = 0.000000;  SYM[0][2] = 0.000000;  SYM[0][3] = 0.00000;
	SYM[1][0] = 0.000000;  SYM[1][1] = 1.000000;  SYM[1][2] = 0.000000;  SYM[1][3] = 196.30000;
	SYM[2][0] = 0.000000;  SYM[2][1] = 0.000000;  SYM[2][2] = -1.000000;  SYM[2][3] = 0.00000;
	SYM[3][0] = 0.000000;  SYM[3][1] = 0.000000;  SYM[3][2] = 0.000000;  SYM[3][3] = 1.00000;


	for (int k = 0; k < MTRIX.size(); k++)
	{
		//Multiply SYM * MTRIX or the otherway around?
		vector<vector<double>> rotation = multi_vec(SYM, MTRIX[k]);

		//Start logic
		vector<vector<double>> new_ATOMs;
		for (int i = 0; i < list_ATOMs.size(); i++)
		{
			new_ATOMs.push_back(vec_multi_matrix(rotation, list_ATOMs[i]));
		}

		//Write them to file
		name_f = "E:/New Matrices/NCS vs MTRIX/";
		string ncsORmtrix = "MTRIX_SYM";
		string writing_f = name_f + "new_" + ncsORmtrix + to_string(k + 1) + "_" + virus;
		ofstream of(writing_f);
		if (!of.is_open()) {
			cout << "error in opeining file for reading \n";
			return;
		}
		//should be a function
		cout << writing_f << endl;
		writing_ATOMS_only(of, new_ATOMs, lables);
		of.close();
	}
}


void write_atoms_tofile(string writing_f, vector<vector<double>> new_ATOMs,
	vector<vector<string>> lables,
	vector<vector<double>> additional_info = vector<vector<double>>())
{
	//Write them to file
	ofstream of(writing_f);
	if (!of.is_open()) {
		cout << "error in opeining file for reading \n";
		return;
	}
	//should be a function
	cout << writing_f << endl;
	for (int i = 0; i < additional_info.size(); i++)
	{
		of << "Remark : ";
		for (int j = 0; j < additional_info[i].size(); j++)
		{
			of <<fixed<<setprecision(6)<< additional_info[i][j] << "\t";
		}
		of << endl;
	}

	writing_ATOMS_only(of, new_ATOMs, lables);
	of.close();
}
//Calculating W - common NCS
void body12() {

	string name_f = "E:/New Matrices/Archive pdb/", virus = "pdb1bev.pdb";
	ifstream rf(name_f + virus);
	if (!rf.is_open()) {
		cout << "error in opeining file for reading \n";
		return;
	}

	vector<vector<string>> lables;
	vector<vector<double>> list_ATOMs = read_ATOM_records(rf, lables);
	rf.close();

	ifstream rf1(name_f + virus);
	if (!rf1.is_open()) {
		cout << "error in opeining file for reading \n";
		return;
	}
	vector<vector<vector<double>>> MTRIX = read_MTRIX_records_vec(rf1);
	rf1.close();
	//Initialize the fixed W which is C*A^-1
	vector<vector<double>> W = MTRIX[1];
	vector<vector<double>> W2 = MTRIX[5];
	vector<vector<double>> A = MTRIX[0];

	/*for (int l = 0; l < 4; l++)
	{
		vector<double> mSYM(4);
		W.push_back(mSYM);
	}
	W[0][0] = -0.174206;   W[0][1] = 0.626685;  W[0][2] = 0.759550;  W[0][3] = 71.372876;
	W[1][0] = -0.303707;  W[1][1] = -0.767936;  W[1][2] = 0.563947;  W[1][3] = -2.419104;
	W[2][0] = 0.936705;  W[2][1] = -0.132438;  W[2][2] = 0.324108;  W[2][3] = -52.960934;
	W[3][0] = 0.000000;  W[3][1] = 0.000000;  W[3][2] = 0.000000;  W[3][3] = 1.000000;
*/
	//Start logic
	//This step is useless if A is identity
	vector<vector<double>> A_ATOM;
	for (int i = 0; i < list_ATOMs.size(); i++)
	{
		A_ATOM.push_back(vec_multi_matrix(A, list_ATOMs[i]));
	}

	name_f = "E:/New Matrices/NCS vs MTRIX/";
	string ncsORmtrix = "W_matrix_2_";
	string writing_f = name_f + "new_" + ncsORmtrix + to_string(0) + "_" + virus;
	
	write_atoms_tofile(writing_f, A_ATOM, lables, A);

	for (int i = 0; i < 5; i++)
	{
		
		vector<vector<double>> new_ATOMs;
		for (int i = 0; i < list_ATOMs.size(); i++)
		{
			new_ATOMs.push_back(vec_multi_matrix(W, A_ATOM[i]));
		}
		//For printing - Its real position
		vector<vector<double>> new_ATOMs_print;
		for (int i = 0; i < new_ATOMs.size(); i++)
		{
			new_ATOMs_print.push_back(vec_multi_matrix(W2, new_ATOMs[i]));
		}
		//Multiplying just to write them in the file
		vector<vector<double>> new_Rot =  multi_vec(W, A);
		vector<vector<double>> new_Rot_print = multi_vec(W2, new_Rot);
		for (int f = 0; f < W2.size(); f++)
		{
			new_Rot_print.push_back(W2[f]);
		}
		
		//Writing to file
		name_f = "E:/New Matrices/NCS vs MTRIX/";
		string ncsORmtrix = "W_matrix_2_";
		string writing_f = name_f + "new_" + ncsORmtrix + to_string(i + 1) + "_" + virus;
		write_atoms_tofile(writing_f, new_ATOMs_print, lables, new_Rot_print);
		//Replace the old matrix with the new one
		A_ATOM = new_ATOMs;
		A = new_Rot;
	}	
}

//Calculating W - common NCS  and generating the whole capsid
// multiplying W2 * W * W*W 
void body12_1() {

	string name_f = "E:/New Matrices/Archive pdb/", virus = "pdb1bev.pdb";
	ifstream rf(name_f + virus);
	if (!rf.is_open()) {
		cout << "error in opeining file for reading \n";
		return;
	}

	vector<vector<string>> lables;
	vector<vector<double>> list_ATOMs = read_ATOM_records(rf, lables);
	rf.close();

	ifstream rf1(name_f + virus);
	if (!rf1.is_open()) {
		cout << "error in opeining file for reading \n";
		return;
	}
	vector<vector<vector<double>>> MTRIX = read_MTRIX_records_vec(rf1);
	rf1.close();
	//Initialize the fixed W which is C*A^-1
	vector<vector<double>> W = MTRIX[1];
	vector<vector<double>> W2_g = MTRIX[5];
	vector<vector<double>> A = MTRIX[0];

	/*for (int l = 0; l < 4; l++)
	{
	vector<double> mSYM(4);
	W.push_back(mSYM);
	}
	W[0][0] = -0.174206;   W[0][1] = 0.626685;  W[0][2] = 0.759550;  W[0][3] = 71.372876;
	W[1][0] = -0.303707;  W[1][1] = -0.767936;  W[1][2] = 0.563947;  W[1][3] = -2.419104;
	W[2][0] = 0.936705;  W[2][1] = -0.132438;  W[2][2] = 0.324108;  W[2][3] = -52.960934;
	W[3][0] = 0.000000;  W[3][1] = 0.000000;  W[3][2] = 0.000000;  W[3][3] = 1.000000;
	*/
	
	//This step is useless if A is identity
	vector<vector<double>> A_ATOM;
	for (int i = 0; i < list_ATOMs.size(); i++)
	{
		A_ATOM.push_back(vec_multi_matrix(A, list_ATOMs[i]));
	}

	name_f = "E:/New Matrices/NCS vs MTRIX/";
	string ncsORmtrix = "W_matrix_2_";
	string writing_f = name_f + "new_" + ncsORmtrix + to_string(0) + "_" + virus;

	write_atoms_tofile(writing_f, A_ATOM, lables, A);
	
	//The logic
	vector<vector<double>> W2 = identity_i_vec(4);
	for (int k = 0; k < 6; k++)
	{
		W2 = multi_vec(W2, W2_g);
				
		for (int i = 0; i < 5; i++)
		{

			vector<vector<double>> new_ATOMs;
			for (int i = 0; i < list_ATOMs.size(); i++)
			{
				new_ATOMs.push_back(vec_multi_matrix(W, A_ATOM[i]));
			}
			//For printing - Its real position
			vector<vector<double>> new_ATOMs_print;
			for (int i = 0; i < new_ATOMs.size(); i++)
			{
				new_ATOMs_print.push_back(vec_multi_matrix(W2, new_ATOMs[i]));
			}
			//Multiplying just to write them in the file
			vector<vector<double>> new_Rot = multi_vec(W, A);
			vector<vector<double>> new_Rot_print = multi_vec(W2, new_Rot);
			//Writing to file
			name_f = "E:/New Matrices/NCS vs MTRIX/";
			string ncsORmtrix = "W_matrix_"+ to_string(k + 1) +"_";
			string writing_f = name_f + "new_" + ncsORmtrix + to_string(i + 1) + "_" + virus;
			write_atoms_tofile(writing_f, new_ATOMs_print, lables, new_Rot_print);
			//Replace the old matrix with the new one
			A_ATOM = new_ATOMs;
			A = new_Rot;
		}
	}
}

//Calculating W - common NCS  and generating the whole capsid
//This function is equivlant to the previous one but more effecient
//Multiplying W2 to the star 
void body12_2() {

	string name_f = "E:/New Matrices/Archive pdb/", virus = "pdb1bev.pdb";
	ifstream rf(name_f + virus);
	if (!rf.is_open()) {
		cout << "error in opeining file for reading \n";
		return;
	}

	vector<vector<vector<double>>> SYM = read_smtry_records_vec(rf);
	vector<vector<string>> lables;
	vector<vector<double>> list_ATOMs = read_ATOM_records(rf, lables);
	rf.close();

	ifstream rf1(name_f + virus);
	if (!rf1.is_open()) {
		cout << "error in opeining file for reading \n";
		return;
	}
	vector<vector<vector<double>>> MTRIX = read_MTRIX_records_vec(rf1);
	rf1.close();
	//Initialize the fixed W which is C*A^-1
	vector<vector<double>> W = MTRIX[7];
	vector<vector<double>> W2_g = MTRIX[1];
	vector<vector<double>> A = MTRIX[0];

	/*for (int l = 0; l < 4; l++)
	{
	vector<double> mSYM(4);
	W.push_back(mSYM);
	}
	W[0][0] = -0.174206;   W[0][1] = 0.626685;  W[0][2] = 0.759550;  W[0][3] = 71.372876;
	W[1][0] = -0.303707;  W[1][1] = -0.767936;  W[1][2] = 0.563947;  W[1][3] = -2.419104;
	W[2][0] = 0.936705;  W[2][1] = -0.132438;  W[2][2] = 0.324108;  W[2][3] = -52.960934;
	W[3][0] = 0.000000;  W[3][1] = 0.000000;  W[3][2] = 0.000000;  W[3][3] = 1.000000;
	*/

	//This step is useless if A is identity
	vector<vector<double>> A_ATOM;
	for (int i = 0; i < list_ATOMs.size(); i++)
	{
		A_ATOM.push_back(vec_multi_matrix(A, list_ATOMs[i]));
	}

	name_f = "E:/New Matrices/NCS vs MTRIX/";
	string ncsORmtrix = "W_matrix_0_";
	string writing_f = name_f + "new_" + ncsORmtrix + to_string(0) + "_" + virus;

	write_atoms_tofile(writing_f, A_ATOM, lables, A);

	//Creating the star look
	vector<vector<double>> star_ATOMs;
	vector<vector<string>> star_lables;
	for (int i = 0; i < 5; i++)
	{

		vector<vector<double>> new_ATOMs;
		for (int i = 0; i < list_ATOMs.size(); i++)
		{
			auto g = vec_multi_matrix(W, A_ATOM[i]);
			new_ATOMs.push_back(g);
			star_ATOMs.push_back(g);
		}
		
		//Multiplying just to write them in the file
		vector<vector<double>> new_Rot = multi_vec(W, A);
		vector<vector<double>> new_Rot_print = multi_vec(W, new_Rot);
		//Writing to file
		name_f = "E:/New Matrices/NCS vs MTRIX/";
		string ncsORmtrix = "W_matrix_" + to_string( 0) + "_";
		string writing_f = name_f + "new_" + ncsORmtrix + to_string(i + 1) + "_" + virus;
		write_atoms_tofile(writing_f, new_ATOMs, lables, new_Rot_print);
		//Replace the old matrix with the new one
		A_ATOM = new_ATOMs;
		A = new_Rot;
	}

	//The logic of W2
	vector<vector<double>> W2_i= identity_i_vec(4);
	vector<vector<double>> ring;
	vector<vector<string>> ring_lables;
	vector<vector<double>> star_cpy = star_ATOMs;
	for (int k = 0; k < 5; k++)
	{
		vector<vector<double>>  star_print;

		W2_i = multi_vec(W2_g, W2_i);

		for (int i = 0; i <star_cpy.size(); i++)
		{
			auto g = vec_multi_matrix(W2_g, star_cpy[i]);
			star_print.push_back(g);
			ring.push_back(g);
			ring_lables.push_back(star_lables[i]);
		}
		//print Atoms star Atoms and matrix W2
		//Writing to file
		name_f = "E:/New Matrices/NCS vs MTRIX/";
		string ncsORmtrix = "W_matrix_star_";
		string writing_f = name_f + "new_" + ncsORmtrix + to_string(k + 1) + "_" + virus;
		write_atoms_tofile(writing_f, star_print, star_lables, W2_i);
		//copy 
		star_cpy = star_print;
	}
	//////Edting SYM
	SYM[1] = MTRIX[1];
	////SYM[1][2][2] = -1;
	//Ring is saved
	//Apply SYM now
	vector<vector<double>> ring_SYM;
	for (int i = 0; i < ring.size(); i++)
	{
		ring_SYM.push_back(vec_multi_matrix(SYM[1], ring[i]));
	}
	//print Atoms ring 
	//Writing to file
	name_f = "E:/New Matrices/NCS vs MTRIX/";
	ncsORmtrix = "W_matrix_ring_";
	writing_f = name_f + "new_" + ncsORmtrix + to_string( 1) + "_" + virus;
	write_atoms_tofile(writing_f, ring_SYM, ring_lables);

}

//Creating Star only using W
void  body12_3(int L) {

	string name_f = "E:/New Matrices/Archive pdb/", virus = "pdb1e57.pdb";
	ifstream rf(name_f + virus);
	if (!rf.is_open()) {
		cout << "error in opeining file for reading \n";
		return ;
	}

	//vector<vector<vector<double>>> SYM = read_smtry_records_vec(rf);
	vector<vector<string>> lables;
	vector<vector<double>> list_ATOMs = read_ATOM_records(rf, lables);
	rf.close();

	ifstream rf1(name_f + virus);
	if (!rf1.is_open()) {
		cout << "error in opeining file for reading \n";
		return;
	}
	vector<vector<vector<double>>> MTRIX = read_BIOMT_records_vec(rf1);
	rf1.close();
	//Initialize the fixed W which is C*A^-1
	vector<vector<double>> W = MTRIX[L];
	vector<vector<double>> A = identity_i_vec(4);//MTRIX[0];

	////This step is useless if A is identity
	vector<vector<double>> A_ATOM = list_ATOMs;
	//for (int i = 0; i < list_ATOMs.size(); i++)
	//{
	//	A_ATOM.push_back(vec_multi_matrix(A, list_ATOMs[i]));
	//}

	//name_f = "E:/New Matrices/NCS vs MTRIX/";
	//string ncsORmtrix = "W_matrix_0_";
	//string writing_f = name_f + "new_" + ncsORmtrix + to_string(0) + "_" + virus;

	//write_atoms_tofile(writing_f, A_ATOM, lables, A);

	//Creating the star look
	vector<vector<double>> star_ATOMs;
	vector<vector<string>> star_lables;
	for (int i = 0; i < 5; i++)
	{

		vector<vector<double>> new_ATOMs;
		for (int i = 0; i < list_ATOMs.size(); i++)
		{
			auto g = vec_multi_matrix(W, A_ATOM[i]);
			new_ATOMs.push_back(g);
			star_ATOMs.push_back(g);
		}

		//Multiplying just to write them in the file
		vector<vector<double>> new_Rot = multi_vec(W, A);
		vector<vector<double>> new_Rot_print = multi_vec(W, new_Rot);
		//Writing to file
		name_f = "E:/New Matrices/NCS vs MTRIX/";
		string ncsORmtrix = "W_matrix_" + to_string(L) + "_";
		string writing_f = name_f + "new_" + ncsORmtrix + to_string(i + 1) + "_" + virus;
		write_atoms_tofile(writing_f, new_ATOMs, lables, new_Rot_print);
		//Replace the old matrix with the new one
		A_ATOM = new_ATOMs;
		A = new_Rot;
	}

	//return star_ATOMs;
}

// Creating Star only using ncs from ncs.ofm
//one based
void  body12_4(int L) {

	string name_f = "E:/New Matrices/Archive pdb/", virus = "pdb1bev.pdb";
	ifstream rf(name_f + virus);
	if (!rf.is_open()) {
		cout << "error in opeining file for reading \n";
		return;
	}

	//vector<vector<vector<double>>> SYM = read_smtry_records_vec(rf);
	vector<vector<string>> lables;
	vector<vector<double>> list_ATOMs = read_ATOM_records(rf, lables);
	rf.close();
	
	name_f = "E:\\New Matrices\\Fortran\\out\\", virus = "ncs.ofm";
	ifstream rf1(name_f + virus);
	if (!rf1.is_open()) {
		cout << "error in opeining file for reading \n";
		return;
	}
	vector<vector<double>> MTRIX = read_ROTATION_records(rf1, L);
	rf1.close();
	//Initialize the fixed W which is C*A^-1
	vector<vector<double>> W = MTRIX;
	vector<vector<double>> A = identity_i_vec(4);//MTRIX[0];

												 ////This step is useless if A is identity
	vector<vector<double>> A_ATOM = list_ATOMs;
	//for (int i = 0; i < list_ATOMs.size(); i++)
	//{
	//	A_ATOM.push_back(vec_multi_matrix(A, list_ATOMs[i]));
	//}

	//name_f = "E:/New Matrices/NCS vs MTRIX/";
	//string ncsORmtrix = "W_matrix_0_";
	//string writing_f = name_f + "new_" + ncsORmtrix + to_string(0) + "_" + virus;

	//write_atoms_tofile(writing_f, A_ATOM, lables, A);

	//Creating the star look
	vector<vector<double>> star_ATOMs;
	vector<vector<string>> star_lables;
	for (int i = 0; i < 5; i++)
	{

		vector<vector<double>> new_ATOMs;
		for (int i = 0; i < list_ATOMs.size(); i++)
		{
			auto g = vec_multi_matrix(W, A_ATOM[i]);
			new_ATOMs.push_back(g);
			star_ATOMs.push_back(g);
		}

		//Multiplying just to write them in the file
		vector<vector<double>> new_Rot = multi_vec(W, A);
		vector<vector<double>> new_Rot_print = multi_vec(W, new_Rot);
		//Writing to file
		name_f = "E:/New Matrices/NCS vs MTRIX/";
		string ncsORmtrix = "ncs_matrix_" + to_string(L) + "_";
		string writing_f = name_f + "new_" + ncsORmtrix + to_string(i + 1) + "_.pdb";
		write_atoms_tofile(writing_f, new_ATOMs, lables, new_Rot_print);
		//Replace the old matrix with the new one
		A_ATOM = new_ATOMs;
		A = new_Rot;
	}

	//return star_ATOMs;
}
int main()
{
	//edit_file();
	//generate_files("2btv");
	//body();
	//sort_space_group();
	//count_space_group();
	//body3();
	//body_1();
	//body_11();
	//body5();
	//body6();
	//body7();
	//body8();
	//body9();
	//body10();
	//preprocessing_input();
	//extract_ranges();
	//body();
	//body0();
	//body12_2();
	for (size_t i = 1; i <= 10; i++)
	{
		body12_4(i);
	}
	//system("pause");
	return 0;
}