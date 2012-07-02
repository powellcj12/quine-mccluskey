/*
 *  QM.cpp
 *
 *  Computes cononical expressions for functions of minterms and don't cares
 *  Able to produce in both SOP and POS format
 *  Limited to number of literals assigned to MAX_LITERALS (10 default)
 *
 *  To run the program, first compile on a UNIX based machine (or Windows with a UNIX compiler):
 *       g++ QM.cpp -o QM
 *       ./QM
 *
 *  If using an input.txt file, there should NOT be a new line at the bottom.
 *
 *  Created by Charlie Powell on 4/15/11.
 *  Copyright 2011 Tufts University. All rights reserved.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <bitset>
#include <cmath>
#include <algorithm>

using namespace std;

const int MAX_LITERALS = 10; //Maximum number of literals allowed for expressions
const int MAX_TERMS = 1024; //Maximum number of terms (equals 2^MAX_LITERALS)

struct cube
{
	vector< bitset<MAX_LITERALS> > bits; //All the minterms associated with the cube
	bitset<MAX_LITERALS> d; //Indicates which literals are not in the cube
	bool used; //Marks whether or not a cube was used in PI table
};

void setupQM(string input);
vector<cube*> runQM(vector<int> &minterms, vector<int> &dontcares, vector<int> &terms, int literals);
void takeInput(vector<int> &v1, vector<int> &v2, string in);
void takeInput2(vector<int> &v, int &i, string s);
vector< bitset<MAX_LITERALS> > mergeTerms(vector< bitset<MAX_LITERALS> > &v1, vector< bitset<MAX_LITERALS> > &v2);
vector< vector<cube*> > createCubes(vector<int> &v, int literals);
void combineCubes(vector< vector<cube*> > &v, vector<cube*> &p);
void makeTable(vector< vector<bool> > &t, vector<cube*> &p, vector<int> &m);
int reduceTable(vector< vector<bool> > &t, vector<cube*> &p, vector<cube*> &e);
int dominatedRows(vector< vector<bool> > &t, vector<cube*> &p);
void petrick(vector< vector<bool> > &table, vector<cube*> &primes, vector<cube*> &essentials, int literals);
void convertSOPCube(cube* c, int l);
void convertPOSCube(cube* c, int l);

int main()
{
	int input;
	cout << "Please enter '1' to take input from a file (must be called \"input.txt\", ";
	cout << "located in the same directory as this file, and must NOT have an empty line at the end)";
	cout << ", or '2' to manually type input: ";
	cin >> input;
	cout << endl;
	
	if(input == 1)
	{
		ifstream fin;
		fin.open("input.txt");
		string tempIn;
		
		while(!fin.eof())
		{
			fin >> tempIn;
			setupQM(tempIn);
		}
	}
	else
	{
		cout << "Please input in the form 'm(#,#..)+d(#,#..)'. Don't Cares may be ommitted:" << endl;
		string tempIn;
		cin >> tempIn;
		setupQM(tempIn);
	}
	
	return 0;
}

/*
 * Input - string of the form: m(#,#...)+d(#,#...)
 * Computes terms and dontcares from the input to prepare for main computation
 */
void setupQM(string input)
{
	vector<int> onSet;
	vector<int> offSet;
	vector<int> dontcares;
	
	takeInput(onSet, dontcares, input);
	
	//Merges and sorts vector ot terms and don't cares in order to create zero cubes
	vector<int> onTerms(onSet.size() + dontcares.size());
	merge(onSet.begin(), onSet.end(), dontcares.begin(), dontcares.end(), onTerms.begin());
	sort(onTerms.begin(), onTerms.end());
	
	int literals = (int)ceil(log(onTerms.back() + 1) / log(2));
	
	//Computes the offSet in order to determine cononical POS expression
	//It checks for every minterm and adds it to the offSet if it is not in the onSet or dontcares
	for(int i = 0; i < (int)pow(2.0,literals); i++)
	{
		bool found = false;
		
		for(int j = 0; j < (int)onSet.size(); j++)
		{
			if(onSet[j] == i)
				found = true;
		}
		for(int j = 0; j < (int)dontcares.size(); j++)
		{
			if(dontcares[j] == i)
				found = true;
		}
		
		if(!found)
			offSet.push_back(i);
	}
	
	//Similarly merge and sort vector of terms for offSet zero cubes
	vector<int> offTerms(offSet.size() + dontcares.size());
	merge(offSet.begin(), offSet.end(), dontcares.begin(), dontcares.end(), offTerms.begin());
	sort(offTerms.begin(), offTerms.end());
	
	//Only run the program if not all the terms were given
	if(offSet.size() != 0)
	{
		vector<cube*> SOP_essentials = runQM(onSet, dontcares, onTerms, literals);
		vector<cube*> POS_essentials = runQM(offSet, dontcares, offTerms, literals);
		
		//Printout SOP expression
		cout << input << endl << "= ";
		for(int i = 0; i < (int)SOP_essentials.size(); i++)
		{
			convertSOPCube(SOP_essentials[i],literals);
			
			if(i != (int)SOP_essentials.size() - 1)
				cout << " + ";
		}
		cout << endl << "= ";
	
		//Printout POS expression
		for(int i = 0; i < (int)POS_essentials.size(); i++)
		{
			cout << "(";
			convertPOSCube(POS_essentials[i], literals);
			cout << ")";
		}
		cout << endl << endl;
	}
	else
		cout << input << endl << "= 1" << endl << "= 1" << endl << endl;
}

/*
 * Input - vectors of minterms, dontcares, and terms (combined vectors and minterms)
 *       - number of literals
 *
 * Output - vector containing the essential prime implicants
 */
vector<cube*> runQM(vector<int> &minterms, vector<int> &dontcares, vector<int> &terms, int literals)
{
	vector<cube*> primes;
	vector< vector<bool> > table;
	vector<cube*> essentials;
	
	vector< vector<cube*> > zeroCubes = createCubes(terms, literals);
	combineCubes(zeroCubes, primes);
	makeTable(table,primes,minterms);
	
	//Repeat until Petrick's Method is done or until the table is fully reduced
	bool pDone = false;
	while(!pDone && reduceTable(table,primes,essentials))
	{
		//Reduction of cyclic tables
		if(!dominatedRows(table, primes))
		{
			//Determine if there is a cube larger than all the others left
			int maxIndex = 0;
			bool maxCube = true;
			for(int i = 1; i < (int)primes.size(); i++)
			{
				if(primes[i] -> bits.size() > primes[maxIndex] -> bits.size())
					maxIndex = i;
			}
			for(int i = 0; i < (int)primes.size(); i++)
			{
				if(maxIndex != i && primes[maxIndex] -> bits.size() == primes[i] -> bits.size())
					maxCube = false;
			}
			
			//If so, make it essential and continue to reduce the table
			if(maxCube)
			{
				essentials.push_back(primes[maxIndex]);
				
				for(int j = 0; j < (int)table[maxIndex].size(); j++)
				{
					if(table[maxIndex][j])
					{
						for(int i = 0; i < (int)table.size(); i++)
							table[i].erase(table[i].begin() + j);
						j--;
					}
				}
				
				table.erase(table.begin()+maxIndex);
				primes.erase(primes.begin()+maxIndex);
			}
			//If not, use Petrick's method to determine which combination of remaining terms
			//will have the lowest cost.
			else
			{
				pDone = true;
				petrick(table, primes, essentials, literals);
			}
		}
	}
	
	return essentials;
}

/*
 * Input - vectors v1 and v2 to fill from string in
 * Fills both vectors with use of the helper function
 */
void takeInput(vector<int> &v1, vector<int> &v2, string in)
{
	int i = 2;
	takeInput2(v1, i, in);
	
	if(i != (int)in.length() - 1)
	{
		i += 3;
		takeInput2(v2, i, in);
	}
	
	//Sort vectors for convenience
	sort(v1.begin(), v1.end());
	sort(v2.begin(), v2.end());
}

//Input helper function to fill a single vector, v
void takeInput2(vector<int> &v, int &i, string s)
{
	while(s[i] != '+' && i < (int)s.length())
	{
		int temp = i;
		
		while(s[i] != ',' && s[i] != ')')
			i++;
		
		char* temp2 = new char[i - temp];
		
		for(int j = 0; j < (i - temp); j++)
			temp2[j] = s[temp + j];

		v.push_back(atoi(temp2));
		i++;
	}
}

/*
 * Input - 2 bitset vectors, v1 and v2, to be merged
 *
 * Output - vector containing the combination of vectors v1 and v2 (non-destructive) in sorted order
 */
vector< bitset<MAX_LITERALS> > mergeTerms(vector< bitset<MAX_LITERALS> > &v1, vector< bitset<MAX_LITERALS> > &v2)
{
	vector< bitset<MAX_LITERALS> > ans;
	vector< bitset<MAX_LITERALS> >::iterator i1 = v1.begin(), i2 = v2.begin();
	
	//Repeat until one of the iterators reaches the end
	while(i1 < v1.end() && i2 < v2.end())
	{
		//Push the next least item onto the output vector
		if(i1 -> to_ulong() < i2 -> to_ulong())
		{
			ans.push_back(*i1);
			i1++;
		}
		else
		{
			ans.push_back(*i2);
			i2++;
		}
	}
	
	//Put all remaining elements of the non-empty vector into the output vector
	if(i1 == v1.end())
	{
		while(i2 < v2.end())
		{
			ans.push_back(*i2);
			i2++;
		}
	}
	else
	{
		while(i1 < v1.end())
		{
			ans.push_back(*i1);
			i1++;
		}
	}
	
	return ans;
}

/*
 * Input - vector of terms, v, from which to create 0-cubes with a maximum number of literals
 *
 * Output - vectors containing the 0-cubes sorted by, and in order of, number of 1's
 */
vector< vector<cube*> > createCubes(vector<int> &v, int literals)
{
	vector< bitset<MAX_LITERALS> > terms;
	vector< vector<cube*> > ans;
	
	//Initialize cube vector
	for(int i = 0; i < (int)v.size(); i++)
		terms.push_back(v[i]);
	
	//Add cubes sorted by, and in order of, number of 1's in each cube, starting with 0
	for(int i = 0; i <= literals; i++)
	{
		vector<cube*> temp;
		
		for(int j = 0; j < (int)terms.size(); j++)
		{
			//Add a new cube if it has the correct number of 1's
			if((int)terms[j].count() == i)
			{
				vector< bitset<MAX_LITERALS> > temp2 (1, terms[j]);
				bitset<MAX_LITERALS> temp3 (0);
				
				cube* c = new cube;
				c -> bits = temp2;
				c -> d = temp3;
				c -> used = false;
				temp.push_back(c);
			}
		}
		
		//Add vector whose cubes have i 1's each in term
		ans.push_back(temp);
	}
	
	return ans;
}

/*
 * Input - vector of cubes v, which may be reduced
 *       - vector of prime implicants, p
 *
 * Creates cubes one size large (if possible) and puts all unused cubes onto p
 */
void combineCubes(vector< vector<cube*> > &v, vector<cube*> &p)
{
	vector< vector<cube*> > nextV;

	for(int i = 0; i < (int)v.size() - 1; i++)
	{
		vector<cube*> v1 = v[i], v2 = v[i + 1];
		vector<cube*> temp;
		
		for(int j = 0; j < (int)v1.size(); j++)
		{
			cube* c1 = v1[j];
			vector< bitset<MAX_LITERALS> > b1 = c1 -> bits;
			bitset<MAX_LITERALS> test1 = b1[0];
			unsigned long temp1 = test1.to_ulong();
			
			for(int k = 0; k < (int)v2.size(); k++)
			{
				cube* c2 = v2[k];
				vector< bitset<MAX_LITERALS> > b2 = c2 -> bits;
				bitset<MAX_LITERALS> test2 = b2[0];
				unsigned long temp2 = test2.to_ulong();
				
				/* Test if two cubes can be combined into a large one:
					1. temp2 > temp1 (Lower cube is greater)
					2. c1 -> d == c2 -> d (Same don't cares)
					3. (test2^test1).count() == 1 (Differ by power of 2)
				*/
				if(temp2 > temp1 && (c1 -> d == c2 -> d) && (test1^test2).count() == 1)
				{
					cube* c = new cube;
					c -> bits = mergeTerms(b1, b2);

					//Determine if new cube already exists
					bool found = false;
					for(int l = 0; l < (int)temp.size(); l++)
					{
						if(c -> bits == temp[l] -> bits)
							found = true;
					}

					//If not, add it to the new vector
					if(!found)
					{
						c -> d = c1 -> d | (test1 ^ test2);
						c -> used = false;
						temp.push_back(c);
					}
					//Otherwise, free its memory
					else
						delete c;

					c1 -> used = true;
					c2 -> used = true;
				}
			}
		}

		//Check for prime implicants in the first vector analyzed since it won't be analyzed again
		for(int j = 0; j < (int)v1.size(); j++)
		{
			if(!(v1[j] -> used))
				p.push_back(v1[j]);
		}

		nextV.push_back(temp);
	}
	
	//Check for prime implicants in the last vector since it cannot be compared to anything else
	for(int i = 0; i < (int)v[v.size() - 1].size(); i++)
	{
		if(!(v[v.size() - 1][i] -> used))
			p.push_back(v[v.size() - 1][i]);
	}

	//Attempt to recombine cubes if there are at least 2 sets left to compare
	if(v.size() - 1 > 0)
		combineCubes(nextV, p);
	//Otherwise add them all as prime implicants since they can't be combined with anything
	else
	{
		for(int i = 0; i < (int)v[0].size(); i++)
			p.push_back(v[0][i]);
	}
}

/*
 * Input - vectors for the tablt t, prime implicants p, and minterms m
 *
 * Constructs the PI table using only the minterms specified (omitting dont cares)
 * and the prime implicants computed
 */
void makeTable(vector< vector<bool> > &t, vector<cube*> &p, vector<int> &m)
{
	for(int i = 0; i < (int)p.size(); i++)
	{
		//Row vector for a given cube
		vector<bool> temp;
		
		//Determine which minterms a certain cube covers
		for(int j = 0; j < (int)m.size(); j++)
		{
			bool temp2 = false;
			
			for(int k = 0; k < (int)p[i] -> bits.size(); k++)
			{
				if(p[i] -> bits[k] == m[j])
					temp2 = true;
			}
			//Add value to the row vector
			temp.push_back(temp2);
		}
		//Add row to the table
		t.push_back(temp);
	}
}

/*
 * Input - table t, prime implicants p, and essential prime implicants e
 *
 * Output - the remaining size of the table
 *
 * This reduces the table when a minterm is covered by exactly 1 cube. That cube is
 * added as an essential prime implicant. This is called recursively until the table is
 * empty, or no more essential prime implicants can be determined.
 */
int reduceTable(vector< vector<bool> > &t, vector<cube*> &p, vector<cube*> &e)
{
	//Determine rows which are the only ones covering at least 1 minterm
	vector<int> tempEIndex;
	for(int j = 0; j < (int)t[0].size(); j++)
	{
		//Count the number of minterms covered at the index of the last one covered
		int count = 0, eIndex = -1;
		
		for(int i = 0; i < (int)t.size(); i++)
		{
			if(t[i][j])
			{
				count++;
				eIndex = i;
			}
		}
		
		//Keeps track of new essential prime implicants without duplicating results
		if(count == 1 && !p[eIndex] -> used)
		{
			e.push_back(p[eIndex]);
			tempEIndex.push_back(eIndex);
			p[eIndex] -> used = true;
		}
	}
	
	//Erase each column covered by each of the new prime implicants
	for(int k = 0; k < (int)tempEIndex.size(); k++)
	{
		for(int j = 0; j < (int)t[tempEIndex[k]].size(); j++)
		{
			if(t[tempEIndex[k]][j])
			{
				for(int i = 0; i < (int)t.size(); i++)
					t[i].erase(t[i].begin() + j);
				j--;
			}
		}
	}
	
	sort(tempEIndex.begin(), tempEIndex.end());
	
	//Erase the rows from the table corresponding to the new essential prime implicants
	for(int i = 0; i < (int)tempEIndex.size(); i++)
	{
		t.erase(t.begin() + tempEIndex[i]);
		p.erase(p.begin() + tempEIndex[i]);
		
		//Adjust index values of ramining essential prime implicants to be removed from table
		for(int j = i + 1; j < (int)tempEIndex.size(); j++)
			tempEIndex[j]--;
	}
	
	if(t.size() == 0 || t[0].size() == 0)
		return 0; //Check if the table is empty
	else if(tempEIndex.size() == 0)
		return t[0].size(); //Return size of table if no essential PIs could be determined
	else
		return reduceTable(t,p,e); //Attempt to reduce table now with fewer rows and columns
}

/*
 * Input - table t, prime implicants p, and essential prime implicants e
 *
 * Output - the number of dominated rows that were removed
 *
 * Determines how many rows are dominated by another in the table, and removes those
 * rows from the table.
 */
int dominatedRows(vector< vector<bool> > &t, vector<cube*> &p)
{
	bitset<MAX_TERMS> temp, temp2;
	int count = 0;
	
	for(int i = 0; i < (int)t.size(); i++)
	{
		//Determines if rowCheck is dominated by another row in the table
		vector<bool> rowCheck = t[i];
		bool dominated = false;
		
		for(int j = 0; j < (int)t.size(); j++)
		{
			if(i != j)
			{
				bool tempDom = true;
				temp.reset();
				temp2.reset();

				for(int k = 0; k < (int)rowCheck.size(); k++)
				{
					if(rowCheck[k])
						temp.set(k);
					if(t[j][k])
						temp2.set(k);
					//If row being checked has a 1 but the other one does not, it cannot be dominated
					if(rowCheck[k] && !t[j][k])
						tempDom = false;
				}
				
				//Also check the rows are not just weakly dominated
				if(tempDom && temp != temp2)
					dominated = true;
			}
		}
		
		if(dominated)
		{
			count++;
			t.erase(t.begin()+i);
			p.erase(p.begin()+i);
			i--;
		}
	}
	
	return count;
}

//Performs Petrick's method on the PI table
void petrick(vector< vector<bool> > &table, vector<cube*> &primes, vector<cube*> &essentials, int literals)
{
	//Sets up the vector "final" for distribution
	//It assigns a number, starting from 1, to each prime implicant left
	vector< vector< vector<int> > > final;
	for(int j = 0; j < (int)table[0].size(); j++)
	{
		vector< vector<int> > temp;
		for(int i = 0; i < (int)table.size(); i++)
		{
			if(table[i][j])
			{
				vector<int> temp2(1,i+1);
				temp.push_back(temp2);
			}
		}
		
		final.push_back(temp);
	}
	
	//Distribute terms until there is a single sequence to compare
	while(final.size() > 1)
	{
		//Focus on first 2 sets
		vector< vector<int> > temp1 = final[0], temp2 = final[1];
		final.erase(final.begin(),final.begin()+2);
		vector< vector<int> > newTemp;
		
		//Distribute the sets
		for(int i = 0; i < (int)temp1.size(); i++)
		{
			vector<int> smallTemp1 = temp1[i];
			
			for(int j = 0; j < (int)temp2.size(); j++)
			{
				vector<int> smallTemp2 = temp2[j];
				vector<int> newSmallTemp(smallTemp1.size() + smallTemp2.size());
				set_union(smallTemp1.begin(), smallTemp1.end(), smallTemp2.begin(), smallTemp2.end(), newSmallTemp.begin());
				
				//Removes extraneous 0's from the union in the case of overlap
				int lastIndex = -1;
				for(int k = 0; k < (int)newSmallTemp.size(); k++)
				{
					if(lastIndex == -1 && newSmallTemp[k] == 0)
						lastIndex = k;
				}
				
				if(lastIndex != -1 && lastIndex < (int)newSmallTemp.size())
					newSmallTemp.resize(lastIndex);
				
				newTemp.push_back(newSmallTemp);
			}
		}
		
		//Reduce the distribution by applying the absorption property
		//This will occur if the intersection of 2 sets is equal to one of them
		for(int i = 0; i < (int)newTemp.size(); i++)
		{
			for(int j = i+1; j < (int)newTemp.size(); j++)
			{
				vector<int> intTest(min(newTemp[i].size(),newTemp[j].size()));
				set_intersection(newTemp[i].begin(),newTemp[i].end(),newTemp[j].begin(),newTemp[j].end(),intTest.begin());
				
				if(intTest == newTemp[i])
				{
					newTemp.erase(newTemp.begin()+j);
					j--;
				}
				else if(intTest == newTemp[j])
				{
					newTemp.erase(newTemp.begin()+i);
					j--;
				}
			}
		}
		
		//Push the new set back onto the vector to be distributed later (if necessary)
		final.push_back(newTemp);
	}
	
	//Determine which set costs the least to implement in practice
	//minIndex determines which term it is
	//minValue determines how much it costs (used in determining the minimum cost)
	int minIndex = 0, minValue = INT_MAX;
	for(int i = 0; i < (int)final[0].size(); i++)
	{
		//Cost of a cube the number of literals minus its size (number of don't cares)
		int temp = (int)final[0][i].size();
		for(int j = 0; j < (int)final[0][i].size(); j++)
			temp += literals - (primes[final[0][i][j]-1] -> d.count());
		
		if(temp < minValue)
		{
			minValue = temp;
			minIndex = i;
		}
	}
	
	//Add terms from least costing implementation as essential prime implicants
	vector<int> finalTerms = final[0][minIndex];
	for(int i = 0; i < (int)finalTerms.size(); i++)
		essentials.push_back(primes[finalTerms[i] - 1]);
}

//Output a cube in SOP format for l literals
void convertSOPCube(cube* c, int l)
{
	//Determine which literals to use
	char chars[] = {'A','B','C','D','E','F','G','H','I','J'};
	vector<char> terms;
	for(int i = l - 1; i >= 0; i--)
		terms.push_back(chars[i]);

	for(int i = l - 1; i >= 0; i--)
	{
		//Only print when the literal is not a don't care
		if(c -> d[i] == 0)
		{
			cout << terms[i];
			if(c -> bits[0][i] == 0)
				cout << "'";
		}
	}
}

//Output a cube in POS format for l literals
void convertPOSCube(cube* c, int l)
{
	vector<char> terms;
	
	//Determine which literals to use
	for(int i = l - 1; i >= 0; i--)
		terms.push_back('A' + i);
	
	int count = 0;
	for(int i = l - 1; i >= 0; i--)
	{
		//Only print when the literal is not a don't care
		if(c -> d[i] == 0)
		{
			count++;
			cout << terms[i];
			if(c -> bits[0][i] == 1)
				cout << "'";
			
			if(count < l - (int)c -> d.count())
				cout << "+";
		}
	}
}
