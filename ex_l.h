#ifndef LI_H
#define LI_H

#include<iostream>
#include<math.h>
#include<stdlib.h>
#include<sstream>
#include <map>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <string.h>
#include "loadreads.h"
#include "GeneralSet.h"
#include"kmer.h"

using namespace std;

//这里map用以存储kmer
string Extend_L (int kmer_length, map<string, vector<pair<int, int> > >k_r_p,map<string,int >kmer_map,int peak_start,int peak_end,vector<int >data_dep,vector<string>data_str,int seed_read,int kmer_avr,int dep_search);
string final_contig_Left(map<string,int >kmer_map,int kmer_length,string contig,int kmer_avr,int dep_search);
#endif
