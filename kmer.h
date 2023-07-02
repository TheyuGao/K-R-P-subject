#ifndef KMER_H
#define KMER_H

#include "loadreads.h"
#include "GeneralSet.h"
#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <list>
#include <boost/unordered_map.hpp>

using namespace std;
extern std::map<std::string,vector<std::string> >K_R_map;
extern std::map<std::string,vector<std::string> >R_K_map;
//这里map用以存储kmer
bool contains_non_gatc (string kmer);
bool Find_peak(map<string,int > kmer_map,map<string,int >data,int& peak_start,int& peak_end,vector<int >data_seq,vector<string >data_str,int dep,int& dep_search);
void get_kmer(int kmer_length,map<string,int >& data,vector<int >data_seq,vector<string >data_str,int& dep_search);
#endif
