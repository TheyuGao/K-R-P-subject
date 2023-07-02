// This file was modified from a file named common.h in Binpack.
// The original copyright info is Copyright (c) 2013, The Broad Institute, Inc.
// Distributed under the  Distributed under Binpack Software LICENSE.



#include "GeneralSet.h"
#include <cstring>
#include <cstdlib>
#include <errno.h>
#include <map>

using namespace std;

//k-mer
int Seq_Dep=0;
int g_kmer_length = 21;
int g_min_kmer_coverage = 1;//1;
float g_min_kmer_entropy = 0.0f;
int g_min_seed_coverage = 1;
float g_min_seed_entropy = 1.5f;
float g_min_ratio_non_error = 0.04f;
std::map<std::string,int>kmer_map;
std::map<string, vector<pair<int, int> > >k_r_p;
int first_p;
int second_p;
int dep_search;
//reads
map<string,int> data;
vector<string> data_str;
vector<int >data_dep;
vector<string> candi_K;
int g_mid_read_id = 0;
bool g_is_paired_end = true; 
int g_fr_strand = 2;
bool g_double_stranded_mode = false; 
int g_pair_gap_length = 200;
int g_max_pair_gap_length = 500;


//contig
extern std::map<std::string,vector<std::string> >K_R_map;
extern std::map<std::string,vector<std::string> >R_K_map;
//file
string g_reads_file = "";
string g_left_file = "";
string g_right_file = "";
string out_dir = ""; 
string g_file_type = "";
//others
bool g_help = false;


