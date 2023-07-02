#ifndef LOADREADS_H
#define LOADREADS_H


#include <vector>
#include <map>
#include <string>


using namespace std;


typedef vector<unsigned long long> read_int_type;

void load_reads(string file, map<string,int>& input_data, bool rev,vector<int >& data_dep,vector<string >& data_str,int fr,int& kmer_avr);
void load_reads_fa(string file, map<string,int>& input_data, bool rev);
void load_reads_fq(string file, map<string,int>& input_data, bool rev,vector<int >& data_dep,vector<string >& data_str,int fr,int& kmer_avr);

#endif