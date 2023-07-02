#include "loadreads.h"
#include "kmer.h"
#include "GeneralSet.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include<algorithm>
using namespace std;

const int MAX_STR = 1024;


string revcomp (const string kmer) {

    string revstring;

    for(int i = kmer.size() -1; i >= 0;i--) {

        char c = kmer[i];
        char revchar;

        switch (c) {

            case 'g':
                revchar = 'c';
                break;

            case 'G':
                revchar = 'C';
                break;

            case 'a':
                revchar = 't';
                break;

            case 'A':
                revchar = 'T';
                break;

            case 't':
                revchar = 'a';
                break;

            case 'T':
                revchar = 'A';
                break;

            case 'c':
                revchar = 'g';
                break;

            case 'C':
                revchar = 'G';
                break;

            default:
                revchar = 'N';
        }

        revstring += revchar;

    }

    return (revstring);

}

void load_reads(string file, map<string,int >& input_data, bool rev,vector<int >& data_dep,vector<string >& data_str,int fr,int& kmer_avr) {

	if (g_file_type == "fa")
		load_reads_fa(file, input_data, rev);
	else
		load_reads_fq(file, input_data, rev,data_dep,data_str,fr,kmer_avr);
}



void load_reads_fa(string file, map<string,int >& input_data, bool rev) {

    time_t s_time = time(NULL);

    fstream in;
    in.open(file.c_str(), fstream::in);

    if (!in.is_open()) {

        cout << "Error! Can't open file " << file << endl;
        exit(1);

    }

    cout << "Loading reads from file " << file << " ..." << endl;

    char temp[MAX_STR];
    string read;
    in.getline(temp, MAX_STR);

    while(getline(in,read))
    {
        if(read[0]=='+')
        {
            getline(in,read);
            if(!contains_non_gatc(read)) { continue; }
            //if(iter==input_data.end()) { input_data.push_back(read); }
            if (rev)
                read = revcomp(read);
            if(input_data.find(read)==input_data.end()) {
                input_data.insert(make_pair(read, 1));
            }
            else
                input_data[read]=input_data[read]+1;
        }

    }


    cout <<"Total load " << input_data.size() << " reads!" << endl;
    in.close();

    time_t e_time = time(NULL);
    cout << "Success! (total cost time: " << (e_time - s_time) << " s)" << endl;

}


void load_reads_fq(string file, map<string,int >& input_data, bool rev,vector<int >& data_dep,vector<string >& data_str,int fr,int& kmer_avr) {

	time_t s_time = time(NULL);
	
	fstream in;
	in.open(file.c_str(), fstream::in);

	if (!in.is_open()) {

		cout << "Error! Can't open file " << file << endl;
		exit(1);

	}

	cout << "Loading reads from file " << file << " ..." << endl;

	char temp[MAX_STR];
	string read;
	in.getline(temp, MAX_STR);

    while(getline(in,read))
{
    if(read[0]=='@')
    {
        getline(in,read);
        if(!contains_non_gatc(read)) { continue; }
        //if(iter==input_data.end()) { input_data.push_back(read); }
        if (rev)
            read = revcomp(read);
        if(input_data.find(read)==input_data.end()) {
            input_data.insert(make_pair(read, 1));
        }


        else
            input_data[read]=input_data[read]+1;
    }

}
    int k=0;
    int kmer_dep;
    map<string,int >::iterator it;
    if(fr==1&&rev== false){
        for(it =input_data.begin();it!=input_data.end();it++){
            k+=it->second;
            data_str.push_back(it->first);
            data_dep.push_back(it->second);
            kmer_dep+=it->first.size()-25;
        }
    }
    if(fr==2&&rev== true){
        for(it =input_data.begin();it!=input_data.end();it++){
            k+=it->second;
            data_str.push_back(it->first);
            data_dep.push_back(it->second);
            kmer_dep+=it->first.size()-25;
        }
    }
    //kmer_dep=kmer_dep/input_data.size();
    //kmer_avr=kmer_dep;
	cout <<"Total load " << k<< " reads!" << endl;
    cout<<"Average kmer dep "<<kmer_dep<<endl;
    cout<<"不重复read个数"<<input_data.size()<<endl;
	in.close();
    cout<<"data_str"<<data_str.size()<<endl;
    cout<<"data_dep"<<data_dep.size()<<endl;
    if(fr==1&&rev== false){
        cout<<data_str[100]<<"   "<<data_dep[100]<<endl;
        cout<<input_data[data_str[100]]<<endl;
    }
	time_t e_time = time(NULL);
	cout << "Success! (total cost time: " << (e_time - s_time) << " s)" << endl;

}


