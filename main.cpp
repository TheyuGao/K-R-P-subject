#include "GeneralSet.h"
#include "loadreads.h"
#include "kmer.h"
#include "splice.h"
#include "ex_r.h"
#include "ex_l.h"
#include <iostream>
#include <ctime>
#include <errno.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <time.h>
#include <getopt.h>
#include <map>

using namespace std;


struct option opts[] = {
        {"kmer_length",   required_argument,   0,   'k'},
        {"out_dir",        required_argument,   0,   'o'},
        {"file_type",          required_argument,         0,   't'},
        {"help",          no_argument,         0,   'h'},
        {"double_stranded_mode", no_argument,  0,   OPT_DOUBLE_STRANDED_MODE},
        {"fr",     required_argument,   0,   OPT_FR_STRAND},
        {"left",          required_argument,   0,   OPT_LEFT},
        {"right",         required_argument,   0,   OPT_RIGHT},
        {"singlefile",    required_argument,   0,   OPT_SINGLEFILE},
        {0,0,0,0}

};


string usage() {

    stringstream usage_info;
    usage_info
            << endl
            << "===============================================================================" << endl
            << " IsoTree Usage " << endl
            << "===============================================================================" << endl
            << " ** Options: **" <<endl
            << "  -k <int>: length of kmer, default 21. " << endl
            << "  -o <string>: output directory. " << endl
            << "  -p : paired-end reads. " << endl
            << "  -t <string>: type of file, fa or fq. " << endl
            << "  -h : help information. " << endl
            << " If pair end reads: " << endl
            << "  --left <string>: left reads file name (.fasta). " << endl
            << "  --right <string>: right reads file name (.fasta). " << endl
            << " If single end reads: " << endl
            << "  -singlefile <string>: reads file name (.fasta). " << endl
            << "  --double_stranded_mode: indicate the pair-end read is double stranded mode" << endl
            << "  --fr <int>: only used for pair-end reads. 1: --1--> <--2--  2: <--1-- --2-->  3: --1--> --2--> or <--1-- <--2--, default 1. " << endl
            << "===============================================================================" << endl
            << endl;

    return usage_info.str();

}


int parse_options(int argc, char* argv[]) {

    int option_index = 0;
    int next_option;
    do {
        next_option = getopt_long(argc, argv, "k:o:t:h", opts, &option_index);
        switch (next_option) {
            case -1:
                break;
            case 'k':
                g_kmer_length = atoi(optarg);
                break;
            case 'o':
                out_dir = optarg;
                break;
            case 'h':
                g_help = true;
                break;
            case 't':
                g_file_type = optarg;
                break;
            case OPT_DOUBLE_STRANDED_MODE:
                g_double_stranded_mode = true;
                break;
            case OPT_FR_STRAND:
                g_fr_strand = atoi(optarg);
                break;
            case OPT_LEFT:
                g_left_file = optarg;
                break;
            case OPT_RIGHT:
                g_right_file = optarg;
                break;
            case OPT_SINGLEFILE:
                g_reads_file = optarg;
                break;
            default:
                exit(1);
        }

    } while (next_option != -1);

    if (g_help) {
        cout << usage();
        exit (1);
    }



    if (g_kmer_length > 32) {
        cout << "Error: the kmer length should shorter than 32." << endl;
        exit(1);
    }

    if (g_reads_file.length() > 1)
        g_is_paired_end = false;

    if (g_fr_strand != 1 && g_fr_strand != 2 && g_fr_strand != 3) {
        cout << "Error: --fr can only be 1, 2 or 3" << endl;
        exit(1);
    }

    return 0;

}



int main(int argc, char* argv[]){
int kmer_avr;
    time_t s_time = time(NULL);
    data.clear();
    int parse_ret = parse_options(argc,argv);
    if (parse_ret)
        return parse_ret;
    if (!g_is_paired_end) {
        load_reads(g_reads_file,data, false,data_dep,data_str,g_fr_strand,kmer_avr);
    } else {
        if (g_double_stranded_mode) {
            load_reads(g_left_file, data, false,data_dep,data_str,g_fr_strand,kmer_avr);
            load_reads(g_right_file, data,false,data_dep,data_str,g_fr_strand,kmer_avr);
        } else {
            if (g_fr_strand == 2) {//--1-->  <--2--
                load_reads(g_left_file, data, false,data_dep,data_str,g_fr_strand,kmer_avr);
                load_reads(g_right_file, data, true,data_dep,data_str,g_fr_strand,kmer_avr);
            }
            if (g_fr_strand == 1) {//<--1-- --2-->
                load_reads(g_left_file, data, true,data_dep,data_str,g_fr_strand,kmer_avr);
                load_reads(g_right_file, data, false,data_dep,data_str,g_fr_strand,kmer_avr);
            }
            if (g_fr_strand == 3) {//--1--> --2-->  or <--1-- <--2--
                load_reads(g_left_file, data, false,data_dep,data_str,g_fr_strand,kmer_avr);
                load_reads(g_right_file, data, false,data_dep,data_str,g_fr_strand,kmer_avr);
            }
        }
    }

    if (g_is_paired_end)
        g_mid_read_id = data.size() / 2;
    else
        g_mid_read_id = data.size();

 

    cout<<kmer_avr<<" kmer_avr"<<endl;
    cout<<data.size()<<" reads have found!"<<endl;
    get_kmer(g_kmer_length,data,data_dep,data_str,dep_search);
   
    cout<<"SeqDep is "<<Seq_Dep<<endl;

//获取种子
    int seed=get_seed(kmer_map,data_str,g_kmer_length,first_p);

//向右拓展

    string contig_r=Extend_R(g_kmer_length,k_r_p,kmer_map,first_p,second_p,data_dep,data_str,seed,kmer_avr,dep_search);

    cout<<"向右拓展完毕，现在开始向左拓展"<<endl;
    cout<<"当前右侧contig为"<<contig_r<<endl;
//向左拓展
   string contig_l=Extend_L(g_kmer_length,k_r_p,kmer_map,first_p,second_p,data_dep,data_str,seed,kmer_avr,dep_search);
    cout<<"当前左侧contig为"<<contig_l<<endl;
string contig_exr=final_contig(kmer_map,g_kmer_length,contig_r,kmer_avr,first_p,dep_search);
    if(dep_search==3){
        cout<<"合并完成，环形基因组contig为"<<contig_exr<<endl;
        time_t a_time = time(NULL);
        cout << "Success! (elapsed time: " << (a_time - s_time) << " s)" << endl;
        return 1;
    }
string Contig_Ex_L = final_contig_Left(kmer_map,g_kmer_length,contig_l,kmer_avr,dep_search);
cout<<"合并完成，contig为"<<Contig_Ex_L+contig_exr<<endl;
    time_t e_time = time(NULL);
    cout << "Success! (elapsed time: " << (e_time - s_time) << " s)" << endl;
    return 1;
}

