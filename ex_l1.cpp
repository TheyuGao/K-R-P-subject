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
//问题1：map有序，vector无序，是否影响？
//
//向右拓展，传入参数k_r_p，

string Extend_L (int kmer_length, map<string, vector<pair<int, int> > >k_r_p,map<string,int >kmer_map,int peak_start,int peak_end,vector<int >data_dep,vector<string>data_str,int seed_read,int kmer_avr,int dep_search){
    //找出一条公共k-mer数量最多的read
    int dep=Seq_Dep;
    int cut_num=0;
    int max_read=seed_read;
    int max_read_num=0;
    //中间变量，存储当前合并的两条k-mer与其位置，与现在的contig
    int first_kmer_p;
    int second_kmer_p=0;
    string first_kmer,second_kmer,last_first_kmer,last_second_kmer,com_read,temp_contig="";
    //统计两个k-mer中间碱基个数的辅助数组，并用以固定长度
    //用以统计用过的read，用过后不再选取
    map<string,bool>had_used_kmer;
    //遍历所有map<string, int >类型的指针
    map<string, int >::iterator l_it;
map<string, int >::iterator l_it1;
    //遍历所有vector<int>类型的指针
    vector<int>::iterator it;
    vector<int>::iterator it1;
    //中间变量，记录一条read中的公共k-mer位置
    vector<int>kmer_position;
    vector<bool> had_used_read(data_str.size(),false);
    //帮助决策的两个变量
    int star_position=0;
map<int,vector<pair<int,int > > >candi_read;
int Ch_read=55;
if(dep_search==1)
    Ch_read=5;
cout<<"Ch_read value"<<Ch_read<<endl;


int last_max_read=0;
bool ex_finish=false;




//for(int i=0;i<read_kmer_num.size();i++){
//        if(read_kmer_num[i]>0)
//            read_kmer_num_count++;
//    }
//    cout<<"有数量的read_kmer_num"<<read_kmer_num_count<<endl;
//
////    找出包含公共kmer个数最多的read
//    for(int i=0;i<read_kmer_num.size();i++){
//        if(max_read_num<read_kmer_num[i]) {
//            max_read = i;
//            max_read_num=read_kmer_num[i];
//        }
//    }
    ////这条read标记为已经使用改为数组初始化为0或者map,这里需要更改为

    //统计其所有公共k-mer位置信息
    //先把read拿出来，定义一个map<string,int >::iterator it=data.begin();在定义一个sting read=(it+max_read)->first
    //把for循环包含到while里面，下面找一个新的read，包含最后两个（优先），也可以最后一个进行拓展
  //  for(int i=0;i<=read.length-kmer_length;i++){
        //拿出来所有kmer
        //判断是不是公共kmer(大于多少小于多少)
        //判断是公共kmer:到k_r_p把其read拿出来，存起来
        //存一下前一个公共kmer包含的read位置，当前公共kmer 的read 信息，两个read信息里面找相同read,距离可以确定
        //距离中间定义一个大小数组，ATGC每个位置
    //}

map<string, vector<pair<int, int> > >::iterator k_r_it;

    //现在最多kmer的read已经找到，为it+=read_kmer_num.begin(),it->first


    string contig="";
    star_position=data_str[max_read].length()-kmer_length;
    //结束条件为contig不再更新
    while(max_read!=-1){
        last_max_read=max_read;

        //cout<<data_str[max_read]<<"is max read"<<endl;
        //截取kmer并保存位置
        for (int j = star_position; j >= 0;j--) {
            const string& kmer = data_str[max_read].substr(j, kmer_length);
            if(kmer_map[kmer]>=peak_start&&kmer_map[kmer]<=peak_end){
                kmer_position.push_back(j);
            }
        }
        //cout<<"已经使用kmer个数"<<had_used_kmer.size()<<endl;
        //cout<<"当前read k-mer个数"<<kmer_position.size()<<endl;




    //现在开始寻找不重叠(包括直接相邻的kmer,即两kmer中间0个字符)的kmer，如果重叠继续寻找别的kmer，直到找到两条不重叠的
    for(int ii= 0;ii<kmer_position.size()-1;ii++){
        //记录两条k-mer在该read上的位置
        first_kmer_p=kmer_position[ii];
        second_kmer_p=kmer_position[ii+1];
        //cout<<"总共"<<kmer_position.size()<<"个kmer，当前比对到"<<ii+1<<"个"<<endl;
        //cout<<"first_kmer"<<first_kmer_p<<"second_kmer"<<second_kmer_p<<endl;
        //用指针找出具体read是什么，并截取出两条k-mer

        first_kmer=data_str[max_read].substr(first_kmer_p,kmer_length);
        second_kmer=data_str[max_read].substr(second_kmer_p,kmer_length);
        //cout<<"first_kmer"<<first_kmer<<"second_kmer"<<second_kmer<<endl;
        //如果两个kmer重叠:直接把这部分加入contig中
        if(kmer_position[ii] - kmer_position[ii+1]<=kmer_length) {

            //cout<<"两kmer重叠"<<endl;
            if(ii==data_str[max_read].length()&&contig.length()==0)
                temp_contig = first_kmer+temp_contig;
            vector<int >Chimera;

            for(int i=0;i<k_r_p[first_kmer].size();i++){
                if(Chimera.size()>Ch_read)
                    break;
                int search_read=k_r_p[first_kmer][i].first;
                int search_read_p=k_r_p[first_kmer][i].second;
                for(int j=search_read_p-1;j>=0;j--){

                    string kmer =data_str[search_read].substr(j,kmer_length);
                    if(kmer==second_kmer){
                        Chimera.push_back(search_read);
                        break;
                    }
                }
            }





            if(Chimera.size()>Ch_read || ex_finish==true){
                if(ii==0&&contig.length()==0)
                    temp_contig += first_kmer;

                temp_contig = data_str[max_read].substr(second_kmer_p,first_kmer_p-second_kmer_p)+temp_contig;
            }
            if(Chimera.size()<Ch_read && ex_finish==false){
                second_kmer=last_second_kmer;
                first_kmer=last_first_kmer;
                for(int l=0;l<Chimera.size();l++){
                    had_used_read[Chimera[l]]=true;
                }
                Chimera.clear();
                temp_contig="";
                cout<<"c-b小于Ch_read"<<endl;
                break;
            }
            //cout<<"重叠部分截取"<<temp_contig<<endl;
        }

        //如果找到两个kmer不重叠，先固定长度，去k_r_p中寻找
        else{
//if(contig.length()==1980){
//cout<<"断点read为"<<data_str[max_read];
//cout<<"当前temp_contig为"<<temp_contig<<endl;}
            candi_read.clear();
            //cout<<"kmer不重叠"<<endl;
            //找到同时拥有这两个kmer的read,然后确定其长度，统计所有拥有这两条kmer的read中间的长度，这个可以写一个单独的函数
            time_t start_time_length = time(NULL);
            int count_break=0;
            vector<int >Chimera;
            for(int i=0;i<k_r_p[first_kmer].size();i++){
                for(int j=0;j<k_r_p[second_kmer].size();j++){
                    if(k_r_p[first_kmer][i].first==k_r_p[second_kmer][j].first){
                        //统计这条read在两条k-mer之间的长度
                        if(k_r_p[first_kmer][i].second>k_r_p[second_kmer][j].second){
                            //存储长度信息与read
                            if(count_break==1000)
                                break;
                            count_break++;
                            candi_read[k_r_p[first_kmer][i].second-k_r_p[second_kmer][j].second].push_back(make_pair(k_r_p[first_kmer][i].first,k_r_p[second_kmer][j].second+kmer_length));
                            Chimera.push_back(k_r_p[first_kmer][i].first);
                        }

                    }

                }
                if(count_break==1000)
                    break;
            }
if(contig.length()==1980){
cout<<"断点read为"<<data_str[max_read]<<endl;
cout<<"first_kmer and seocnd_kmer"<<first_kmer<<"  "<<second_kmer<<endl;
cout<<"ch_read"<<Ch_read<<endl;
cout<<"chimera"<<Chimera.size()<<endl;
cout<<Chimera[1]<<endl;
cout<<Chimera[2]<<endl;
cout<<Chimera[3]<<endl;
}
            time_t end_time_length = time(NULL);
            //cout<<"统计长度用时"<<(end_time_length-start_time_length)<<endl;
            //上面循环结束，已经统计出所有拥有这两条k-mer的read的长度情况
            //找出该长度最多的一条read的两k-mer中间部分直接加入contig中
            const int max_value_p=candi_read.begin()->first;
            //cout<<"中间需要决策碱基"<<max_value_p<<endl;

            //首先contig加入第一条k-mer
            if(ii==0&&contig.length()==0){
                    temp_contig=first_kmer+temp_contig;
                   // cout<<"重叠kmer直接加入"<<temp_contig<<endl;
            }



            if(Chimera.size()<Ch_read && ex_finish==false){
                first_kmer=last_first_kmer;
                second_kmer=last_second_kmer;
                temp_contig="";
                for(int t=0;t<Chimera.size();t++){
                    had_used_read[Chimera[t]]=true;
                }
                break;
            }


            //长度符合要求的read一次性统计所有中间碱基
            if(max_value_p<=kmer_length){
                temp_contig=second_kmer+temp_contig;
            }
            else{


            int count_A[max_value_p-kmer_length]= {0};int count_T[max_value_p-kmer_length]= {0};int count_G[max_value_p-kmer_length]= {0};int count_C[max_value_p-kmer_length]= {0};

            time_t start_time_position = time(NULL);
            vector<pair<int,int > >temp_vec=candi_read.begin()->second;
            for(int itx=0;itx<temp_vec.size();itx++){
                int read=temp_vec[itx].first;
                for(int k=0;k<max_value_p-kmer_length;k++){
                 if(data_str[read][temp_vec[itx].second+k]=='A')
                     count_A[k]++;
                 if(data_str[read][temp_vec[itx].second+k]=='T')
                     count_T[k]++;
                 if(data_str[read][temp_vec[itx].second+k]=='G')
                     count_G[k]++;
                 if(data_str[read][temp_vec[itx].second+k]=='C')
                     count_C[k]++;
                }
            }
    for(int i=max_value_p-kmer_length;i>=0;i--){
        if(max({count_A[i],count_T[i],count_C[i],count_G[i]})==count_A[i])
            temp_contig="A"+temp_contig;
        if(max({count_A[i],count_T[i],count_C[i],count_G[i]})==count_T[i])
            temp_contig="T"+temp_contig;
        if(max({count_A[i],count_T[i],count_C[i],count_G[i]})==count_G[i])
            temp_contig="G"+temp_contig;
        if(max({count_A[i],count_T[i],count_C[i],count_G[i]})==count_C[i])
            temp_contig="C"+temp_contig;
    }





            //最后是加入后一条kmer
            temp_contig=second_kmer+temp_contig;
        }
            }
        //
    }
    //cout<<temp_contig<<endl;

//现在first_kmer与second_kmer存放数据为最后两个kmer，现在找到拥有这两个kmer且右侧公共Kmer最多的read,前两个for循环用以寻找拥有这两条k-mer的read
//后两个for循环找出其公共k-mer并排序，循环体需要优化
//cout<<"当前第一个kmer"<<first_kmer<<endl;
//cout<<"当前第二个kmer"<<second_kmer<<endl;

    max_read=-1;
int kmer_num=0;
ex_finish=false;
vector<pair<int,int > >first_kmer_read=k_r_p[first_kmer];
    for(int i=0;i<k_r_p[first_kmer].size();i++){
        if(had_used_read[first_kmer_read[i].first]==true)
            continue;
         int search_read=k_r_p[first_kmer][i].first;
         int search_read_p=k_r_p[first_kmer][i].second;
         int start_p;
         for(int j=search_read_p-1;j>=0;j--){
             string kmer =data_str[search_read].substr(j,kmer_length);
             if(kmer==second_kmer){
                 start_p=j;
                 int kmer_num=0;
                 for(int k=j-1;k>=0;k--){
                     if(had_used_kmer[data_str[first_kmer_read[i].first].substr(k,kmer_length)]==true){
                         kmer_num=0;
                         break;
                     }
                     if(kmer_map[data_str[first_kmer_read[i].first].substr(k,kmer_length)]>=peak_start
                     &&kmer_map[data_str[first_kmer_read[i].first].substr(k,kmer_length)]<=peak_end){
                         star_position=k;
                         kmer_num++;
                     }
                 }
                 if(kmer_num>0){
                     max_read=first_kmer_read[i].first;
                     star_position=start_p;
                     break;
                 }
             }
         }
        if(max_read!=-1)
            break;
    }


    if(max_read==-1&&last_first_kmer==first_kmer&&last_max_read!=-1){
        ex_finish=true;
        max_read=last_max_read;
    }
    last_first_kmer=first_kmer;
    last_second_kmer=second_kmer;





    string finish_str="";

    if (max_read==-1){
        cout<<"拓展完毕"<<endl;
        contig.erase(contig.length()-cut_num,cut_num);
        return contig;
    }
    if(max_read==-1){
        int max_length = 0;
        string max_length_read;
        for(int j = 0;j<k_r_p[second_kmer].size();j++){
            if(data_str[k_r_p[second_kmer][j].first].find(second_kmer)>max_length){
                max_length =data_str[k_r_p[second_kmer][j].first].find(second_kmer);
                max_length_read = data_str[k_r_p[second_kmer][j].first];
            }
        }
        for(int j =0 ;j<max_length_read.find(second_kmer)-1;j++){
            finish_str+=max_length_read[j];
        }
    }


    //candi_read中存储了所有包含最后两个k-mer的read，拿出最多的一个
    if(contig.length()==0)
        cut_num=temp_contig.length();
    contig=temp_contig+contig;
    for(int i=0;i<temp_contig.length();i++){
        if((contig.length()-i)>kmer_length){
            string kmer=contig.substr(i,kmer_length);
            if(kmer_map[kmer]>=peak_start&&kmer_map[kmer]<=peak_end){
                if(kmer_map[kmer]-dep/2>peak_start){
                    kmer_map[kmer]-=dep/2;
                    had_used_kmer[kmer]=true;
                }
                else{
                    had_used_kmer[kmer]=true;
                }
            }
        }
    }
    if(max_read==-1)
        contig=finish_str+contig;
    cout<<"当前contig长度为"<<contig.length()<<endl;
    temp_contig="";
    kmer_position.clear();
    //这层for结束了，这个read的所有kmer中间位置信息也结束了
    }  //接下来挑选下一条read，目前想法为同时拥有上一条read最右侧两个公共kmer且右侧至少还有一个公共kmer，直到遍历结束，向右拓展结束
    contig.erase(contig.length()-cut_num,cut_num);
return contig;



}
string final_contig_Left(map<string,int >kmer_map,int kmer_length,string contig,int kmer_avr,int dep_search) {
    //cout << "拓展前contig为" << contig << endl;
    string temp_contig = contig.substr(0, kmer_length - 1);
    string temp_A = "A" + temp_contig;
    string temp_T = "T" + temp_contig;
    string temp_G = "G" + temp_contig;
    string temp_C = "C" + temp_contig;
    int min_fre = kmer_avr/3;
    if(dep_search== 1)
        min_fre = 0;
    else
        min_fre=90;
    cout<<temp_A<<kmer_map[temp_A]<<endl;
    cout<<temp_T<<kmer_map[temp_T]<<endl;
    cout<<temp_G<<kmer_map[temp_G]<<endl;
    cout<<temp_C<<kmer_map[temp_C]<<endl;
    while (kmer_map[temp_A] > min_fre || kmer_map[temp_T] > min_fre || kmer_map[temp_G] > min_fre ||kmer_map[temp_C] > min_fre) {
        int max_cov = max({kmer_map[temp_A], kmer_map[temp_T], kmer_map[temp_G], kmer_map[temp_C]});

        if (max_cov == kmer_map[temp_A]) {
            contig = "A"+contig;
            temp_contig = "A"+temp_contig;
            kmer_map[temp_A]=0;
            kmer_map[temp_T]=0;
            kmer_map[temp_C]=0;
            kmer_map[temp_G]=0;
        }
        if (max_cov == kmer_map[temp_T]) {
            contig = "T"+contig;
            temp_contig = "T"+temp_contig;
            kmer_map[temp_A]=0;
            kmer_map[temp_T]=0;
            kmer_map[temp_C]=0;
            kmer_map[temp_G]=0;
        }
        if (max_cov == kmer_map[temp_G]) {
            contig = "G"+contig;
            temp_contig = "G"+temp_contig;
            kmer_map[temp_A]=0;
            kmer_map[temp_T]=0;
            kmer_map[temp_C]=0;
            kmer_map[temp_G]=0;
        }
        if (max_cov == kmer_map[temp_C]) {
            contig = "C"+contig;
            temp_contig = "C"+temp_contig;
            kmer_map[temp_A]=0;
            kmer_map[temp_T]=0;
            kmer_map[temp_C]=0;
            kmer_map[temp_G]=0;
        }
        temp_contig = temp_contig.substr(0, kmer_length-1);
        temp_A = "A" + temp_contig;
        temp_T = "T" + temp_contig;
        temp_G = "G" + temp_contig;
        temp_C = "C" + temp_contig;
    }
    //cout << "右侧补充拓展完毕,contig为" << contig << endl;
    return contig;
}
