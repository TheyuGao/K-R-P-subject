
#include "loadreads.h"
#include "GeneralSet.h"
#include"kmer.h"
#include<iostream>
#include<math.h>
#include<stdlib.h>
#include<sstream>
#include <map>
#include <fstream>
#include <algorithm>
#include <cmath>
using namespace std;
//此文件完成了一下工作：
//1.生成了kmer
//2.生成kmer频次信息
//3.初始化read-kmer映射信息
//4.提取出公共kmer并建立com_kmer-read映射信息（com_kmer为公共kmer）


//这里map用以存储kmer
//判断是否具有非ATGC碱基
bool contains_non_gatc (string kmer) {
//从头到尾遍历字符串，判断是否具有非ATGC的碱基
    for(unsigned int i =0; i < kmer.size(); ++i) {

        unsigned char c = kmer[i];

        if(c!='A'&&c!='T'&&c!='G'&&c!='C'&&c!='a'&&c!='t'&&c!='g'&&c!='c'&&c!='N') {
            // cout<<c<<endl;
            return (false);

        }
    }

    return(true);
}
//寻找峰值
bool Find_peak(map<string,int > kmer_map,map<string,int >data,int& peak_start,int& peak_end,vector<int >data_seq,vector<string >data_str,int dep,int& dep_search){
    //先要建立频次数据
    //找峰值范围为测序深度的四分之一到测序深度的四分之三
    //Seq_Dep=reads个数x(read长度-k+1)/实际上产生的kmer个数
    //1.首先算出来最大的Kmer频次 数组大小与最大频次相关 窗口和步长都更改：窗口最大频次的1/100，步长1/1000，min_kmerfre要更改定义，最小值为平均值Seq_Dep,max_kmerfre改为最大频率.
    //

    int max_kmerfre=0;
    int dep_judge=0;

    for(map<string,int >::iterator it=kmer_map.begin();it!=kmer_map.end();it++){
        if(max_kmerfre<it->second) { max_kmerfre = it->second; }
    }
    cout<<max_kmerfre<<"is max_kmer_fre"<<endl;


    //求出平均read长度
    int count_kmer[max_kmerfre]= {0};

    for (map<string,int >::iterator it= kmer_map.begin();it!=kmer_map.end();it++){
        //范围内寻找
        //符合范围，直接在统计的数组中加一 it->second为频率
        count_kmer[it->second]++;
    }
    int U=0;
    Seq_Dep=0;
    for (int i=0;i<max_kmerfre;i++){
        if(count_kmer[i]>0&&i>10) {
            Seq_Dep+=count_kmer[i];
            U++;
        }
    }
    Seq_Dep=Seq_Dep/U;
    dep=Seq_Dep;
    cout<<"Seq dep is"<<Seq_Dep<<endl;

    //现在有了频率数据,生成每个区间的数据
    //找到最大peak和平均peak厚度
    int peak[max_kmerfre]= {0};
    int count_peak=0;
    int step_length=max_kmerfre/1000;
    if (max_kmerfre<1000){
        step_length=1;
    }
    int win_length=step_length*10;
    //非常大的值123444468 100
    //count_kmer[]

    //找最后一个非零区间，确定频次图最右侧的位置,用mark标记
    int mark=0;
    for(int i=0;i<max_kmerfre-win_length;i+=(step_length)){
        int none_zero=0;
        for(int j=0;j<(win_length);j++){

            if(count_kmer[i+j]!=0)
            { none_zero++;
                peak[i]+=count_kmer[i+j];
            }
        }
        if(none_zero>0){
            peak[i]=peak[i]/none_zero;
            mark=i;
        }
        count_peak++;
    }


    cout<<"tttest!!!!!!!!!!!!! "<<endl;

    //有了每个区间的kmer总数/非0个数，开始寻找峰值，条件为比上一个1000和下一个1000都大
    //是峰值则大于Seq_Dep,最后peak[i]最后一个值往前找，如果大于Seq_Dep则认为可能是一个峰，
    map<int,int>count_map;
    //1000和100需要更改为1/100，和1/1000
    // for(int i =0;i<max_kmerfre;i++){
    //     if(peak[i]>0)
    //         cout<<"当前Peak"<<i<<"   当前peak平均值"<<peak[i]<<endl;
    // }

    //cout<<"最后一个非0区间"<<mark<<"    每次左移长度"<<step_length<<endl;

     for(int i=mark;i>=0;i-=step_length){
         dep_judge=0;
         int j=i-step_length;

             while(j>=0&&peak[j]>=peak[j+step_length]){
                 j-=step_length;
             }
             if(j>0){
                 j-=step_length;
                 while(j>=0&&peak[j]>=peak[j+step_length]){
                     j-=step_length;
                 }
                 while(j>=Seq_Dep&&peak[j]<=peak[j+step_length]){
                     j-=step_length;
                 }
                 if(j>=Seq_Dep){
                    j-=step_length;
                     while(j>=Seq_Dep&&peak[j]<=peak[j+step_length]){
                         j-=step_length;
                     }
                 }
             }

         //i-j大于等于多少具体值需要调一下
         if(i-j>=3*step_length){
             //1
             int judge=0;
             int kmer_dep=0;
             int count_dep_kmer=0;
             for(int k=max_kmerfre;k>2;k--){
                 if(count_kmer[k]>0&&count_kmer[k]<max_kmerfre){
                     kmer_dep+=count_kmer[k];
                     count_dep_kmer++;
                 }
             }
             cout<<"去除超低复杂度后k-mer总量为"<<kmer_dep<<endl;
             cout<<"去除低频k-mer后的平均k-mer厚度为"<<kmer_dep/count_dep_kmer<<endl;


             //1
             if ((j-(max_kmerfre-j))>0)
                 peak_start=j-(max_kmerfre-j);
             else
                 peak_start=j;
             //先到max_fre
             peak_end=max_kmerfre;
             first_p=j-(max_kmerfre-j);
             second_p=max_kmerfre;
             dep_search= 0;
             cout<<"寻找到的最后一个峰值区间峰顶为"<<first_p<<endl;
             //2
             for(map<string,int >::iterator kmer_it=kmer_map.begin();kmer_it!=kmer_map.end();kmer_it++){
                 if((kmer_it->second)>first_p){
                     dep_judge++;
                 }
             }
             cout<<"寻峰k-mer数量"<<dep_judge<<endl;
             cout<<"寻峰k-mer数量"<<kmer_dep/count_dep_kmer<<endl;
             if(max_kmerfre<5000&&(kmer_dep/count_dep_kmer)<25){
                 dep_search= 1;
                 for(int k=max_kmerfre;k>2;k--){
                     if(count_kmer[k]>0&&count_kmer[k]<max_kmerfre)
                         judge+=count_kmer[k];
                     if(judge>(kmer_dep/4*3)){
                         cout<<"识别到低测序深度，使用备选方案选取公共kmer,起始位置为"<<k<<",当前选取数量"<<judge<<endl;
                         peak_start=k;
                         first_p=k;

//                         peak_start=33;
//                         first_p=33;

                         peak_end=max_kmerfre;
                         second_p=max_kmerfre;
                         return true;
                     }
                 }

             }
             //2

             //3
             if(max_kmerfre>10000)
                 dep_search= 2;
             //3

             cout<<"识别到高测序深度，寻峰K-mer起始位置"<<peak_start<<endl;
             return true;
         }
     }
return false;
}
//生成Kmer
void get_kmer(int kmer_length,map<string,int > & data,vector<int >data_seq,vector<string >data_str,int& dep_search) {

    int data_size = data_str.size();
    //it1迭代器用于寻找map中是否含有该kmer
    map<string ,int>::iterator it1;
    cout << "kmer map contructing ..." << endl;
    time_t start_time = time(NULL);
//.empty（）为C++自带函数，判断是否为空，为空返回1,不空返回0
    if (data_str.empty()) {
        cout << "data empty!!" << endl;
        return;
    }
    //生成kmer
    //i控制当前遍历到哪一个vector中的元素，即控制读取哪一条read
    //j控制当前read生成kmer，使用了.substr(起始位置,长度)方法截取Kmer
    for ( int num=0;num<data_str.size();num++) {
        //cout<<"kmer"<<endl;
        const string& read = data_str[num];
        {
            if (read.length() <= kmer_length)
                continue;
        }
        for (int j = 0; j <= read.length()-kmer_length;j++) {
            const string& kmer = read.substr(j, kmer_length);


            if (!contains_non_gatc(kmer)) {
                continue;
            }
            //判断map容器中是否含有该kmer，如果没有则添加该元素，并初始化值为1，如果有，则map_second+1;

            if (kmer_map.find(kmer)==kmer_map.end())
            {
                kmer_map.insert(make_pair(kmer,1));
            }
            else
            {
                kmer_map[kmer]+=data_dep[num];
            }
        }
    }
    time_t kmer_map_time = time(NULL);
    cout << " kmer_map" << (kmer_map_time - start_time) << " s)" << endl;
    cout<<"k-mer频次数据k-mer-map生成完毕,k-mer个数为"<<kmer_map.size()<<endl;
    int aver_kmer=0;
    int read_num=0;
    for(int x=0;x<data_str.size();x++){
        aver_kmer+=data_str[x].size()*data_seq[x];
    }
    aver_kmer=aver_kmer/kmer_map.size();
    cout<<"k-mer平均厚度！！！！！！！！！！！！！！！！！！"<<aver_kmer<<endl;
    cout<<"公共k-mer-read映射关系构建中"<<endl;
//单独写一个函数
//阈值为频次取峰值，超过多大频率的kmer为公共kmer，某个范围之内。写代码求出来，最后一个峰值的区间，左右多少的取值需要调整，算出来之后就能确定
    //寻找公共kmer,寻找这些公共kmer出现过的reads
    //建立新的映射信息
    //现在有了所有公共kmer，存储在candi_K这个容器中，现在建立映射关系
    //存储com_kmer-read映射关系

    //中间变量存储read信息
    //vector<int >read_information;
    //存储每条read中包含的公共kmer个数,初始化read个数，每个都初始化值为0
//    vector<int>read_kmer_num (data.size(),0);
    //统计包含公共kmer个数大于2的read数量
    int more_2k_read=0;

//更改遍历条件

    int peak_start=0;
    int peak_end=0;
    int dep;
    if(!Find_peak(kmer_map,data,peak_start,peak_end,data_dep,data_str,dep,dep_search)){
        cout<<"can not find peak!"<<endl;
        return;
    }
        time_t find_peak_time = time(NULL);
    cout << " kmer_map" << (find_peak_time - start_time) << " s)" << endl;
    Seq_Dep=dep;

    cout<<"峰值区间起始位置"<<peak_start<<"峰值区间结束位置"<<peak_end<<endl;
    cout<<"平均k-mer厚度"<<dep<<endl;

    cout<<"将要遍历"<<data_str.size()<<"条read"<<endl;
    for ( int i=0;i<data_str.size();i++) {
        // if(i%10000==0)
        //     cout<<"当前进行到"<<i<<endl;
        //cout<<"kmer"<<endl;
        const string& read = data_str[i];//更改

            if (read.length() <= kmer_length)
                continue;

        for (int j = 0; j <= read.length()-kmer_length;j++) {
            const string& kmer = read.substr(j, kmer_length);
            if(kmer_map[kmer]>=peak_start&&kmer_map[kmer]<=peak_end){//自动识别一个区间如果不好确定，首先去除频次很高的k-mer 去除频次为1的k-mer 去除后取总k-mer个数的前多少。kmer平均多厚：read个数X长度/
                //read_information.push_back(j);//更改  在这个for循环前定义一个i，每次循环i++
                k_r_p[kmer].push_back(make_pair(i,j));
            }
        }
        // if(read_information.size()>=2){
        //     for(int k=0;k<read_information.size();k++){
        //         k_r_p[read.substr(read_information[k],kmer_length)].push_back(make_pair(i,read_information[k]));

        //     }
        //     more_2k_read++;
        // }
        // read_information.clear();
    }
cout<<"k_r_p success found! size is "<<k_r_p.size()<<endl;
    time_t end_time = time(NULL);
    cout << "Kmer_map has been constructed, total " << kmer_map.size() << " kmers! time total: " << (end_time - start_time) << " s)" << endl;

}



//已经获得了所有拥有两条公共kmer的read及其公共kmer
//1.先找出拥有公共Kmer最多的一条Read作为模板
//从左到右遍历该条read，如果kmer重叠则不进行比对，直到找到两条不重叠的k-mer
//传入数据为k_r_p信息，kmer频次信息，read中包含公共kmer个数,所有read






