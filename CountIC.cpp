#define _WIN32_WINNT 0x0A00
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#pragma comment(lib, "psapi.lib")
#else
#include <sys/resource.h>
#include <unistd.h>
#endif
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <set>
#include <map>
#include <ctime>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <bitset>
#include "read_dataset.cpp"
#include <random>

using namespace std;


////////////////////////////////////////////////////      helper function   from line     //////////////////////////////////////////////////////////////

struct Compare {
    const std::vector<std::vector<int>>& hyperedge2node;

    Compare(const std::vector<std::vector<int>>& hyperedge2node)
            : hyperedge2node(hyperedge2node) {}

    bool operator()(int lhs, int rhs) const {
        if (hyperedge2node[lhs].size() == hyperedge2node[rhs].size()) {
            return lhs < rhs;
        }
        return hyperedge2node[lhs].size() < hyperedge2node[rhs].size();
    }
};


map<int,vector<vector<int>>> result;
vector<long long> h_motif;
vector< vector<int> > node2hyperedge;

vector< vector<int> > hyperedge2node;
vector< unordered_set<int> > hyperedge2node_set;
vector< unordered_map<int,int>> PairTypeInc; //<<id1, id2>,intersection>, <<id1, id2>, inclusion>.... note that we assume the pair.first < pair.second
vector< unordered_map<int,vector<int>> > PairTypeInt; //<<id1, id2>,intersection>, <<id1, id2>, inclusion>.... note that we assume the pair.first < pair.second
vector<vector<int>> HyperEdgeConnectionInc; //<<id2>.....> HyperEdge id1 and id2 are connect through inclusion type
vector<vector<int>> HyperEdgeConnectionInt; //<<id2>.....> HyperEdge id1 and id2 are connect through intersection type
vector< vector<int>> InclusionList;
vector< unordered_set<int> >  TotalSame;
vector< vector<int> > HyperEdgeConnectionRefinedInc;
vector< vector<int> > HyperEdgeConnectionRefinedInt;
vector< vector<int> > HyperEdgeConnectionRefinedInt2;
vector< vector<int>> EnclosureList;
vector< vector<int>> IntersectionList;
//vector< unordered_set<int>> IntersectionSet;
vector<int> sortedHyperEdge;
//vector< unordered_map<int,int>> PairCheck;
vector<unordered_map<int,vector<bool>>> PairCount;
vector<unordered_map<int,int>> PairId;
int TotalPair=0;
vector<bool> checkBefore;
vector<unordered_set<int>> NotConnected;

vector< unordered_set<int>> IncludePair; //<<id1, id2>,intersection>, <<id1, id2>, inclusion>.... note that we assume the pair.first < pair.second
vector< unordered_map<int,vector<int>> > IntersecPair; //<<id1, id2>,intersection>, <<id1, id2>, inclusion>.... note that we assume the pair.first < pair.second

long long inclusion_pairs = 0;
long long intersection_pairs = 0;
vector<int> findCommonElements(vector<int>& vec1, vector<int>& vec2) {
    vector<int> result;
    int i = 0;
    int j = 0;

    while (i < vec1.size() && j < vec2.size()) {
        if (vec1[i] < vec2[j]) {
            ++i;
        } else if (vec1[i] > vec2[j]) {
            ++j;
        } else {
            result.push_back(vec1[i]);
            ++i;
            ++j;
        }
    }

    return result;
}

bool compareHyperEdge(const std::pair<int, std::vector<int>>& a, const std::pair<int, std::vector<int>>& b) {
    if (a.second.size() != b.second.size()) {
        return a.second.size() < b.second.size();
    }
    return a.first < b.first;
}

bool AltCompareHyperEdge(int a, int b) {
    if (hyperedge2node[a].size() == hyperedge2node[b].size()) {
        return a < b;
    }
    return hyperedge2node[a].size() < hyperedge2node[b].size();
}

std::vector<int> sortHyperedge2Node(const std::vector<std::vector<int>>& hyperedge2node) {
    std::vector<std::pair<int, std::vector<int>>> indexedHyperedge2Node;

    // 将每个vector<int>与其index存储到一个pair中
    for (int i = 0; i < hyperedge2node.size(); ++i) {
        indexedHyperedge2Node.emplace_back(i, hyperedge2node[i]);
    }

    // 使用compare函数进行排序
    std::sort(indexedHyperedge2Node.begin(), indexedHyperedge2Node.end(), compareHyperEdge);

    // 构建排序后的结果
    std::vector<int> result;
    for (const auto& pair : indexedHyperedge2Node) {
        result.push_back(pair.first);
    }

    return result;
}

vector<int>* findPair(vector< unordered_map<int,vector<int>> >& PairTypeInt, int id1, int id2){
    if(id1<id2){
        auto it = PairTypeInt[id1].find(id2);
        if(it!=PairTypeInt[id1].end()){
            return &it->second;
        }

    }else{
        auto it = PairTypeInt[id2].find(id1);
        if(it!=PairTypeInt[id2].end()){
            return &it->second;
        }

    }
    return NULL;
}

int findPairSize(vector< unordered_map<int,vector<int>> >& PairTypeInt, int id1, int id2){
    if(id1<id2){
        auto it = PairTypeInt[id1].find(id2);
        if(it!=PairTypeInt[id1].end()){
            return it->second.size();
        }

    }else{
        auto it = PairTypeInt[id2].find(id1);
        if(it!=PairTypeInt[id2].end()){
            return it->second.size();
        }

    }
    return 0;
}

vector<int>* findPair2(vector< unordered_map<int,int >>& PairTypeInc, vector< vector<int> >& hyperedge2node, int id1, int id2){
    if(id1<id2){
        auto it = PairTypeInc[id1].find(id2);
        if(it==PairTypeInc[id1].end()){
            return NULL;
        }
        int includId=it->second;
        return &hyperedge2node[includId];
    }else{
        auto it = PairTypeInc[id2].find(id1);
        if(it==PairTypeInc[id2].end()){
            return NULL;
        }
        int includId=it->second;
        return &hyperedge2node[includId];
    }
}

int findPairSize2(vector< unordered_map<int,int >>& PairTypeInc, vector< vector<int> >& hyperedge2node, int id1, int id2){
    if(id1<id2){
        auto it = PairTypeInc[id1].find(id2);
        if(it!=PairTypeInc[id1].end()){
            int includId=it->second;
            return hyperedge2node[includId].size();
        }

    }else{
        auto it = PairTypeInc[id2].find(id1);
        if(it!=PairTypeInc[id2].end()){
            int includId=it->second;
            return hyperedge2node[includId].size();
        }

    }
    return 0;

}

int DistinguishAAA(int part1, int part2, int part3, int part12, int part13, int part23, int part123){
    int count=0;
    if(part12!=0){
        count+=1;
    }
    if(part23!=0){
        count+=1;
    }
    if(part13!=0){
        count+=1;
    }

    if(part1!=0 && part2!=0 && part3!=0){//motif 2 6 12 16 26
        if(part123==0){
            return 26;
        }
        if(count==0){
            return 2;
        }else if(count==1){
            return 6;
        }else if(count==2){
            return 12;
        }else if(count==3){
            return 16;
        }
    }else if((part1==0 && part2!=0 && part3!=0) || (part1!=0 && part2==0 && part3!=0) ||(part1!=0 && part2!=0 && part3==0) ){//11 15 25
        if(part123==0){
            return 25;
        }
        if(count==2){
            return 11;
        }else if(count==3){
            return 15;
        }
    }else if(part1==0 && part2==0 && part3==0){//13 23
        if(part123==0){
            return 23;
        }else{
            return 13;
        }
    }else{// 14 24
        if(part123==0){
            return 24;
        }else{
            return 14;
        }
    }
    return -1;
}

bool CheckSame(vector< unordered_set<int> >& TotalSame, int id1, int id2){
    if(id1<id2){
        auto it = TotalSame[id1].find(id2);
        if(it==TotalSame[id1].end()){
            return true;
        }
        return false;
    }else{
        auto it = TotalSame[id2].find(id1);
        if(it==TotalSame[id2].end()){
            return true;
        }
        return false;
    }

}

////////////////////////////////////////////////////   end of helper function     //////////////////////////////////////////////////////////////


void preprocess(int currentId,vector< unordered_set<int> >& TotalSame, vector< vector<int>>& InclusionList, vector< vector<int>>& HyperEdgeConnectionRefinedInc, vector< vector<int>>& HyperEdgeConnectionRefinedInt, vector< vector<int>>& HyperEdgeConnectionRefinedInt2, vector< vector<int> >& hyperedge2node, vector< vector<int> >& HyperEdgeConnectionInc, vector< vector<int> >& HyperEdgeConnectionInt, vector< unordered_map<int,int> >& PairTypeInc, vector< unordered_map<int,vector<int>> >& PairTypeInt, vector<vector<int>>& Hyperedge2Node2Hyperedge, vector<int>& NodeIdList){
    while(true){
        int id=-1;
        bool check=false;
        for(size_t i=0; i<Hyperedge2Node2Hyperedge.size();i++){
            if(Hyperedge2Node2Hyperedge[i].size()>0){
                check=true;
                if(id==-1){
                    id=Hyperedge2Node2Hyperedge[i][0];
                }
                if(id>Hyperedge2Node2Hyperedge[i][0]){
                    id=Hyperedge2Node2Hyperedge[i][0];
                }
            }
        }
        if(!check){
            break;
        }
        int TtoalNumber=0;
        for(size_t i=0; i<Hyperedge2Node2Hyperedge.size();i++){
            if(Hyperedge2Node2Hyperedge[i].size()>0){
                if(id==Hyperedge2Node2Hyperedge[i][0]){
                    //Hyperedge2Node2Hyperedge[i].erase(Hyperedge2Node2Hyperedge[i].begin());
                    TtoalNumber+=1;
                }
            }
        }

        int MinimumId=id;
        int EdgeNodeNum1=hyperedge2node[currentId].size();
        int EdgeNodeNum2=hyperedge2node[MinimumId].size();
        int InclusionNum=TtoalNumber;

        if((EdgeNodeNum1>EdgeNodeNum2 && EdgeNodeNum2==InclusionNum) || (EdgeNodeNum1<EdgeNodeNum2 && EdgeNodeNum1==InclusionNum) ){//hyperedge2 is included in hyperedge1 or hyperedge1 is included in hyperedge2
            pair<int, int> pair2 = make_pair(currentId, MinimumId);
            //PairTypeInc.insert(make_pair(pair2,pair1.second));
            //StoreEdgeType(MinimumId, currentId, HyperEdgeType, "inclusion");
            for(size_t i=0; i<Hyperedge2Node2Hyperedge.size();i++){
                if(Hyperedge2Node2Hyperedge[i].size()>0){
                    if(id==Hyperedge2Node2Hyperedge[i][0]){
                        Hyperedge2Node2Hyperedge[i].erase(Hyperedge2Node2Hyperedge[i].begin());
                    }
                }
            }
            if(EdgeNodeNum1>EdgeNodeNum2){
                PairTypeInc[currentId].insert(make_pair(MinimumId,MinimumId));
                InclusionList[MinimumId].push_back(currentId);
                EnclosureList[currentId].push_back(MinimumId);
                IncludePair[MinimumId].insert(currentId);
            }else{
                PairTypeInc[currentId].insert(make_pair(MinimumId,currentId));
                InclusionList[currentId].push_back(MinimumId);
                EnclosureList[MinimumId].push_back(currentId);
                IncludePair[currentId].insert(MinimumId);
            }
            inclusion_pairs++;
            HyperEdgeConnectionInc[currentId].push_back(MinimumId);
            HyperEdgeConnectionInc[MinimumId].push_back(currentId);
            HyperEdgeConnectionRefinedInc[currentId].push_back(MinimumId);
        }else if(InclusionNum!=EdgeNodeNum2 && InclusionNum!=EdgeNodeNum1){//intersection
            PairTypeInt[currentId].insert(make_pair(MinimumId,vector<int>()));
            TotalPair+=1;
            intersection_pairs++;
            //IntersectionSet[currentId].insert(MinimumId);


            //PairCheck[currentId].insert(make_pair(MinimumId,0));

            if(AltCompareHyperEdge(currentId, MinimumId)){
                IntersectionList[currentId].push_back(MinimumId);
                IntersecPair[currentId].insert(make_pair(MinimumId,vector<int>()));
                //PairCount[currentId].insert(make_pair(MinimumId,vector<bool>(hyperedge2node.size()+1, false)));

                //PairId[currentId].insert({MinimumId, TotalPair-1});

                IntersecPair[currentId][MinimumId].push_back(TotalPair-1);
            }else{
                IntersectionList[MinimumId].push_back(currentId);
                IntersecPair[MinimumId].insert(make_pair(currentId,vector<int>()));
                IntersecPair[MinimumId][currentId].push_back(TotalPair-1);
            }

            for(size_t i=0; i<Hyperedge2Node2Hyperedge.size();i++){
                if(Hyperedge2Node2Hyperedge[i].size()>0){
                    if(id==Hyperedge2Node2Hyperedge[i][0]){
                        Hyperedge2Node2Hyperedge[i].erase(Hyperedge2Node2Hyperedge[i].begin());
                        PairTypeInt[currentId][MinimumId].push_back(NodeIdList[i]);
                        if(AltCompareHyperEdge(currentId, MinimumId)){
                            IntersecPair[currentId][MinimumId].push_back(NodeIdList[i]);
                        }else{
                            IntersecPair[MinimumId][currentId].push_back(NodeIdList[i]);
                        }
                    }
                }
            }
            HyperEdgeConnectionInt[currentId].push_back(MinimumId);
            HyperEdgeConnectionInt[MinimumId].push_back(currentId);
            HyperEdgeConnectionRefinedInt[currentId].push_back(MinimumId);
            HyperEdgeConnectionRefinedInt2[currentId].push_back(MinimumId);

        }else if(InclusionNum==EdgeNodeNum2 && InclusionNum==EdgeNodeNum1){

            for(size_t i=0; i<Hyperedge2Node2Hyperedge.size();i++){
                if(Hyperedge2Node2Hyperedge[i].size()>0){
                    if(id==Hyperedge2Node2Hyperedge[i][0]){
                        Hyperedge2Node2Hyperedge[i].erase(Hyperedge2Node2Hyperedge[i].begin());
                    }
                }
            }
            TotalSame[currentId].insert(MinimumId);
        }else{
            for(size_t i=0; i<Hyperedge2Node2Hyperedge.size();i++){
                if(Hyperedge2Node2Hyperedge[i].size()>0){
                    if(id==Hyperedge2Node2Hyperedge[i][0]){
                        Hyperedge2Node2Hyperedge[i].erase(Hyperedge2Node2Hyperedge[i].begin());
                    }
                }
            }
        }
    }
}



int CountIII(){
    for(int i=0; i<sortedHyperEdge.size();i++){
        int currentEdge=sortedHyperEdge[i];

        for(int j=0;j<IntersectionList[currentEdge].size();j++){
            int secondEdge=IntersectionList[currentEdge][j];
            vector<int>* CommonVer=&IntersecPair[currentEdge][secondEdge];
            int FirstId=(*CommonVer)[0];
            for(int k=0;k<IntersectionList[secondEdge].size();k++){
                int thirdEdge=IntersectionList[secondEdge][k];
                if(AltCompareHyperEdge(IntersectionList[currentEdge][IntersectionList[currentEdge].size()-1], thirdEdge)){
                    break;
                }
                auto it2 = IntersecPair[currentEdge].find(thirdEdge);
                if(it2== IntersecPair[currentEdge].end() ){
                    continue;
                }
                vector<int>* firstCommon=&it2->second;
                auto it3 = IntersecPair[secondEdge].find(thirdEdge);
                vector<int>* secondCommon=&it3->second;
                int TripleCommon=0;
                int ii = 1, jj = 1;
                int n1 = firstCommon->size(), n2 = secondCommon->size();
                while (ii < n1 && jj < n2 ) {
                    if ((*firstCommon)[ii] == (*secondCommon)[jj] ) {
                        TripleCommon+=1;
                        ii++;
                        jj++;

                    }
                    else if ((*firstCommon)[ii] < (*secondCommon)[jj]) {
                        ii++;
                    }
                    else {
                        jj++;
                    }

                }
                int part123=TripleCommon;
                int part12= CommonVer->size()-1-part123;
                int part13=firstCommon->size()-1-part123;
                int part23=secondCommon->size()-1-part123;
                int part1=hyperedge2node[currentEdge].size()-part12-part13-part123;
                int part2=hyperedge2node[secondEdge].size()-part12-part23-part123;
                int part3=hyperedge2node[thirdEdge].size()-part23-part13-part123;
                int motifNum=DistinguishAAA(part1,part2, part3, part12, part13, part23, part123);
                h_motif[motifNum-1]+=1;
            }
        }
    }
    return 0;
}

int CountIIC(){

    for (int ii=0; ii<PairTypeInc.size();ii++) {
        for(const auto& pair : PairTypeInc[ii]){

            int firstId=ii;
            int secondId=pair.first;
            int includeId=pair.second;
            vector<int> CommonVer = hyperedge2node[includeId];
            if(HyperEdgeConnectionInt[firstId].size()==0 && HyperEdgeConnectionInt[secondId].size()==0 ){//this two hyperedges only contian in inclusion type
                continue;
            }
            int i = 0;
            int j = 0;
            while (i < HyperEdgeConnectionInt[firstId].size() || j < HyperEdgeConnectionInt[secondId].size()) {
                if(i == HyperEdgeConnectionInt[firstId].size()){
                    int thirdId=HyperEdgeConnectionInt[secondId][j];
                    ++j;
                    if(!CheckSame(TotalSame, firstId, thirdId)){
                        continue;
                    }
                    if(thirdId==firstId){
                        continue;
                    }
                    int Common2=findPairSize2(PairTypeInc, hyperedge2node, thirdId, firstId);
                    if(Common2!=0){
                        continue;
                    }
                    int firstCommonSize=findPairSize2(PairTypeInc,hyperedge2node, firstId, secondId);
                    int secondCommonSize=findPairSize(PairTypeInt, thirdId, secondId);
                    if(firstCommonSize+secondCommonSize==hyperedge2node[secondId].size()){
                        h_motif[19-1]+=1;
                    }else{
                        h_motif[20-1]+=1;
                    }
                }else if(j == HyperEdgeConnectionInt[secondId].size()){
                    int thirdId=HyperEdgeConnectionInt[firstId][i];
                    ++i;
                    if(thirdId==secondId){
                        continue;
                    }

                    if(!CheckSame(TotalSame, secondId, thirdId)){
                        continue;
                    }
                    int Common2=findPairSize2(PairTypeInc, hyperedge2node, thirdId, secondId);
                    if(Common2!=0){
                        continue;
                    }
                    int firstCommonSize=findPairSize2(PairTypeInc,hyperedge2node, firstId, secondId);
                    int secondCommonSize=findPairSize(PairTypeInt, thirdId, firstId);
                    if(firstCommonSize+secondCommonSize==hyperedge2node[firstId].size()){
                        h_motif[19-1]+=1;
                    }else{
                        h_motif[20-1]+=1;
                    }
                }else{
                    if (HyperEdgeConnectionInt[firstId][i] < HyperEdgeConnectionInt[secondId][j]) {
                        int thirdId=HyperEdgeConnectionInt[firstId][i];
                        ++i;
                        if(thirdId==secondId){
                            continue;
                        }
                        if(!CheckSame(TotalSame, secondId, thirdId)){
                            continue;
                        }
                        int Common2=findPairSize2(PairTypeInc, hyperedge2node, thirdId, secondId);
                        if(Common2!=0){
                            continue;
                        }
                        int firstCommonSize=findPairSize2(PairTypeInc,hyperedge2node, firstId, secondId);
                        int secondCommonSize=findPairSize(PairTypeInt, thirdId, firstId);
                        if(firstCommonSize+secondCommonSize==hyperedge2node[firstId].size()){
                            h_motif[19-1]+=1;
                        }else{
                            h_motif[20-1]+=1;
                        }

                    } else if (HyperEdgeConnectionInt[firstId][i] > HyperEdgeConnectionInt[secondId][j]) {
                        int thirdId=HyperEdgeConnectionInt[secondId][j];
                        ++j;
                        if(!CheckSame(TotalSame, firstId, thirdId)){
                            continue;
                        }
                        if(thirdId==firstId){
                            continue;
                        }
                        int Common2=findPairSize2(PairTypeInc, hyperedge2node, thirdId, firstId);
                        if(Common2!=0){
                            continue;
                        }
                        int firstCommonSize=findPairSize2(PairTypeInc,hyperedge2node, firstId, secondId);
                        int secondCommonSize=findPairSize(PairTypeInt, thirdId, secondId);
                        if(firstCommonSize+secondCommonSize==hyperedge2node[secondId].size()){
                            h_motif[19-1]+=1;
                        }else{
                            h_motif[20-1]+=1;
                        }
                    } else {
                        int commonHyperedgeId=HyperEdgeConnectionInt[firstId][i];
                        vector<int>* firstCommon=findPair(PairTypeInt, commonHyperedgeId, firstId);
                        vector<int>* secondCommon=findPair(PairTypeInt, commonHyperedgeId, secondId);
                        if(firstCommon->size()!=secondCommon->size()){//belong to motif type 9, 10
                            int ABC;
                            int BC;
                            int tempId;
                            if(firstCommon->size()<secondCommon->size()){
                                ABC=firstCommon->size();
                                tempId=secondId;
                                BC= secondCommon->size();
                            }else{
                                ABC=secondCommon->size();
                                tempId=firstId;
                                BC=firstCommon->size();
                            }
                            if(BC+CommonVer.size()-ABC==hyperedge2node[tempId].size()){//belong to motif type 9
                                //result[8].push_back(motifEdges);
                                h_motif[8]+=1;
                            }else{//belong to motif type 10
                                //result[9].push_back(motifEdges);
                                h_motif[9]+=1;
                            }
                        }else {
                            //result[4].push_back(motifEdges);
                            h_motif[4]+=1;
                        }
                        ++i;
                        ++j;
                    }
                }
            }
        }
    }
    return 0;
}

int CountCCC() {
    for (int ii = 0; ii < PairTypeInc.size(); ii++) {
        for (const auto& pair : PairTypeInc[ii]) {
            int firstId = ii;
            int secondId = pair.first;
            for (int i = 0; i < HyperEdgeConnectionRefinedInc[secondId].size(); i++) {
                int commonHyperedgeId = HyperEdgeConnectionRefinedInc[secondId][i];
                if (commonHyperedgeId >
                    HyperEdgeConnectionRefinedInc[firstId][HyperEdgeConnectionRefinedInc[firstId].size() - 1]) {
                    break;
                }
                auto it = PairTypeInc[firstId].find(commonHyperedgeId);
                if (it != PairTypeInc[firstId].end()) {
                    h_motif[2] += 1; // 模体3
                }
            }
        }
    }
    return 0;
}

int CountICC(){

    for (int ii=0; ii<PairTypeInt.size();ii++) {
        for(const auto& pair : PairTypeInt[ii]){
            int firstId=ii;
            int secondId=pair.first;
            vector<int> CommonVer = pair.second;
            if(HyperEdgeConnectionInc[firstId].size()==0 || HyperEdgeConnectionInc[secondId].size()==0 ){//this two hyperedges only contian in inclusion type
                continue;
            }
            if(HyperEdgeConnectionInc[firstId][HyperEdgeConnectionInc[firstId].size()-1]<HyperEdgeConnectionInc[secondId][0] || HyperEdgeConnectionInc[firstId][0]>HyperEdgeConnectionInc[secondId][HyperEdgeConnectionInc[secondId].size()-1]){
                continue;
            }
            int i = 0;
            int j = 0;
            while (i < HyperEdgeConnectionInc[firstId].size() && j < HyperEdgeConnectionInc[secondId].size()) {
                if (HyperEdgeConnectionInc[firstId][i] < HyperEdgeConnectionInc[secondId][j]) {
                    ++i;
                } else if (HyperEdgeConnectionInc[firstId][i] > HyperEdgeConnectionInc[secondId][j]) {
                    ++j;
                } else {
                    int commonHyperedgeId=HyperEdgeConnectionInc[firstId][i];
                    vector<int>* firstCommon=findPair2(PairTypeInc, hyperedge2node, commonHyperedgeId, firstId);
                    vector<int>* secondCommon=findPair2(PairTypeInc, hyperedge2node,commonHyperedgeId, secondId);
                    vector<int> motifEdges;
                    motifEdges.push_back(firstId);
                    motifEdges.push_back(secondId);
                    motifEdges.push_back(commonHyperedgeId);
                    if(hyperedge2node[firstId].size()>hyperedge2node[commonHyperedgeId].size() && hyperedge2node[secondId].size()>hyperedge2node[commonHyperedgeId].size()){//h-motif 1 or 4
                        if(CommonVer.size()==firstCommon->size()){
                            h_motif[0]+=1;
                        }else{
                            h_motif[3]+=1;
                        }
                    }else{//h-motif 7 or 8
                        int countNum=firstCommon->size()+secondCommon->size()-CommonVer.size();
                        if(countNum==hyperedge2node[commonHyperedgeId].size()){
                            h_motif[6]+=1;
                        }else{
                            h_motif[7]+=1;
                        }
                    }
                    ++i;
                    ++j;
                }
            }
        }
    }
    return 0;
}

int CountII()
{
    for (int a = 0; a < (int)PairTypeInt.size(); ++a)
    {
        for (const auto& kv : PairTypeInt[a])
        {
            int b = kv.first;
            if (b <= a) continue;
            for (int c : IntersectionList[b])
            {
                if (c <= b) continue;
                if (findPairSize(PairTypeInt, a, c) != 0) continue;
                int ab = (int)kv.second.size();
                int bc = findPairSize(PairTypeInt, b, c);
                int sz_b = (int)hyperedge2node[b].size();
                if (sz_b == ab + bc)
                    h_motif[21 - 1] += 1;
                else if (sz_b > ab + bc)
                    h_motif[22 - 1] += 1;
            }
        }
    }
    return 0;
}


int CountCC()
{

    for (int parent = 0; parent < (int)PairTypeInc.size(); ++parent)
    {
        const vector<int>& children = EnclosureList[parent];
        if (children.size() < 2) continue;
        for (size_t idx1 = 0; idx1 < children.size(); ++idx1)
        {
            int b = children[idx1];
            for (size_t idx2 = idx1 + 1; idx2 < children.size(); ++idx2)
            {
                int c = children[idx2];
                bool disjoint = true;
                const auto& set_b = hyperedge2node_set[b];
                const auto& vec_c = hyperedge2node[c];
                for (int v : vec_c)
                    if (set_b.count(v)) { disjoint = false; break; }
                if (!disjoint) continue;
                int sz_parent = (int)hyperedge2node[parent].size();
                int sz_b       = (int)hyperedge2node[b].size();
                int sz_c       = (int)hyperedge2node[c].size();
                int sum_bc     = sz_b + sz_c;
                if (sz_parent == sum_bc)
                    h_motif[17 - 1] += 1;
                else if (sz_parent > sum_bc)
                    h_motif[18 - 1] += 1;
            }
        }
    }
    return 0;
}
int Inclusion_Motif(){

    CountICC();

    CountCCC();

    CountCC();

    return 1;
}
int Intersection_Motif(){

    CountIII();

    CountIIC();

    CountII();

    return 1;
}


inline long long convert_id(int hyperedge_a, int hyperedge_b){
    return hyperedge_a * (1LL << 31) + hyperedge_b;
}

int main(int argc, char *argv[])
{
    clock_t start;
    clock_t run_start;
    int progress;

    //for test
    map<int,vector<vector<int>>> result2;
    vector<long long> h_triangle2;

    for(int i=0; i<26; i++){
        vector<vector<int>> temp;
        result.insert(make_pair(i,temp));
        h_motif.push_back(0);
        vector<vector<int>> temp2;
        result2.insert(make_pair(i,temp2));
        h_triangle2.push_back(0);
    }

    string graphFile = (argc > 2) ? argv[2] : "email-Enron-full.txt";

    // Read data
    start = clock();

    read_data(graphFile, node2hyperedge, hyperedge2node, hyperedge2node_set);


    int V = node2hyperedge.size(), E = hyperedge2node.size();
    cout << "# of nodes: " << V << '\n';
    cout << "# of hyperedges: " << E << '\n';
    cout << "Reading data done: "
         << (double)(clock() - start) / CLOCKS_PER_SEC << " sec" << endl;
    cout << "------------------------------------------" << endl << endl;
    /* =====  重叠度 Overlapness = sum(|e|) / |V|  ===== */
    long long sumEdgeSize = 0;
    for (const auto &e : hyperedge2node) sumEdgeSize += e.size();
    cout << "Overlapness: " << fixed << setprecision(6)
         << (double)sumEdgeSize / V << endl;
    cout << "Density: " << fixed << setprecision(6)
         << (double)E / V << endl;
    run_start = clock();

    vector<vector<int>> copyVec=node2hyperedge;
    vector<int> NodeIdList;

    for(int i=0; i<hyperedge2node.size();i++){
        PairTypeInc.push_back(unordered_map<int,int>());
        PairTypeInt.push_back(unordered_map<int,vector<int>> ());
        HyperEdgeConnectionInc.push_back(vector<int>());
        HyperEdgeConnectionInt.push_back(vector<int>());
        TotalSame.push_back(unordered_set<int> ());
        HyperEdgeConnectionRefinedInc.push_back(vector<int>());
        HyperEdgeConnectionRefinedInt.push_back(vector<int>());
        HyperEdgeConnectionRefinedInt2.push_back(vector<int>());
        InclusionList.push_back(vector<int>());
        EnclosureList.push_back(vector<int>());
        IntersectionList.push_back(vector<int>());
        PairCount.push_back(unordered_map<int,vector<bool>>());
        IncludePair.push_back(unordered_set<int>());
        IntersecPair.push_back(unordered_map<int,vector<int>> ());
        PairId.push_back(unordered_map<int,int>());
    }
    start = clock();

    for(int id=0; id<hyperedge2node.size();id++){
        /*if (id % 10 == 0)
            cout << "Hyperedge: " << id << " / " << hyperedge2node.size() << endl;*/
        vector<int> connectNode=hyperedge2node[id];
        vector<vector<int>> Hyperedge2Node2Hyperedge;
        for(size_t i=0; i<connectNode.size();i++){
            int NodeId=connectNode[i];
            if(copyVec[NodeId].size()==1){//this node should be removed
                copyVec[NodeId][0]=-1;
            }else{
                copyVec[NodeId].erase(copyVec[NodeId].begin());
                Hyperedge2Node2Hyperedge.push_back(copyVec[NodeId]);
                NodeIdList.push_back(NodeId);
            }
        }

        preprocess(id,TotalSame,InclusionList, HyperEdgeConnectionRefinedInc, HyperEdgeConnectionRefinedInt, HyperEdgeConnectionRefinedInt2, hyperedge2node, HyperEdgeConnectionInc, HyperEdgeConnectionInt, PairTypeInc, PairTypeInt, Hyperedge2Node2Hyperedge, NodeIdList);
        Hyperedge2Node2Hyperedge.clear();
        NodeIdList.clear();
    }

    clock_t time1=clock() - start;
    cout << "Adjacency list construction done: "
         << (double)(time1) / CLOCKS_PER_SEC << " sec" << endl;
    cout << "------------------------------------------" << endl << endl;
    clock_t start2 = clock();

    ////////////////////////////////////     build     hyperEdge    network      /////////////////////////
    sortedHyperEdge=sortHyperedge2Node(hyperedge2node);

    Compare comp(hyperedge2node);
    for (auto& vec : IntersectionList) {
        std::sort(vec.begin(), vec.end(), comp);
    }
    for (auto& vec : InclusionList) {
        std::sort(vec.begin(), vec.end(), comp);
    }

    clock_t time2=clock() - start2;
    cout << "build trees done: "
         << (double)(time2) / CLOCKS_PER_SEC << " sec" << endl;
    cout<<'\n';

    ///////////////////////// ///////////////////////// ///////////////////////// ///////////////////////// /////////////////////////

    start = clock();
    cout << "------------------------------------------" << endl;
    cout << "Number of inclusion hyperedge pairs: " << inclusion_pairs << endl;
    cout << "Number of intersection hyperedge pairs: " << intersection_pairs << endl;

    if (argc >= 2) {
        string flag = argv[1];
        if (flag == "Inclusion") {
            Inclusion_Motif();
        } else if (flag == "Intersection") {
            Intersection_Motif();
        } else {
            Inclusion_Motif();
            Intersection_Motif();
        }
    } else {
        Inclusion_Motif();
        Intersection_Motif();
    }
    double runtime = (double)(clock() - run_start) / CLOCKS_PER_SEC;
    clock_t time3=clock() - start;
    vector<long long> h_triangle_final(26, 0);
    long long total_motif=0;
    int index = 0;
    for (int i = 1; i <= 26; i++){
        cout << fixed << "h-motif " << ++index << ": " << fixed << h_motif[i-1] << endl;

        total_motif+=h_motif[i-1];

    }
    cout << "counting motif data done: "
         << (double)(time3) / CLOCKS_PER_SEC << " sec" << endl;
    cout << "Total runtime: " << runtime << endl;
    /* 峰值内存（KB） */
    auto getPeakRSS = []() -> long long {
#ifdef _WIN32
        PROCESS_MEMORY_COUNTERS pmc;
        return GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)) ?
               pmc.PeakWorkingSetSize / 1024 : 0;
#else
        struct rusage ru;
    return getrusage(RUSAGE_SELF, &ru) == 0 ? ru.ru_maxrss : 0;
#endif
    };
    cout << "Peak memory: " << getPeakRSS() << " KB" << endl;
    return 0;
}
