#ifdef __GNUC__
#include <unistd.h>
#endif

#include<iostream>
#include <stdio.h>
#include <list>
#include<sstream>
#include<cmath>
#include<algorithm>
#include<utility>
#include<numeric>

//headers
#include "timetrack.h"
#include "calcdb.h"
#include "eqclass.h"
#include "stats.h"
#include "chashtable.h"
#include "eclat.h"

using namespace std;

//global vars
char infile[300];
Dbase_Ctrl_Blk *DCB;
Stats stats;
unsigned long iCnt=0;
unsigned int TOP_K_EQ =50;// added by Karlton
int dist_type = 0;        // added by Karlton
int iNoOfDatasets = 1;    // added by Karlton
double dR = 1.0;          // added by Karlton
double dRho = 0.8;        // added by Karlton
bool bFindCentres=true;   // added by Karlton
int iXi=10;               // added by Karlton- # bins per dimension
double dAvgCluDim=0.15;    // added by Karlton

double MINSUP_PER = 0.05;
int MINSUPPORT=-1;
int DBASE_MAXITEM;
//default flags
bool output = false; //don't print freq itemsets
bool output_idlist = false; //don't print idlist

diff_vals diff_type = diff2; //default is diffset2 mode
diff_vals max_diff_type = nodiff; //default is to use tidsets for maxtest
sort_vals sort_type = incr; //default is to sort in increasing order
alg_vals alg_type = basicmax; //default is to find all freq patterns
closed_vals closed_type = chash; //default is not to eliminate non-closed
prune_vals prune_type = EXPONENTIAL; // Initialised to SCHISM

//extern functions
extern void enumerate_max_freq(vector<int>& itcnt,vector<Eqclass>& eqBest,Eqclass *eq, int iter, idlist &newmax,long DBASE_NUM_TRANS,double dThreshold);
extern void enumerate_max_closed_freq(vector<int>& itcnt,Eqclass *eq, int iter, idlist &newmax);
extern void enumerate_closed_freq(vector<int>& itcnt,Eqclass *eq, int iter, idlist &newmax);
extern void enumerate_freq(vector<int>& itcnt,Eqclass *eq, int iter);
extern void form_closed_f2_lists(Eqclass *eq);
extern void form_f2_lists(Eqclass *eq);
extern cHashTable hashtest; //for closed (with chash) testing

bool pruneTest(vector<int>& itcnt,int iSup, int iDim,long DBASE_NUM_TRANS=1000,Eqclass* neq=NULL, int iDiv=iXi) { // added by Karlton
  double dSup = double(iSup)/DBASE_NUM_TRANS;
  if(iDim==1) return (dSup >= 1.0/iDiv);
  double c=6.8114,a=.2325;
  int iDim2 = DBASE_MAXITEM/iDiv;
  double dDim = log(double(DBASE_NUM_TRANS))/log(double(iDiv));
  double dMean = 1.0;
  switch (dist_type) {
    case 0:
    dMean = pow(double(iDiv),double(-iDim));
    break;

    case 1:
    if(neq != (Eqclass *)NULL) {
      vector<int> vecPrefix=neq->prefix();
      for (unsigned int i=0; i < vecPrefix.size();){
        int iTmpMulti=(vecPrefix[i]/iDiv + 1)*iDiv;
        int iSum = 0;
        while(i < vecPrefix.size() && vecPrefix[i]<iTmpMulti) {
          iSum += (itcnt[Dbase_Ctrl_Blk::FreqIdx[vecPrefix[i]]]);
          ++i;
        }
        if(iSum != 0) dMean *= (double(iSum)/DBASE_NUM_TRANS);
//        dMean *= (double(itcnt[Dbase_Ctrl_Blk::FreqIdx[vecPrefix[i]]])/DBASE_NUM_TRANS);
      }
    } else dMean = pow(double(iDiv),double(-iDim));
    break;

    default:
    cerr << " INVALID dist_type = " << dist_type << endl;
    exit(0);
  }
  double dSupport=dMean+sqrt(log(dR)/(2*DBASE_NUM_TRANS));
  if(prune_type==CLIQUE) dSupport=MINSUP_PER;
  if(prune_type==MAFIA) dSupport = log(dR)*dMean/log(DBASE_NUM_TRANS*1.0);
  if(iDim < dDim) {
    switch(prune_type) {
      case PARABOLIC:
        dSupport=min(MINSUP_PER,dMean+(c-a*iDim*iDim)*sqrt(log(dR)/(2*DBASE_NUM_TRANS)));
        break;
      case EXPONENTIAL:
        dSupport=min(MINSUP_PER,dMean+sqrt(log(dR)*iDim2/(2*DBASE_NUM_TRANS*iDim)));
 //       dSupport=min(MINSUP_PER,dSupport);
        break;
      case CLIQUE:
        dSupport = MINSUP_PER;
        break;
      case MAFIA:
        dSupport = log(dR)*dMean/log(DBASE_NUM_TRANS*1.0);
        break;
    }
  }
  return (dSup >= dSupport);
}

double probOverlapRandom(vector<int> const &vec1, vector<int> const &vec2,int iDiv) {
  double dProb = 1.0;
  unsigned int i1=0;
  for(unsigned int i2=0;i2<vec1.size();++i2) {
    if(vec2[i1]/iDiv > vec1[i2]/iDiv)
      continue;
    while(vec2[i1]/iDiv != vec1[i2]/iDiv && i1 < vec2.size())
      ++i1;
    if(i1 < vec2.size()) {
      dProb *= (1+abs(vec2[i1]%iDiv - vec1[i2]%iDiv))/double(iDiv);
    }
  }
  return dProb;
}

void parse_args(int argc, char **argv)
{
   extern char * optarg;
   int c;

   if (argc < 2){
      cout << "usage: eclat\n";
      cout << "\t-a<algtype> (<0=eclat>, 1=charm, 2=basicmax, 3=maxcharm)  (default=" << alg_type << ")\n";
      cout << "\t-b (binary input file?)                                   (default=" << Dbase_Ctrl_Blk::binary_input << ")\n";
      cout << "\t-c<closedtype> (0=cnone, 1=chash, 2=cmax)                 (default=" << closed_type << ")\n";
      cout << "\t-d<difftype> (0=nodiff, 1=diff, 2=diff2, 3=diffin)        (default=" << diff_type << ")\n";
      cout << "\t-D<assumed distribution> (0=uniform, 1=estimated)         (default=" << dist_type << ")\n";
      cout << "\t-e<avg dim fraction>                                      (default=" << dAvgCluDim << ")\n";
      cout << "\t-f<find ?> (0=cluster,1=centre)                           (default=" << bFindCentres << ")\n";
      cout << "\t-i<infile>\n";
      cout << "\t-k<#(eq classes retained)>                                (default=" << TOP_K_EQ << ")\n";
      cout << "\t-l (output + output_tidlist?)                             (default=" << output << ' ' << output_idlist << ")\n";
      cout << "\t-m<maxdifftype> (<0=nodiff>, 1=diff, 2=diff2, 3=diffin)   (default=" << max_diff_type << ")\n";
      cout << "\t-n<#(datasets)>                                           (default=" << iNoOfDatasets << ")\n";
      cout << "\t-o (print patterns?)                                      (default=" << output << ")\n";
      cout << "\t-p<prune type> (0=PARABOLIC, 1=EXPONENTIAL, 2=CONSTANT)   (default=" << prune_type << ")\n";
      cout << "\t-r<interestingness=> p=O(n^(-r))>                         (default=" << dR << ")\n";
      cout << "\t-R<overlap>                                               (default=" << dRho << ")\n";
      cout << "\t-s<support>                                               (default=" << MINSUP_PER << ")\n";
      cout << "\t-S<absolute support> (specify -s or -S)                   (default=" << MINSUPPORT << ")\n";
      cout << "\t-x<# intervals/dimension>                                 (default=" << iXi << ")\n";
      cout << "\t-y (no offset?)                                           (default=" << Dbase_Ctrl_Blk::nooffset << ")\n";
      cout << "\t-z<sorttype> (0=nosort, <1=incr>, 2=decr, 3=incr_noclass) (default=" << sort_type << ")\n";
      exit(-1);
   }
   else{
     while ((c=getopt(argc,argv,"a:bc:d:D:e:f:i:k:lm:n:o:p:r:R:s:S:x:yz:"))!=-1){
       switch(c){
       case 'a':
          alg_type = (alg_vals) atoi(optarg);
          break;
       case 'b':
          Dbase_Ctrl_Blk::binary_input = true;
          break;
       case 'c':
          closed_type = (closed_vals) atoi(optarg);
          break;
       case 'd':
          diff_type = (diff_vals) atoi(optarg);
          break;
       case 'D':
          dist_type = atoi(optarg);
          break;
       case 'e':
          dAvgCluDim = atof(optarg);
          break;
       case 'f':
          bFindCentres = (bool) atoi(optarg);
          break;
       case 'i': //input files
	 sprintf(infile,"%s",optarg);
	 break;
       case 'k': //added by Karlton
	 TOP_K_EQ = atoi(optarg);
	 break;
       case 'l': //print idlists along with freq subtrees
          output=true;
          output_idlist = true;
          break;
       case 'm':
          max_diff_type = (diff_vals) atoi(optarg);
          break;
       case 'n': //added by Karlton
	 iNoOfDatasets = atoi(optarg);
	 break;
       case 'o': //print freq subtrees
	 output = true;
	 break;
       case 'p': //print freq subtrees
	 prune_type = (prune_vals)atoi(optarg);
	 break;
       case 'r':
          dR = atof(optarg);
          break;
       case 'R':
          dRho = atof(optarg);
          break;
       case 's': //support value for L2
	 MINSUP_PER = atof(optarg);
	 break;
       case 'S': //absolute support
	 MINSUPPORT = atoi(optarg);
	 break;
       case 'x':
          iXi = atoi(optarg);
          break;
       case 'y':
          Dbase_Ctrl_Blk::nooffset = true;
          break;
       case 'z':
          sort_type = (sort_vals) atoi(optarg);
          break;
       }               
     }
   }
   //cout << "ENUM VALS " << alg_type << " " << diff_type << " " 
   //     << max_diff_type << " " << closed_type << " " << endl;
}


ostream & operator<<(ostream& fout, idlist &vec){
  for (unsigned int i=0; i < vec.size(); ++i)
     fout << vec[i] << " ";
  return fout;
}


void get_F1(vector<int>& itcnt,long& DBASE_NUM_TRANS)
{
  TimeTracker tt;
  double te;
  int i, j, it;

  tt.Start();

  DBASE_MAXITEM=0;
  DCB->Cidsum = 0;
  
   while(DCB->get_next_trans())   // for each transaction in DB
   {
      for (i=0; i < DCB->TransSz; ++i){  // for each item in the transaction
         it = DCB->TransAry[i];

         if (it >= DBASE_MAXITEM){
            itcnt.resize(it+1,0);
            DBASE_MAXITEM = it+1;
         }
         ++itcnt[it];               // update count for that item
      }
      
      if (DCB->MaxTransSz < DCB->TransSz) DCB->MaxTransSz = DCB->TransSz;     
      ++DBASE_NUM_TRANS;
      DCB->Cidsum += DCB->Cid; //used to initialize hashval for closed set mining
   }
   
   //set the value of MINSUPPORT
   if (MINSUPPORT == -1)
     MINSUPPORT = (int) (MINSUP_PER*DBASE_NUM_TRANS+0.5);
   dR = pow(double(DBASE_NUM_TRANS),dR);     // added by Karlton
   if(dR < 3) {
     cerr << "ERROR: THRESHOLD r SET TOO LOW " << endl;
     exit(0);
   }
   
   if (MINSUPPORT<1) MINSUPPORT=1;
/*   cout<<"DBASE_NUM_TRANS : "<< DBASE_NUM_TRANS << endl;
   cout<<"DBASE_MAXITEM : "<< DBASE_MAXITEM << endl;
   cout<<"MINSUPPORT : "<< MINSUPPORT << " (" << MINSUP_PER << ")" << endl;*/

   //count number of frequent 1-itemsets
   DCB->NumF1 = 0;
   bool bOneFreq[DBASE_MAXITEM];
   for (i=0; i < DBASE_MAXITEM; ++i) {
     bOneFreq[i]=false;
     if (pruneTest(itcnt,itcnt[i],1,DBASE_NUM_TRANS,NULL,iXi)) {
       ++DCB->NumF1;
       bOneFreq[i]=true;
     }
   }

   //construct forward and reverse mapping from items to freq items
   DCB->FreqIdx = new int [DCB->NumF1];
   DCB->FreqMap = new int [DBASE_MAXITEM];
   for (i=0,j=0; i < DBASE_MAXITEM; ++i) {
      if (bOneFreq[i]==true) {
//       if (output && alg_type == eclat)cout << i << " - " << itcnt[i] << endl;
         DCB->FreqIdx[j] = i;
         DCB->FreqMap[i] = j;
         ++j;
      }
      else DCB->FreqMap[i] = -1;
   }
   
   //cout<< "F1 - " << DCB->NumF1 << " " << DBASE_MAXITEM << endl;  

   DCB->alloc_ParentClass(itcnt); // create Eqnodes for each frequent 1-itemset
   
   te = tt.Stop();
   stats.add(DBASE_MAXITEM, DCB->NumF1, te);
}

list<Eqclass *> * get_F2(vector<int>& itcnt,long DBASE_NUM_TRANS=1000)
{
  int i,j;
  int it1, it2;
  double te;
  TimeTracker tt;
  tt.Start();

  list<Eqclass *> *F2list = new list<Eqclass *>;

  //itcnt2 is a matrix of pairs p, p.first is count, p.second is flag
  int **itcnt2 = new int*[DCB->NumF1];

  //unsigned int **itcnt2 = new unsigned int *[DCB->NumF1];
  for (i=0; i < DCB->NumF1; ++i){  // for each frequent 1-itemset
    itcnt2[i] = new int [DCB->NumF1];
    //cout << "alloc " << i << " " << itcnt2[i] << endl;
    for (j=0; j < DCB->NumF1; ++j){
      itcnt2[i][j] = 0;
    }
  }
    
   while(DCB->get_next_trans())  // somehow the pointer is set to start of file
   {
      DCB->get_valid_trans();  // remove infrequent items from transaction
      DCB->make_vertical();
      //DCB->print_trans();
      //count a pair only once per cid
      for (i=0; i < DCB->TransSz; ++i){
         it1 = DCB->TransAry[i];
         for (j=i+1; j < DCB->TransSz; ++j){
            it2 = DCB->TransAry[j];
            ++itcnt2[it1][it2];
         }
      }
   }                           

   //compute class size ie # frequent 2-itemsets each frequent item is a part of
   DCB->class_sz = new int[DCB->NumF1]; 
   DCB->F2sum = new int[DCB->NumF1]; // sum of cnts of freq 2-itemsets
   for (i=0; i < DCB->NumF1; ++i){
      DCB->class_sz[i] = 0;
      DCB->F2sum[i] = 0;
   }
   
   for (i=0; i < DCB->NumF1; ++i){      // check each pair of frequent items
      for (j=i+1; j < DCB->NumF1; ++j){  // for frequent 2-itemsets
         if (pruneTest(itcnt,itcnt2[i][j],2,DBASE_NUM_TRANS,NULL,iXi)){
//            cerr << "TEST " << DCB->FreqIdx[i] << " " << DCB->FreqIdx[j] << " - " << itcnt2[i][j] << endl;
            ++DCB->class_sz[i];
            ++DCB->class_sz[j];
            DCB->F2sum[i] += itcnt2[i][j];
            DCB->F2sum[j] += itcnt2[i][j];
         }
      }
   }
   
   DCB->sort_ParentClass();
   int F2cnt=0;

   // count frequent patterns and generate eqclass
   Eqclass *eq;
   int sup;
   for (i=0; i < DCB->NumF1; ++i) {  // for each frequent item
      eq = new Eqclass;              // create an Eqclass
      eq->prefix().push_back(i);
      //if (alg_type==charm) eq->closedsupport()=DCB->ParentClass[i]->support();
      it1 = DCB->ParentClass[i]->val;
      for (j=i+1; j < DCB->NumF1; ++j) {  // containing frequent items based on 
         it2 = DCB->ParentClass[j]->val;
         if (it1 < it2) sup = itcnt2[it1][it2]; // ordering w/ which strong 
         else sup = itcnt2[it2][it1];
         if (pruneTest(itcnt,sup,2,DBASE_NUM_TRANS,NULL,iXi)){ // corr.exists
            ++F2cnt;
	    eq->add_node(j);
//	    if (output && alg_type == eclat) cout << DCB->FreqIdx[it1] << " " << DCB->FreqIdx[it2] << " - " << sup << endl;
         }
      }   
      F2list->push_back(eq);
   }
   
   //remap FreqIdx, FreqMap and ParentClass vals in sorted order
   //cout << "FREQMAP :\n";
   for (i=0; i < DCB->NumF1; ++i){
      DCB->FreqMap[DCB->FreqIdx[DCB->ParentClass[i]->val]] = i;
   }

   for (i=0; i < DBASE_MAXITEM; ++i)
      if (DCB->FreqMap[i] != -1){
         DCB->FreqIdx[DCB->FreqMap[i]] = i;
      }


   //cout << "ORDER :";
   for (i=0; i < DCB->NumF1; ++i){
      DCB->ParentClass[i]->val = i;
      //cout << " " << DCB->FreqIdx[i];
   }
   //cout << endl;
   
   for (i=0; i < DCB->NumF1; ++i) {
      delete [] itcnt2[i];
   }
   delete [] itcnt2;
   delete [] DCB->class_sz;
   delete [] DCB->F2sum;
   
   //cout << "F2 - " << F2cnt << " " << DCB->NumF1 * DCB->NumF1 << endl;
   te = tt.Stop();
   stats.add(DCB->NumF1 * DCB->NumF1, F2cnt, te);

   return F2list;
}

//performs l1 intersect l2
subset_vals get_intersect(idlist *l1, idlist *l2, idlist *join, int &idsum, 
                          int minsup=0)
{
   int diffmax1, diffmax2;
   diffmax1 = l1->size() - minsup;
   diffmax2 = l2->size() - minsup;
   
   int diffcnt1 = 0, diffcnt2 = 0;
   int n1, n2;
   unsigned int i1 = 0, i2 = 0;

   idsum = 0;
   while (i1 < l1->size() && i2 < l2->size()){
	// && diffcnt1 <= diffmax1 && diffcnt2 <= diffmax2)
      n1 = (*l1)[i1];
      n2 = (*l2)[i2];

      //look for matching cids
      if (n1 < n2){
         ++i1;
         ++diffcnt1;
      }
      else if (n1 > n2){
         ++i2;
         ++diffcnt2;
      }
      else{
         join->push_back(n1);
         idsum += n1;
         ++i1;
         ++i2;
      }
   }

   if (i1 < l1->size()) ++diffcnt1;
   if (i2 < l2->size()) ++diffcnt2;
   
   if (diffcnt1 == 0 && diffcnt2 == 0) return equals;
   else if (diffcnt1 == 0 && diffcnt2 > 0) return subset;
   else if (diffcnt1 > 0 && diffcnt2 == 0) return superset;
   else return notequal;   
}

//performs l1 - l2
subset_vals get_diff (idlist *l1, idlist *l2, idlist *join, 
                      int &idsum, int diffmax=INT_MAX)
{
//     insert_iterator<idlist> differ(*join,join->begin());
//     set_difference(l1->begin(), l1->end(), 
//                      l2->begin(), l2->end(), 
//                      differ);
   int n1, n2;
   int diffcnt1 = 0, diffcnt2 = 0;
   unsigned int i1 = 0, i2 = 0;
   
   idsum = 0;

   //while (i1 < l1->size() && i2 < l2->size() && diffcnt1 <= diffmax)
    while (i1 < l1->size() && i2 < l2->size()){
      n1 = (*l1)[i1];
      n2 = (*l2)[i2];

      if (n1 < n2){
         //implies that n1 is not to be found in n2
         join->push_back(n1);
         ++diffcnt1;
         idsum += n1;
         
         ++i1;
      }
      else if (n1 > n2){
         ++i2;
         ++diffcnt2;
      }
      else{
         ++i1;
         ++i2;
      }
   }   
   
   //add any remaining elements in l1 to join
   while (i1 < l1->size()){
      join->push_back((*l1)[i1]);
      idsum += (*l1)[i1];
      ++i1;
      ++diffcnt1;
   }
   
   if (i2 < l2->size()) ++diffcnt2;

   if (diffcnt1 == 0 && diffcnt2 == 0) return equals;
   else if (diffcnt1 == 0 && diffcnt2 > 0) return superset;
   else if (diffcnt1 > 0 && diffcnt2 == 0) return subset;
   else return notequal;   
}

subset_vals get_join(Eqnode *l1, Eqnode *l2, Eqnode *join, int iter)
{
   int diffmax = l1->support()-MINSUPPORT;
   int idsum;

   subset_vals sval = notequal;
 
   //compute tidset or diffset for join of l1 nd l2
   switch (diff_type){
   case diff2:
      if (iter == 2) sval = get_diff(&l1->tidset, &l2->tidset, 
                                     &join->tidset, idsum, diffmax);
      else sval = get_diff(&l2->tidset, &l1->tidset, &join->tidset, 
                           idsum, diffmax);
      if (sval == subset) sval = superset;
      else if (sval == superset) sval = subset;
      join->support() = l1->support() - join->tidset.size();            
      join->hashval() = l1->hashval() - idsum;
      break;
   case nodiff:
      sval = get_intersect(&l1->tidset, &l2->tidset, 
                           &join->tidset, idsum, MINSUPPORT);
      join->support() = join->tidset.size();      
      join->hashval() = idsum;
      break;
   case diffin:
      sval = get_diff(&l2->tidset, &l1->tidset, &join->tidset, idsum, diffmax);
      if (sval == subset) sval = superset;
      else if (sval == superset) sval = subset;
      join->support() = l1->support() - join->tidset.size();            
      join->hashval() = l1->hashval() - idsum;
      break;
   case diff:
      if (iter == 2){
         sval = get_intersect(&l1->tidset, &l2->tidset,
                              &join->tidset, idsum, MINSUPPORT);
         join->support() = join->tidset.size();      
         join->hashval() = idsum;
      }
      else{
         if (iter == 3) 
            sval = get_diff(&l1->tidset, &l2->tidset, &join->tidset, 
                            idsum, diffmax);
         else
            sval = get_diff(&l2->tidset, &l1->tidset, &join->tidset, 
                            idsum, diffmax);
         if (sval == subset) sval = superset;
         else if (sval == superset) sval = subset;
         join->support() = l1->support() - join->tidset.size();            
         join->hashval() = l1->hashval() - idsum;
      }
      break;
   }

   ++Stats::numjoin;
   return sval;
}

void get_max_join(Eqnode *l1, Eqnode *l2, Eqnode *join, int iter)
{
   int idsum;
   //find local maximal context for join
   //i.e., which maximal sets contain join as a subset
   switch(max_diff_type){
   case nodiff:
      get_intersect(&l1->maxset, &l2->maxset, &join->maxset, idsum);
      join->maxsupport() = join->maxset.size();  
      break;
   case diff2:
      if (iter == 2) get_diff(&l1->maxset, &l2->maxset, &join->maxset, idsum);
      else get_diff(&l2->maxset, &l1->maxset, &join->maxset, idsum);
      join->maxsupport() = l1->maxsupport() - join->maxset.size();
      break;
   case diffin:
      get_diff(&l2->maxset, &l1->maxset, &join->maxset, idsum);
      join->maxsupport() = l1->maxsupport() - join->maxset.size();
      break;
   case diff:
      cout << "diff NOT HANDLED\n";
      exit(-1);
      break;
   }
}

void get_Fk(vector<int>& itcnt,list<Eqclass *> &F2list,vector<Eqclass>& eqBest,long DBASE_NUM_TRANS=1000){
   Eqclass *eq;
   idlist newmax;
   int iter = 2,iDiv=iXi;
   double dThreshold=.05*DBASE_MAXITEM/(iDiv*iDiv)+sqrt(.05*DBASE_MAXITEM*log(dR)/(2*iDiv));
/*   if(dThreshold < 3) {
     cerr << "ERROR: THRESHOLD e=" <<  dThreshold << " SET TOO LOW. ";
     cerr << "Min threshold=" << log(DBASE_MAXITEM) << endl;
     exit(0);
   }*/
   while(!F2list.empty()){  // while there remain equivalence classes
      eq = F2list.front();  // pop the first equivalence class

      switch(alg_type){
      case eclat:
         form_f2_lists(eq);
         enumerate_freq(itcnt,eq, iter+1);
         break;
      case charm:
         form_closed_f2_lists(eq);
         newmax.clear();
         enumerate_closed_freq(itcnt,eq, iter+1, newmax);
         break;
      case basicmax:
         form_f2_lists(eq);
         newmax.clear();
         enumerate_max_freq(itcnt,eqBest,eq, iter+1, newmax,DBASE_NUM_TRANS,dThreshold);
         break;
      case maxcharm:
         form_closed_f2_lists(eq);
         newmax.clear();
         enumerate_max_closed_freq(itcnt,eq, iter+1, newmax);
         break;
      default:
         cerr << " Invalid algorithm type " << alg_type << endl;
         exit(0);
      }      
      delete eq;
      F2list.pop_front();
   }
}

bool readIBMPt(ifstream& fData,vector<int>& vecPt) {
   string sLine;
   getline(fData,sLine);
   if(fData.eof()) {
//     cerr << "ERROR: CANNOT READ ANOTHER POINT FROM FILE" << endl;
     return false;
   }
   istringstream istrs((char *)sLine.c_str());
   int iCid,iTid,iTransSz;
   istrs >> iCid;istrs >> iTid;istrs >> iTransSz;
   for(int iDim=0; iDim < iTransSz; ++iDim) {
     istrs >> vecPt[iDim];
   }
   return true;
}

double sim(vector<int>& vecPrefix1,vector<int>& vecPrefix2,int iDiv,int iMode=1) {
  double dSim = 0.0;
  vector<int> vecResult(vecPrefix2.size());
  unsigned int j=0,k=0,l=0,m=0,i1Intervals=0,i2Intervals=0;
  int i=0,iTmpMulti=0;
  switch(iMode) {
    case 0:
    dSim=set_intersection(vecPrefix2.begin(),vecPrefix2.end(),vecPrefix1.begin(),vecPrefix1.end(),vecResult.begin())-vecResult.begin();
    break;

    case 1:
    for(i=0;i<DBASE_MAXITEM/iDiv;++i) { // compute Jaccard coefficient
      iTmpMulti += iDiv;
      i1Intervals = 0,i2Intervals=0;
      while(j < vecPrefix1.size() && vecPrefix1[j]<iTmpMulti) {
        ++i1Intervals;
        ++j;
      }
      while(k < vecPrefix2.size() && vecPrefix2[k]<iTmpMulti) {
        ++i2Intervals;
        ++k;
      }
      if(i1Intervals > 0) ++l;
      if(i2Intervals > 0) ++m;
      if(i1Intervals > 0 && i2Intervals > 0) {
        vector<int> vecNum(i2Intervals);
        vector<int> vecDen(i1Intervals + i2Intervals);
        dSim += double(set_intersection(vecPrefix2.begin()+k-i2Intervals,vecPrefix2.begin()+k,vecPrefix1.begin()+j-i1Intervals,vecPrefix1.begin()+j,vecNum.begin())-vecNum.begin())/(set_union(vecPrefix2.begin()+k-i2Intervals,vecPrefix2.begin()+k,vecPrefix1.begin()+j-i1Intervals,vecPrefix1.begin()+j,vecDen.begin())-vecDen.begin());
      }
      if(j>= vecPrefix1.size() || k >= vecPrefix2.size()) break;
    }
    dSim = dSim*DBASE_MAXITEM/(iDiv*min(l,m));
    break;

    default:
      cerr << " Invalid Similarity mode : " << iMode << endl;
      exit(0);
  }
  return dSim;
}

int* classifyDataset(ifstream& fIn,ifstream& fData, int const* iaConfMatrix,double* daProbMatrix,vector<Eqclass>& eqBest,int iDiv=10) {
  int iNoOfDims,iNoOfPts,iNumClusters=0,iCorrect=0;
  double dNoise;
  fIn.seekg(0,ios::beg);
  fIn >> iNoOfPts;
  fIn >> iNoOfDims;
  fIn >> iNumClusters; fIn >> dNoise;
  vector<int> vFracCluPts(iNumClusters);
  fIn >> dNoise;
  vFracCluPts[0]=int(rint(dNoise*iNoOfPts));
  for (int iClassId = 1; iClassId < iNumClusters; ++iClassId) {
    fIn >> dNoise;
    vFracCluPts[iClassId] = int(rint(dNoise*iNoOfPts))+vFracCluPts[iClassId-1];
  }
  int iCurrentClass=0,iPtId=0;
  vector<int> vecPt(iNoOfDims);   // for each planted subspace
  fData.seekg(0,ios::beg);
  while(readIBMPt(fData,vecPt)) {
    // read each point in it and find the closest subspace
    if(iPtId++ > vFracCluPts[iCurrentClass]) ++iCurrentClass;
    vector<double> vecPtd(eqBest.size());
    vector<double> dProd(iNumClusters);
    for(unsigned int iSubSpcId = 0;iSubSpcId < eqBest.size();++iSubSpcId) {
      double dSim = sim(vecPt,eqBest[iSubSpcId].prefix(),iDiv,1);
      vecPtd[iSubSpcId]=dSim;
    }
    double dMaxSim = 0.0;
    int iNearestClass = eqBest.size();
    for (int iClassId = 0; iClassId < iNumClusters; ++iClassId) {
      dProd[iClassId]=0.0;
      for(unsigned int iSubSpcId = 0;iSubSpcId < eqBest.size();++iSubSpcId)
        dProd[iClassId] += (vecPtd[iSubSpcId]*iaConfMatrix[iClassId*(eqBest.size()+1)+iSubSpcId]);
      dMaxSim = max(dMaxSim,dProd[iClassId]);
      if(dMaxSim == dProd[iClassId]) iNearestClass = iClassId;
    }
    if(iCurrentClass==iNearestClass)++iCorrect;
  }
  cerr << " Correct = " << iCorrect << " Accuracy = " << double(iCorrect)/vFracCluPts[iCurrentClass] << endl;
  return NULL;
}

int* assignPoints(ifstream& fIn,ifstream& fData,int& iNumClusters, double const* daProbMatrix,vector<Eqclass>& eqBest,vector<int>& vAssign,double& dNoise,int iDiv=10) {
   int iNoOfDims,iNoOfPts, iNoOfCluPts;
   double dFracCluPts;
   fIn.seekg(0,ios::beg);
   fIn >> iNoOfPts;
   fIn >> iNoOfDims;
   fIn >> iNumClusters;
   fIn >> dNoise;
   int* iaConfMatrix=new int[(iNumClusters+1)*(eqBest.size()+1)];
   vector<int> vecPt(iNoOfDims);   // for each planted subspace
   double dThreshold=dAvgCluDim*iNoOfDims/iDiv+sqrt(dAvgCluDim*iNoOfDims*log(dR)/2);
   if(dThreshold < 2) {
     cerr << "ERROR: THRESHOLDS e/r SET TOO LOW. ";
     cerr << "Min threshold=" << log(iNoOfDims) << endl;
     exit(0);
   }
   vector<vector<int> > vvAssign(eqBest.size());
#ifdef DNA
/*   char *caGeneNames="yeast_genes";
   ifstream fName(caGeneNames);
   if(fName.bad() || fName.eof()) {
     cerr << "ERROR: READING FILE " << caGeneNames << endl;
     exit(0);
   }
   vector<string> vecNames(iNoOfPts);
   for(int iPtId = 0;iPtId < iNoOfPts;++iPtId) 
     fName >> vecNames[iPtId];*/
#endif
   int iRealPtId=-1;
   bool baEmpty[eqBest.size()+1];
   for(unsigned int j=0;j < eqBest.size()+1;++j) baEmpty[j]=true;
   baEmpty[eqBest.size()+1]=false;
   fData.seekg(0,ios::beg);
   for (int iClassId = 0; iClassId < iNumClusters+1; ++iClassId) {
     for(unsigned int j=0;j < eqBest.size()+1;++j)
       iaConfMatrix[iClassId*(eqBest.size()+1)+j]=0;
     if(iClassId < iNumClusters) {
       if(!(fIn >> dFracCluPts)) {
         cerr << "ERROR: in config file" << endl;
         exit(0);
       }
       iNoOfCluPts = int(rint(dFracCluPts*iNoOfPts));
     } else {
       if(dNoise==0)break;
       iNoOfCluPts = int(rint(dNoise*iNoOfPts));
     }
     for(int iPtId = 0;iPtId < iNoOfCluPts;++iPtId) {
       // read each point in it and find the closest subspace
       if(readIBMPt(fData,vecPt)==false) break;
       double dMaxSig = 0.0,dMaxSim = 0.0;
       int iNearestSubspace = eqBest.size();
       ++iRealPtId;
       for(unsigned int iSubSpcId = 0;iSubSpcId < eqBest.size();++iSubSpcId) {
         double dSim = sim(vecPt,eqBest[iSubSpcId].prefix(),iDiv,1);
         double dSig = (dSim-eqBest[iSubSpcId].prefix().size()/iDiv)*(dSim-eqBest[iSubSpcId].prefix().size()/iDiv)/eqBest[iSubSpcId].prefix().size();
	 dMaxSig = max(dSig,dMaxSig);
	 dMaxSim = max(dSim,dMaxSim);
 //        if(dMaxSig == dSig) iNearestSubspace = iSubSpcId;
         if(dMaxSim == dSim) iNearestSubspace = iSubSpcId;
       }
 //      if(exp(-2*dMaxSig)*dR < 1.0) {
       if(dMaxSim > dThreshold) {
         iaConfMatrix[iClassId*(eqBest.size()+1)+iNearestSubspace] += 1;
         baEmpty[iNearestSubspace]=false;
         vvAssign[iNearestSubspace].push_back(iRealPtId);
         vAssign[iRealPtId]=iNearestSubspace;
       } else {
         iaConfMatrix[iClassId*(eqBest.size()+1)+eqBest.size()] += 1;
         vAssign[iRealPtId]=-1;
       }
     }
   }
   if(false) {
     for (int iClassId = 0; iClassId < iNumClusters; ++iClassId) {
       cerr << endl;
       for(unsigned int iSubSpcId = 0;iSubSpcId < eqBest.size()+1;++iSubSpcId) {
         if(baEmpty[iSubSpcId]==false || iSubSpcId == eqBest.size()) {
           cerr.width(4);cerr.fill(' ');
           cerr << iaConfMatrix[iClassId*(eqBest.size()+1)+iSubSpcId];
         }
       }
     }
     cerr << endl;
     cerr << " Size = " << eqBest.size() << endl;
   }
#ifdef DNA
   for(int i=0;i<eqBest.size();i++) {
     if(vvAssign[i].size()<3) continue;
     for(int j=0;j<vvAssign[i].size();j++) {
       cout << vvAssign[i][j] << ' ';
//       cout << vecNames[vvAssign[i][j]] << ' ';
     }
     cout << endl;
   }
#endif
   return iaConfMatrix;
}

double *computeEntropy(int* iaConfMatrix,int iNumClusters,vector<Eqclass>& eqBest,double tottime,double dNoise,long DBASE_NUM_TRANS=1000) {
  double *daProbMatrix = new double[iNumClusters*eqBest.size()];
  int iaNoOfSubPts[eqBest.size()+1];
  double dEntropy=0.0;
  for(unsigned int iSubSpcId = 0;iSubSpcId < eqBest.size()+1;++iSubSpcId) {
    iaNoOfSubPts[iSubSpcId]=0;
    for (int iClassId = 0; iClassId < iNumClusters; ++iClassId) {
      iaNoOfSubPts[iSubSpcId] += iaConfMatrix[iClassId*(eqBest.size()+1)+iSubSpcId];
    }
    if(iSubSpcId==eqBest.size()) continue;
    double dSum = 0.0;
    for (int iClassId = 0; iClassId < iNumClusters; ++iClassId) {
      if(iaNoOfSubPts[iSubSpcId]!=0) {  // compute p_{iClassId,iSubSpcId}
        daProbMatrix[iClassId*eqBest.size()+iSubSpcId] = double(iaConfMatrix[iClassId*(eqBest.size()+1)+iSubSpcId])/iaNoOfSubPts[iSubSpcId];
        if(daProbMatrix[iClassId*eqBest.size()+iSubSpcId]>0)
          dSum -= (daProbMatrix[iClassId*eqBest.size()+iSubSpcId]*log(daProbMatrix[iClassId*eqBest.size()+iSubSpcId]));
      } else daProbMatrix[iClassId*eqBest.size()+iSubSpcId]=0;
    }  // compute \sum_{iClassId} -p_{iClassId,iSubSpcId}ln(p_{iClassId,iSubSpcId})
    dEntropy += (iaNoOfSubPts[iSubSpcId]*dSum);
  }
  dEntropy /= DBASE_NUM_TRANS;
  cerr << tottime/60 << ' ' << eqBest.size() << ' ' << dEntropy << " " << (1.0-iaNoOfSubPts[eqBest.size()]/((1.0-dNoise)*DBASE_NUM_TRANS)) << endl;
  return daProbMatrix;
}

vector<double>* getSubspaceDistr(int iNumClusters,vector<int>& vAssign,ifstream& fData,int iDiv,long DBASE_NUM_TRANS) {
   vector<int> vecPt(DBASE_MAXITEM/iDiv);
   vector<double>* vecDistr = new vector<double>(iNumClusters*DBASE_MAXITEM);
   fData.seekg(0,ios::beg);
   for(unsigned int i=0;i<vAssign.size();++i) { // for each point in dataset
     if(vAssign[i] == -1) continue;    // which is classified
     if(readIBMPt(fData,vecPt)==false) {
       cerr << " ERROR reading point " << i << " from data file " << endl;
       exit(0);
     }
     for(unsigned int j=0;j<vecPt.size();++j) { // for each dimension
       (*vecDistr)[DBASE_MAXITEM*vAssign[i]+vecPt[j]]++;
     }
   }
   for(int i=0;i<iNumClusters;++i) {
     double dSum = accumulate(vecDistr->begin()+i*DBASE_MAXITEM,vecDistr->begin()+i*DBASE_MAXITEM+iDiv,0.0);
     if(dSum <= log(DBASE_NUM_TRANS)) {
       for(unsigned int k=0;k<vecPt.size()*iDiv;++k)
         ((*vecDistr)[i*DBASE_MAXITEM+k])=0.0;
       continue;
     }
     for(unsigned int j=0;j<vecPt.size();++j) {
       for(int k=0;k<iDiv;++k) 
         ((*vecDistr)[i*DBASE_MAXITEM+j*iDiv+k]) /= dSum;
     }
   }
   return vecDistr;
}

vector<double>* refineSubspaces(vector<Eqclass>& eqBest,vector<int>& vAssign,ifstream& fData,int iDiv,long DBASE_NUM_TRANS) {
   vector<double>* vecDistr=getSubspaceDistr(eqBest.size(),vAssign,fData,iDiv,DBASE_NUM_TRANS);
   for(unsigned int i=0;i<eqBest.size();++i) {
     vector<int> vecNew;
     double dSum = accumulate(vecDistr->begin()+i*DBASE_MAXITEM,vecDistr->begin()+i*DBASE_MAXITEM+iDiv,0.0);
     if(dSum == 0.0) {
       eqBest.erase(eqBest.begin()+i);
       (*vecDistr).erase(vecDistr->begin()+DBASE_MAXITEM*i,vecDistr->begin()+(i+1)*DBASE_MAXITEM);
       --i;
       continue;
     }
     for(int j=0;j<DBASE_MAXITEM/iDiv;++j) {
       for(int k=0;k<iDiv;++k) {
         if((*vecDistr)[i*DBASE_MAXITEM+j*iDiv+k] > (1+3*(1-1.0/iDiv))/iDiv)
           vecNew.push_back(j*iDiv+k);
       }
     }
     eqBest[i].set_prefix(vecNew);
   }
   return vecDistr;
}

void processArgTypes(ofstream& summary) {
   switch(alg_type){
     case basicmax: summary << "MAX "; break;
     case charm: summary << "CHARM ";
        switch(closed_type){
          case cnone: break;
          case chash: summary << "CHASH "; break;
          case cmax: summary << "CMAX "; break;
        }
      break;
     case maxcharm: summary << "MAXCHARM "; break;
     case eclat: summary << "ECLAT "; break;
     case colex: summary << "COLEX "; break;
   }
   switch(diff_type){
     case nodiff: summary << "NODIFF "; break;
     case diff: summary << "DIFF "; break;
     case diff2: summary << "DIFF2 "; break;
     case diffin: summary << "DIFFIN "; break;
   } 
   switch(max_diff_type){
     case nodiff: summary << "MNODIFF "; break;
     case diff: summary << "MDIFF "; break;
     case diff2: summary << "MDIFF2 "; break;
     case diffin: summary << "MDIFFIN "; break;
   }
   switch(sort_type){
     case nosort: summary << "NOSORT "; break;
     case incr: summary << "INCR "; break;
     case incr_noclass: summary << "INCRNOC "; break;
     case decr: summary << "DECR "; break;
   }
/*   switch(prune_type){
     case prune: summary << "PRUNE "; break;
     case noprune: break;
   }*/
}

// added by Karlton to merge similar subspaces in eqBest
void mergeSimilarSubspaces(int iDiv,vector<double>* vecSim,vector<Eqclass>& eqBest) {
   for(unsigned int i=0; i <(*vecSim).size(); ++i) {
     for(unsigned int j=i+1; j <(*vecSim).size(); ++j) {
       int iRowSize = (int)sqrt((*vecSim).size());
       if((*vecSim)[i*iRowSize+j] > .6) {
         eqBest.erase(eqBest.begin()+j);
         (*vecSim).erase((*vecSim).begin()+j*iRowSize,(*vecSim).begin()+(j+1)*iRowSize);
         for(unsigned int i=0; i < (*vecSim).size(); ++i)
           (*vecSim).erase((*vecSim).begin()+i*iRowSize+j);
       }
     }
   }
}

// added by Karlton to merge similar subspaces in eqBest
void mergeSimilarSubspaces(int iDiv,vector<Eqclass>& eqBest,int iMode=1) {
   sort_heap(eqBest.begin(),eqBest.end());
   for(unsigned int i=0; i <eqBest.size(); ++i) {
     vector<int>& vecPrefix1(eqBest[i].prefix());
     for(unsigned int j=i+1; j <eqBest.size(); ++j) {
       vector<int>& vecPrefix2(eqBest[j].prefix());
       double dSim = 0.0;
       switch(iMode) {
         case 0:
           dSim = sim(vecPrefix1,vecPrefix2,iDiv,0)/min(vecPrefix2.size(),vecPrefix1.size());
           break;

         case 1:
           dSim = sim(vecPrefix1,vecPrefix2,iDiv,1)*iDiv/DBASE_MAXITEM;
           break;

         default:
           cerr << "ERROR: Invalid sim mode: " << iMode << endl;
           exit(0);
       }
       if(dSim > dRho) {
         vector<int> vecResult(vecPrefix2.size()+vecPrefix1.size());
         vector<int>::iterator vIter=set_union(vecPrefix2.begin(),vecPrefix2.end(),vecPrefix1.begin(),vecPrefix1.end(),vecResult.begin());
         vecResult.resize(vIter-vecResult.begin());
         vecPrefix1 = vecResult;
         eqBest.erase(eqBest.begin()+j);
         --j;
       }
     }
   }
}

void printSubspaces(int iDiv,vector<Eqclass>& eqBest) {
   for (unsigned int i=0; i < eqBest.size(); ++i) {
     vector<int> vecPrefix=eqBest[i].prefix();
     int iTmpMulti = iDiv;
     for(unsigned int j=0;j<vecPrefix.size();++j) {
       cout << vecPrefix[j] << ' ';
     }
     while(iTmpMulti < DBASE_MAXITEM/iDiv) {
       iTmpMulti += iDiv;
     }
     cout << endl;
   }
}

vector<double>* convertSubspacesToSimMatrix(int iDiv,vector<Eqclass>& eqBest,vector<double>* vecDistr) {
  vector<double>* daSim2 = new vector<double>(eqBest.size()*eqBest.size());
  cout << eqBest.size() << endl;
  int iRowSize = (*vecDistr).size()/eqBest.size();
  for(unsigned int i=0; i < eqBest.size(); ++i) {
//    vector<int> vecPrefix1(eqBest[i].prefix());
    for(unsigned int j=i; j < eqBest.size(); ++j) {
//      vector<int> vecPrefix2(eqBest[j].prefix());
//      (*daSim2)[i*eqBest.size()+j]=(*daSim2)[j*eqBest.size()+i]=sim(vecPrefix1,vecPrefix2,iDiv,1);
      (*daSim2)[i*eqBest.size()+j]=0.0;
      for(int k=0; k<iRowSize; ++k)
        (*daSim2)[i*eqBest.size()+j]+=(*vecDistr)[i*iRowSize+k]*(*vecDistr)[j*iRowSize+k];
      (*daSim2)[j*eqBest.size()+i]=(*daSim2)[i*eqBest.size()+j];
    }
  }
  for(unsigned int i=0; i < eqBest.size(); ++i) {
    for(unsigned int j=0;j < eqBest.size();++j) {
      cout << (*daSim2)[i*eqBest.size()+j]/sqrt((*daSim2)[j*eqBest.size()+j]*(*daSim2)[i*eqBest.size()+i]) << ' ';
    }
    cout << endl;
  }
  return daSim2;
}

int main(int argc, char **argv) {
  TimeTracker tt;
  parse_args(argc, argv); 
  int iDiv = iXi;
  double dNoise = 0.05;
  long DBASE_NUM_TRANS=0;
  
#ifdef OLD_OUTPUT
  ofstream summary("summary.out", ios::app);
  processArgTypes(summary);
#endif

  for(int iFile=0;iFile < iNoOfDatasets; ++iFile) {
    tt.Start(); 
    char cafile[1000];
    if(iNoOfDatasets > 1)
      sprintf(cafile,"%.*ss%.2d%s",int(strlen(infile)-4),infile,iFile,infile+strlen(infile)-4);
    else strcpy(cafile,infile);
    vector<int> itcnt(10,0);
    DCB = new Dbase_Ctrl_Blk(cafile); 
    get_F1(itcnt,DBASE_NUM_TRANS);
    list<Eqclass *> *F2list = get_F2(itcnt,DBASE_NUM_TRANS);

    //DCB->print_vertical();
    vector<Eqclass> eqBest;
    get_Fk(itcnt,*F2list,eqBest,DBASE_NUM_TRANS);
 
    int iNumClusters=0;
    char caCluSizeFile[1000];
    memset(caCluSizeFile,'\0',1000);
    strncpy(caCluSizeFile,cafile,strlen(cafile)-4);
    ifstream fIn(caCluSizeFile);
    if(fIn.bad() || fIn.eof()) {
      cerr << "ERROR: READING FILE " << caCluSizeFile << endl;
      exit(0);
    }
    ifstream fData(cafile);
    if(fData.bad() || fData.eof()) {
      cerr << "ERROR: READING FILE " << cafile << endl;
      exit(0);
    }

    vector<int> vAssign(DBASE_NUM_TRANS);
    int *iaConfMatrix;
//    vector<double>* vecDistr;
//    mergeSimilarSubspaces(iDiv,vecSim,eqBest);
//  vector<double>* vecSim = convertSubspacesToSimMatrix(iDiv, eqBest,vecDistr);
    for(int i=0;i<3;++i) {
//      mergeSimilarSubspaces(iDiv,eqBest,1);
      delete[] assignPoints(fIn,fData,iNumClusters,(double *)NULL,eqBest,vAssign,dNoise,iDiv);
      delete refineSubspaces(eqBest,vAssign,fData,iDiv,DBASE_NUM_TRANS);
    }
    mergeSimilarSubspaces(iDiv,eqBest,1);
    iaConfMatrix=assignPoints(fIn,fData,iNumClusters,(double *)NULL,eqBest,vAssign,dNoise,iDiv);
    double tottime = tt.Stop();
#ifdef OLD_OUTPUT
    for (unsigned int i=0; i < stats.size(); ++i){
       cout << "F" << i+1 << " - " << stats[i].numlarge << " " << stats[i].numcand << " " << stats[i].nummax << endl;
    }
    stats.tottime = tottime;
    summary << cafile << " " << MINSUP_PER << " " << DBASE_NUM_TRANS << " " << MINSUPPORT << " " << stats << endl;
    summary.close();
    exit(0);
#endif
#ifndef DNA
    computeEntropy(iaConfMatrix,iNumClusters,eqBest,tottime,dNoise,DBASE_NUM_TRANS);
#endif
    if(output == true) printSubspaces(iDiv,eqBest);
//    else convertSubspacesToSimMatrix(iDiv, eqBest,vecDistr);
    fIn.close();
    fData.close();
    delete DCB;
  }
}
