#include<iostream>
#include <cmath>
#include <vector>
#include "eclat.h"
#include "timetrack.h"
#include "calcdb.h"
#include "eqclass.h"
#include "stats.h"
#include "maximal.h"
#include "chashtable.h"

typedef vector<bool> bit_vector;

//extern vars
extern Dbase_Ctrl_Blk *DCB;
extern bool bFindCentres;         // added by Karlton
extern unsigned int TOP_K_EQ;     // added by Karlton
extern int iXi;                   // added by Karlton
extern int dist_type;             // added by Karlton
extern double dRho;               // added by Karlton
extern Stats stats;
MaximalTest maxtest; //for maximality & closed (with cmax) testing
cHashTable hashtest; //for closed (with chash) testing
extern unsigned long iCnt;
//extern functions
extern subset_vals get_intersect(idlist *l1, idlist *l2, 
                                 idlist *join, int minsup=0);
extern subset_vals get_diff (idlist *l1, idlist *l2, idlist *join);
extern subset_vals get_join(Eqnode *l1, Eqnode *l2, Eqnode *join, int iter);
extern void get_max_join(Eqnode *l1, Eqnode *l2, Eqnode *join, int iter);
extern bool pruneTest(vector<int>& itcnt,int iSup,int iDim,long DBASE_NUM_TRANS=1000,Eqclass* neq=NULL,int iDiv=iXi);
extern double probOverlapRandom(vector<int> const&,vector<int> const&,int iDiv=iXi);
extern double sim(vector<int>& vecPrefix1,vector<int>& vecPrefix2,int iDiv,int iMode);

void print_tabs(int depth) {
   for (int i=0; i < depth; ++i)
      cout << "\t";
}

static bool notfrequent (vector<int>& itcnt,Eqnode &n,Eqclass *neq, int len,long DBASE_NUM_TRANS=1000){
  //cout << "IN FREQ " << n.sup << endl;
  return !pruneTest(itcnt,n.support(),len,DBASE_NUM_TRANS,neq,iXi);
}

bool updatePrefix(vector<int>& vecPrefix1,vector<int>& vecPrefix2,int iXi) {
  int iMinSize = min(vecPrefix2.size(),vecPrefix1.size());
  if(sim(vecPrefix1,vecPrefix2,iXi,0)/iMinSize > dRho) {
/*    Eqnode* join = new Eqnode ((eq)->val);
    get_join(*ni, *nj, join, iter);
    stats.incrcand(iter-1);
    if (notfrequent(itcnt,*join,neq,iter)) {
      delete join;
      return true;
    } else {
      get_max_join(*ni, *nj, join, iter);
      neq->add_node(join);
      stats.incrlarge(iter-1);
    }*/

    vector<int> vecResult(vecPrefix2.size()+vecPrefix1.size());
    vector<int>::iterator vIter=set_union(vecPrefix2.begin(),vecPrefix2.end(),vecPrefix1.begin(),vecPrefix1.end(),vecResult.begin());
    vecResult.resize(vIter-vecResult.begin());
    vecPrefix1 = vecResult;
    return false;
  }
  return true;
}

bool greedyJoin(Eqclass& eq1,Eqclass& eq2,int iXi,int iter,vector<int>& itcnt) {
  vector<int>& vecPrefix1 = eq1.prefix();
  vector<int>& vecPrefix2 = eq2.prefix();
  return updatePrefix(vecPrefix1,vecPrefix2,iXi);
}

void enumerate_max_freq(vector<int>& itcnt,vector<Eqclass>& eqBest,Eqclass *eq, int iter, idlist &newmax,long DBASE_NUM_TRANS,double dThreshold)
{
   int nmaxpos;
   TimeTracker tt;
   Eqclass *neq;
   list<Eqnode *>::iterator ni, nj;
   Eqnode *join;
   bool extend = false;
   bool subset = false;
   
   eq->sort_nodes();
   nmaxpos = newmax.size(); //initial newmax pos

   //if ni is a subset of maxset then all elements after ni must be
   for (ni = eq->nlist().begin(); ni != eq->nlist().end() && !subset; ++ni){
      tt.Start();

      neq = new Eqclass;
      neq->set_prefix(eq->prefix(),*(*ni));

      subset = false;
      extend = false;
      bool bAdd = true;       // should we add this Eqclass or merge it?
      bool res = maxtest.update_maxset(ni, eq, newmax, nmaxpos, subset);
      if (!subset || res){
            subset = maxtest.subset_test(ni, eq);
      }
      nmaxpos = newmax.size();
      if (!subset){
         nj = ni; ++nj;
         for (; nj != eq->nlist().end(); ++nj){ 
            join = new Eqnode ((*nj)->val);
            
            get_join(*ni, *nj, join, iter);
            //cout << "ISECT " << *join;
            stats.incrcand(iter-1);
            if (notfrequent(itcnt,*join,neq,iter)) delete join;
            else{
               extend = true;
               get_max_join(*ni, *nj, join, iter);
               neq->add_node(join);
               stats.incrlarge(iter-1);
            }
         }

         stats.incrtime(iter-1, tt.Stop());
         bAdd = true;
         if (!neq->nlist().empty()){
           if(bFindCentres) {
             Eqclass tmp;
             tmp.set_prefix(neq->prefix());
             tmp.sort_prefix();
             vector<int> vecPrefix = tmp.prefix();
             for(unsigned int i=0;i<eqBest.size();++i) {
               bAdd = updatePrefix(eqBest[i].prefix(),vecPrefix,iXi);
 //              bAdd = greedyJoin(eqBest[i],*neq,iXi,iter,itcnt);
               if(bAdd == false)break;
             }
           }
           if(bAdd) {
             enumerate_max_freq(itcnt,eqBest,neq, iter+1,newmax,DBASE_NUM_TRANS,dThreshold);
           }
         }
      }
      if (!extend && (*ni)->maxsupport() == 0){
         maxtest.add(eq->prefix(), (*ni)->val, (*ni)->support());
         newmax.push_back(maxtest.maxcount-1);
         stats.incrmax((iter-2));

         // added by Karlton
         double dSup = double((*ni)->support())/DBASE_NUM_TRANS;
         double dMean = 1.0;
         switch(dist_type) {
           case 0:
           dMean = pow(double(iXi),double(-iter));
           break;

           case 1:
           if(neq != (Eqclass *)NULL) {
             vector<int> vecPrefix=neq->prefix();
             for (unsigned int i=0; i < vecPrefix.size(); ++i){
               int iSum = 0;
               int iTmpMulti=(vecPrefix[i]/iXi + 1)*iXi;
               while(i < vecPrefix.size() && vecPrefix[i]<iTmpMulti) {
                 iSum += (itcnt[Dbase_Ctrl_Blk::FreqIdx[vecPrefix[i]]]);
                 ++i;
               }
               if(iSum != 0)dMean *= (double(iSum)/DBASE_NUM_TRANS);
             }
           } else dMean = pow(double(iXi),double(-iter));
           break;

           default:
           cerr << " Invalid dist_type = " << dist_type << endl;
           exit(0);
         }
         neq->set_rarity(-2*DBASE_NUM_TRANS*(dSup-dMean)*(dSup-dMean));
         neq->set_support((*ni)->support()*100.0/DBASE_NUM_TRANS);
         bAdd = true; // should we add neq or merge it w/ an existing one?
//         if(bFindCentres) {
           neq->sort_prefix();
           for(unsigned int i=0;i<eqBest.size();++i) {
             bAdd = updatePrefix(eqBest[i].prefix(),neq->prefix(),iXi);
             if(bAdd == false)break;
           }
//         }
         if(bAdd && neq->prefix().size() >= dThreshold) {
           if(eqBest.size() < TOP_K_EQ) {
             eqBest.push_back(*neq);
             if(eqBest.size()==TOP_K_EQ) make_heap(eqBest.begin(),eqBest.end());
           } else if(*neq < eqBest[0]) {
             pop_heap(eqBest.begin(),eqBest.end());
             eqBest.resize(TOP_K_EQ-1);
             eqBest.push_back(*neq);
             push_heap(eqBest.begin(),eqBest.end());
           }
         }
      }
      delete neq;
   }
}

void enumerate_max_closed_freq(vector<int>& itcnt,Eqclass *eq, int iter, idlist &newmax)
{
   int nmaxpos;

   TimeTracker tt;
   Eqclass *neq;
   list<Eqnode *>::iterator ni, nj;
   Eqnode *join;
   subset_vals sval;

   bool extend = false;
   bool subsetflg = false;
   
   eq->sort_nodes();

   //cout << "CMAX " << *eq;
   
   nmaxpos = newmax.size(); //initial newmax pos

   //if ni is a subset of maxset then all elements after ni must be
   for (ni = eq->nlist().begin(); ni != eq->nlist().end() && !subsetflg; ++ni){
      tt.Start();

      neq = new Eqclass;
      neq->set_prefix(eq->prefix(),*(*ni));

      subsetflg = false;
      extend = false;
      bool res = maxtest.update_maxset(ni, eq, newmax, nmaxpos, subsetflg);
      if (!subsetflg || res){
         subsetflg = maxtest.subset_test(ni, eq);
      }
      
      nmaxpos = newmax.size();

      if (!subsetflg){
         nj = ni; ++nj;
         for (; nj != eq->nlist().end();){ 
            join = new Eqnode ((*nj)->val);
            
            sval = get_join(*ni, *nj, join, iter);
            //cout << "SVAL " << (int)sval;
            //eq->print_node(*(*nj));
            //cout << endl;
            //cout << "ISECT " << *join;
            stats.incrcand(iter-1);
            if (notfrequent(itcnt,*join,neq,iter)){
               delete join;
               ++nj;
            }
            else{
               extend = true;
               //get_max_join(*ni, *nj, join, iter);
               //neq->add_node(join);
               stats.incrlarge(iter-1);

               switch(sval){
               case subset:
                  //add nj to all elements in eq by adding nj to prefix
                  neq->prefix().push_back((*nj)->val);               
                  //neq->closedsupport() = join->support();
                  delete join;
                  ++nj;
                  break;
               case notequal:
                  get_max_join(*ni, *nj, join, iter);
                  neq->add_node(join);               
                  ++nj;
                  break;
               case equals:
                  //add nj to all elements in eq by adding nj to prefix
                  neq->prefix().push_back((*nj)->val); 
                  //neq->closedsupport() = join->support();
                  delete *nj;
                  nj = eq->nlist().erase(nj); //remove nj
                  delete join;
                  break;
               case superset:
                  get_max_join(*ni, *nj, join, iter);
                  delete *nj;
                  nj = eq->nlist().erase(nj); //remove nj
                  //++nj;
                  neq->add_node(join);            
                  break;
               }
            }
         }
      
         stats.incrtime(iter-1, tt.Stop());
         if (neq->nlist().size() > 1){
            enumerate_max_closed_freq(itcnt,neq, iter+1, newmax);
         }
         else if (neq->nlist().size() == 1){
            nj = neq->nlist().begin();
            if ((*nj)->maxsupport() == 0){
               if (false) //output) 
                  neq->print_node(*(*nj)) << endl;
               maxtest.add(neq->prefix(), (*nj)->val);
               newmax.push_back(maxtest.maxcount-1);
               stats.incrmax(neq->prefix().size());
            }
         }
         else if (extend && (*ni)->maxsupport() == 0){
            if (false) //output)
               neq->print_prefix() << endl;
            maxtest.add(neq->prefix(), -1);
            newmax.push_back(maxtest.maxcount-1);
            stats.incrmax(neq->prefix().size()-1);
         }
      }
      
      if (!extend && (*ni)->maxsupport() == 0){
         if (false) //output)
            neq->print_prefix() << endl;
         maxtest.add(eq->prefix(), (*ni)->val);
         newmax.push_back(maxtest.maxcount-1);
         stats.incrmax(neq->prefix().size()-1);
      }

      delete neq;
   }
}

void enumerate_closed_freq(vector<int>& itcnt,Eqclass *eq, int iter, idlist &newmax)
{
   TimeTracker tt;
   Eqclass *neq;
   int nmaxpos;
   bool cflg;
   list<Eqnode *>::iterator ni, nj;
   Eqnode *join;
   subset_vals sval;

   nmaxpos = newmax.size(); //initial newmax pos
   eq->sort_nodes();
   //print_tabs(iter-3);
   //cout << "F" << iter << " " << *eq;
   for (ni = eq->nlist().begin(); ni != eq->nlist().end(); ++ni){
      neq = new Eqclass;
      neq->set_prefix(eq->prefix(),*(*ni));

      //cout << "prefix " << neq->print_prefix() << endl;
      tt.Start();

      if (closed_type == cmax) 
         maxtest.update_maxset(ni, eq, newmax, nmaxpos, cflg);

      nmaxpos = newmax.size(); //initial newmax pos
      nj = ni;
      for (++nj; nj != eq->nlist().end(); ){
         join = new Eqnode ((*nj)->val);
         sval = get_join(*ni, *nj, join, iter);
         stats.incrcand(iter-1);
         if (notfrequent(itcnt,*join,neq, iter)){
            delete join;
            ++nj;
         }
         else{
            stats.incrlarge(iter-1);
            switch(sval){
            case subset:
               //add nj to all elements in eq by adding nj to prefix
               neq->prefix().push_back((*nj)->val);               
               //neq->closedsupport() = join->support();
               //cout << "SUSET " << *join << endl;
               delete join;
               ++nj;
               break;
            case notequal:
               if (closed_type == cmax) get_max_join(*ni, *nj, join, iter);
               neq->add_node(join);               
               ++nj;
               break;
            case equals:
               //add nj to all elements in eq by adding nj to prefix
               neq->prefix().push_back((*nj)->val); 
               //neq->closedsupport() = join->support();
               delete *nj;
               nj = eq->nlist().erase(nj); //remove nj
               //cout << "EQUAL " << *join << endl;
               delete join;
               break;
            case superset:
               if (closed_type == cmax) get_max_join(*ni, *nj, join, iter);
               delete *nj;
               nj = eq->nlist().erase(nj); //remove nj
               //++nj;
               neq->add_node(join);            
               break;
            }
         }
      }
      
      cflg = true;
      if (closed_type == cmax){
         cflg = maxtest.check_closed(*ni);
         if (cflg){
            maxtest.add(neq->prefix(), -1, (*ni)->support());
            newmax.push_back(maxtest.maxcount-1);
         }
      }
      else if (closed_type == chash){
         cflg = hashtest.add(neq->prefix(), -1, (*ni)->support(), 
                             (*ni)->hval);
      }

      if (cflg){
         stats.incrmax(neq->prefix().size()-1);
         if (false) {//output){
            cout << "CLOSEDxy ";
            neq->print_prefix(true);
            cout << endl;
         }
      }
      
      stats.incrtime(iter-1, tt.Stop());
      //if (output) cout << *neq;
      if (neq->nlist().size() > 1){
         enumerate_closed_freq(itcnt,neq, iter+1, newmax);
      }
      else if (neq->nlist().size() == 1){
         cflg = true;
         if (closed_type == cmax){
            cflg = maxtest.check_closed(neq->nlist().front());
            if (cflg){
               maxtest.add(neq->prefix(), neq->nlist().front()->val, 
                           neq->nlist().front()->sup);
               newmax.push_back(maxtest.maxcount-1);
            }
         }
         else if (closed_type == chash){
            cflg = hashtest.add(neq->prefix(), neq->nlist().front()->val, 
                                neq->nlist().front()->sup, 
                                neq->nlist().front()->hval);
         }
         
         if (cflg){
            if (false) /*output*/ cout << "CLOSEDy " << *neq;
            stats.incrmax(neq->prefix().size());
         }
      }
      delete neq;
   }
}

void enumerate_freq(vector<int>& itcnt,Eqclass *eq, int iter)
{
   TimeTracker tt;
   Eqclass *neq;
   list<Eqnode *>::iterator ni, nj;
   Eqnode *join;

   for (ni = eq->nlist().begin(); ni != eq->nlist().end(); ++ni){
      neq = new Eqclass;
      neq->set_prefix(eq->prefix(),*(*ni));
      tt.Start();
      nj = ni;
      for (++nj; nj != eq->nlist().end(); ++nj){ 
         join = new Eqnode ((*nj)->val);
         get_join(*ni, *nj, join, iter);
         stats.incrcand(iter-1);
         if (notfrequent(itcnt,*join,neq, iter)) delete join;
         else{
            neq->add_node(join);
            stats.incrlarge(iter-1);
         }
      }
      stats.incrtime(iter-1, tt.Stop());
//      if (output) cout << *neq;
      if (neq->nlist().size()> 1){
         enumerate_freq(itcnt,neq, iter+1);
      }
      delete neq;
   }
}

void form_closed_f2_lists(Eqclass *eq)
{
   static bit_vector *bvec = NULL;
   if (bvec == NULL){
      bvec = new bit_vector(DCB->NumF1, true);
   }
   
   subset_vals sval;
   list<Eqnode *>::iterator ni;
   Eqnode *l1, *l2;
   int pit, nit;
   TimeTracker tt;
   bool extend = false;
   bool cflg;
   
   tt.Start();
   pit = eq->prefix()[0];
   l1 = DCB->ParentClass[pit];

   if (!(*bvec)[pit]){
      eq->nlist().clear();
     return;
   }
   
   if (alg_type == maxcharm){
      cflg = maxtest.subset_test_f2(eq);
      if (cflg){
         eq->nlist().clear();
         return;
      }
   }
   
   for (ni=eq->nlist().begin(); ni != eq->nlist().end(); ){
      nit = (*ni)->val;
      if (!(*bvec)[nit]){
         ++ni;
         continue;
      }
      l2 = DCB->ParentClass[nit];
      sval = get_join(l1, l2, (*ni), 2);
      //cout << "SVAL " << (int)sval << endl;
      switch(sval){
      case subset:
         //add nj to all elements in eq by adding nj to prefix
         eq->prefix().push_back((*ni)->val);               
         extend = true;
         //eq->closedsupport() = (*ni)->support();
         delete *ni;
         ni = eq->nlist().erase(ni); //remove ni
         //cout << "CAME HERE " << eq->nlist().size() << endl;
         break;
      case notequal:
         if (alg_type == maxcharm || closed_type == cmax) 
            get_max_join(l1, l2, (*ni), 2);
         ++ni;
         break;
      case equals:
         //add nj to all elements in eq by adding nj to prefix
         eq->prefix().push_back((*ni)->val); 
         extend = true;
         //eq->closedsupport() = (*ni)->support();
         delete *ni;
         ni = eq->nlist().erase(ni); //remove ni
         (*bvec)[nit] = false;
         //cout << "CAME HERE " << eq->nlist().size() << endl;
         break;
      case superset:
         (*bvec)[nit] = false;
         if (alg_type == maxcharm || closed_type == cmax) 
            get_max_join(l1, l2, (*ni), 2);
         ++ni;
         break;
      }
   }
   
   if (alg_type == charm){

      cflg = true;
      if (closed_type == cmax){
         cflg = maxtest.check_closed(l1);
         if (cflg) maxtest.add(eq->prefix(), -1, l1->support());
      }
      else if (closed_type == chash){
         cflg = hashtest.add(eq->prefix(), -1, l1->support(), 
                             l1->hval);
      }

      if (cflg){
         stats.incrmax(eq->prefix().size()-1);
         if (false) {
            cout << "CLOSEDzx ";
            eq->print_prefix(true);
            cout << endl;
         }
      }
      
      
      if (eq->nlist().size() == 1){
         cflg = true;
         if (closed_type == cmax){
            cflg = maxtest.check_closed(eq->nlist().front());
            if (cflg){
               maxtest.add(eq->prefix(), eq->nlist().front()->val, 
                           eq->nlist().front()->sup);
            }
         }
         else if (closed_type == chash){
            cflg = hashtest.add(eq->prefix(), eq->nlist().front()->val, 
                                eq->nlist().front()->sup, 
                                eq->nlist().front()->hval);
         }

         if (cflg){
            if (false) {
               cout << "CLOSEDz ";
               cout << *eq;
            }
            
            stats.incrmax(eq->prefix().size());
         }
         eq->nlist().clear();
      }
   }
   else if (alg_type == maxcharm){
      if (eq->nlist().empty() && l1->maxsupport() == 0){
         if (false)
            eq->print_prefix() << endl;
         maxtest.add(eq->prefix(), -1);
         stats.incrmax(eq->prefix().size()-1);
      }
      else if (eq->nlist().size() == 1){
         l1 = eq->nlist().front();
         if (l1->maxsupport() == 0){
            if (false) //output)
               eq->print_node(*l1) << endl;
            maxtest.add(eq->prefix(), l1->val);
            stats.incrmax(eq->prefix().size()); 
         }
         eq->nlist().clear();
      }
   }

   //cout << "F2 " << *eq << endl;
   stats.incrtime(1,tt.Stop());   
   //delete l1;
}

void form_f2_lists(Eqclass *eq)
{
   list<Eqnode *>::iterator ni;
   Eqnode *l1, *l2;
   int pit, nit;
   TimeTracker tt;
   
   tt.Start();
   pit = eq->prefix()[0];
   l1 = DCB->ParentClass[pit];   // for every item w/ which extension is viable
   for (ni=eq->nlist().begin(); ni != eq->nlist().end(); ++ni){
      nit = (*ni)->val;
      l2 = DCB->ParentClass[nit];

      get_join(l1, l2, (*ni), 2);
      if (alg_type == basicmax) 
         get_max_join(l1, l2, (*ni), 2);
   }
   //delete l1;
   stats.incrtime(1,tt.Stop());   
}
