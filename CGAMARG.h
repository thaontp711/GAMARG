#ifndef  __CGAMARG__
#define  __CGAMARG__

#include <stdio.h>
#include <stdlib.h>



#include <string>
#include <string.h>
#include <vector>

#include <queue>
#include <pthread.h>





#include "COneNode.h"
#include "Match.h"
#include "define.h"
#include "CData.h"
#include "utilities.h"


struct onepair
{
  COneNode* node1, *node2;
  int match;
  int mismatchpoint; //mismatch point
  bool matchLeft;
  char node1chardif;
  bool gamete;
  
};

struct oneoption
{
	int id;
    int idvecgs;
	int brpointgs;
    int matchlength;
};

struct gemsite
{
	int gamete; // is 0 for 00; 1 for 01; 2 for 10; 3 for 11
	int sleft;  // siteid has gamete in form [sleft,sright]
    int sright;
};


/*struct gamete
{
	int count00;
	int count01;
	int count10;
	int count11;
};*/

//typedef std::vector <sitepair> vec_sitepair;
typedef std::priority_queue <onepair> prioqueue_onepair;
typedef std::vector <onepair> vec_onepair;
typedef std::priority_queue <gemsite> prioqueue_gemsite;
//typedef std::vector <gemsite> vecgem;
typedef std::vector <oneoption> vec_option;

bool operator < (const onepair &x, const onepair &y);
bool operator < (const gemsite &x, const gemsite &y);

class CGAMARG
{
 private:
  //vecpCOneNode nodes;

  const CData *data;

  vecint activenodes;
  vecint mutlist; //list of possible mutation
  prioqueue_onepair *priqueue;
  prioqueue_gemsite prigemlist;
  vec_onepair coallist;
  vec_option vecop;
    //vecgem gemlist;
  //vec_sitepair sitelist; /////////////new added
  vecint count0, count1;
    vecint count00, count01, count10, count11;
    vecbool sitemarked;
  vecpCOneNode nodetothread;
  bool mu;

  //int countvec;
  int nthread;
    int distgam;
  int nRec;
  //  int nsite, nsp;
  void Add(COneNode* &node);
  void AddNodetoThread();

  void RemoveActiveID(int id);

  void Mutation(int siteid, int id);
  void Coal(COneNode *node1, COneNode* node2);

  bool PossibleMutation(void);
  bool PossibleCoal(void);
  bool PossibleGem(void); // check gamete
  bool Recom(void);
  void CheckGamete(char p, char p1, int idleft, int idright);
  void CheckGameteCoal(char p, char p1, int idleft, int idright);
  void CheckGameteMut1st(char p, char p1, int idleft, int idright);
  void CheckGameteMut2nd(char p, char p1, int idleft, int idright);
    bool IsNodegamete(COneNode *node, gemsite gs);

  void Recom(COneNode *node1, COneNode *node2, int siteid, bool matchLeft);
  void RecomGemnode(COneNode *node,  int matchpoint, bool matchLeft);
  void  CheckCount(void) const;

  void FindPath(int id) const;

  void Check(const COneNode *root, int siteid) const;

 public:

  vecpCOneNode nodes;
  CGAMARG(int nth, int ndist);
  ~CGAMARG();

  void ReadData(const CData *idata);
  void BuildARG(void);
  void View(void) const;

  void ReadrmSNPs      (char *fn);
  void ReadselectedSNPs(char *fn);
  void Write2ARG(FILE *fout);
  
  void Check(void) const;
};


#endif
