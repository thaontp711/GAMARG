#ifndef __CARG__
#define __CARG__

#include <stdio.h>
#include <stdlib.h>



#include <string>
#include <string.h>
#include <vector>
#include "define.h"
#include "COneNode.h"
#include "CData.h"
#include "utilities.h"


struct oneargnode
{
  int id;
  int st, ed;
  int child1, child2;
  int parent1, parent2;
};

class CARG
{

 private:
  //vecstring names;
  int numsequences, nummarkers;
  
 public:
  vecpCOneNode nodes;
  CARG(const CData *idata);

  void Load(vecpCOneNode *nodes, const vecint &poss);
  void getMarginalTrees(FILE *ftrees,int noArg);
  void getMarginalTreeForMarker(FILE *ftrees,COneNode *root, int siteid, int noArg);

};


#endif
