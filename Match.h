#ifndef __MATCH__
#define __MATCH__

#include <string.h>

#include "COneNode.h"
#include "define.h"


struct onethread
{
  vecpCOneNode *nodes;
  COneNode     *onenode;
  vecint activenodes, match, mismatchpoint;
  vecbool isLeft; // get matching direction Left or Right
};

void *LeftMatch(void*ptr);
int LeftMatch4Nodes(const COneNode *node1, const COneNode* node2, int &matchLeft, int &mismatchpointLeft);
int RightMatch4Nodes(const COneNode *node1, const COneNode* node2, int &matchRight, int &mismatchpointRight);

#endif
