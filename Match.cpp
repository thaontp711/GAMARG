#include "Match.h"


void *LeftMatch(void*ptr)
{
  onethread *e = (onethread*) ptr;
  vecpCOneNode *nodes = e->nodes;
  //int match, mismatchpoint;
  int matchLeft, matchRight, mismatchpointLeft, mismatchpointRight;
  for (int i = 0; i != (int) e->activenodes.size(); ++i)
    {
      int id = e->activenodes[i];
      LeftMatch4Nodes(e->onenode, nodes->at(id), matchLeft, mismatchpointLeft);
      RightMatch4Nodes(e->onenode, nodes->at(id), matchRight, mismatchpointRight);
        int j = rand() % 2;
      if((matchLeft > matchRight) || ((matchLeft == matchRight)&&(j==0)))
      {
        e->match[i] = matchLeft;
        e->mismatchpoint[i] = mismatchpointLeft;
        e->isLeft[i] = true;
      }
      else
      {
        //fprintf(stderr, "----- Right match! \n");
        e->match[i] = matchRight;
        e->mismatchpoint[i] = mismatchpointRight;
        e->isLeft[i] = false;
      }

      /*LeftMatch4Nodes(e->onenode, nodes->at(id), match, mismatchpoint);
      e->match[i] = match;
      e->mismatchpoint[i] = mismatchpoint;
      */
    }
  return NULL;
}

int LeftMatch4Nodes(const COneNode *node1, const COneNode* node2, int &matchLeft, int &mismatchpointLeft)
{

  matchLeft = 0;
  mismatchpointLeft = -1;

  int l1, l2, r1, r2;

  node1->GetLR(l1, r1);
  node2->GetLR(l2, r2);

  int ed = MIN (r1, r2);

  int id1 = ed - l1;
  int id2 = ed - l2;

  const char* p1 = node1->GetData();
  const char* p2 = node2->GetData();

  while (id1 >= 0 && id2 >= 0)
    {
      if (NOTMISSING(p1[id1]) && NOTMISSING(p2[id2]))
	{
	  if (p1[id1] != p2[id2])
	    {
	      mismatchpointLeft = id1+l1;
	      return 0;
	    }
	  else matchLeft ++;
	}
      id1 --;
      id2 --;
    }

  return 1;
}

int RightMatch4Nodes(const COneNode *node1, const COneNode* node2, int &matchRight, int &mismatchpointRight)
{

  matchRight = 0;
  mismatchpointRight = -1;

  int l1, l2, r1, r2;

  node1->GetLR(l1, r1);
  node2->GetLR(l2, r2);

  int ed = MAX (l1, l2);

  int id1 = ed-l1;
  int id2 = ed-l2;

  int length1 = r1-l1;
  int lenght2 = r2-l2;

  const char* p1 = node1->GetData();
  const char* p2 = node2->GetData();

  while ((id1 <= length1) && (id2 <= lenght2))
    {
      if (NOTMISSING(p1[id1]) && NOTMISSING(p2[id2]))
	{
	  if (p1[id1] != p2[id2])
	    {
	      mismatchpointRight = id1+l1;
	      return 0;
	    }
	  else matchRight ++;
	}
      id1 ++;
      id2 ++;
    }
  return 1;
}




