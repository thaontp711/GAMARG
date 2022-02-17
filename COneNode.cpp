#include "COneNode.h"



COneNode ::  COneNode(COneNode * ichild1, COneNode * ichild2, int ileft, int iright,
		      char* idata, bool copy)
{
  mut = false;
  mark = false;
  leaf = false;
  child1 = ichild1;
  child2 = ichild2;

  left   = ileft;
  right = iright;

  parent1 = NULL;
  parent2 = NULL;


  if (copy)
    {
      data = (char*)calloc(right - left + 3, sizeof(char));
      strncpy(data, idata, (right-left) + 1);
      data[ right - left +1] = '\0';
    }
  else  {
    data = idata;
    if (data != NULL)	data = (char*)realloc(data, sizeof(char) * (right - left + 3));
  }

  id = -1;
}


COneNode :: ~COneNode()
{
  if (data != NULL) { free(data); data = NULL;}
}


int COneNode :: GetID(void) const
{
  return id;
}


void COneNode ::  SetLeaf(void)
{
  leaf = true;
}

void COneNode :: SetID(int iid)
{
  id = iid;
}


int COneNode :: GetLeft(void) const
{
  return left;
}



int COneNode :: GetRight(void) const
{
  return right;
}



int COneNode :: GetLength(void) const
{
  return (right - left + 1);
}

void COneNode :: GetChild(const COneNode* &c1, const COneNode *&c2) const
{
  c1 = child1;
  c2 = child2;
}

const char* COneNode :: GetData(void) const
{
  return data;
}


char* COneNode :: GetDataOut(void)
{
  char *p = data;
  data = NULL;
  return p;
}


bool COneNode :: GetMark(void) const
{
  return mark;
}

bool COneNode :: Mutation(int site, int id, COneNode *&e)
{

  if (mark) return false;
  if (site < left || site > right) return false;

  int siteid = site - left;

  if (data[ siteid] != id) return false;
  if (data[ siteid] == '0') data[siteid] = '1';
  else
    if (data[ siteid] == '1') data[siteid] = '0';

  e = new COneNode (this, NULL, left, right, data, NOCOPY);
  parent1 = e;
  data = NULL;
  mark = 1;
  e->SetMut();
  return true;
}


void COneNode :: ResetData(void)
{
  data = NULL;
}

/*
void COneNode :: ClearData(void)
{
  if (data != NULL) { free(data); data = NULL;}
}
*/

void COneNode :: Deactivate(void)
{
  mark = true;
}

void COneNode :: SetParents(COneNode *p1, COneNode *p2){
  parent1 = p1;
  parent2 = p2;
}


bool COneNode :: Contain(int k) const
{
  return (k >= left && k <= right);
}

bool COneNode :: ExistSite(int k) const
{
  return (!mark && k <= right && k >= left);
}

char COneNode :: GetSite (int k) const
{
  return (data[k - left]);
}

void COneNode :: View(void) const
{
  fprintf(stderr, "%d\t%d\t%3d\t", left, right, id);
  if (data != NULL)
    {
      for (int i = 0; i <= right; ++i) 
	if (i < left) fprintf(stderr, ".");
	else fprintf(stderr,"%c", data[i-left]);

    }
      fprintf(stderr,"\n");
}


void COneNode :: GetLR(int &l, int &r) const
{
  l = left;
  r = right;
}

void COneNode :: Replace(COneNode *c1, COneNode *c2)
{
  if (c1 == child1) child1 = c2;
  else
    if (c1 == child2) child2 = c2;
    else
      {
	printf("error in replace child\n");
	exit(1);
      }

}


bool COneNode :: MutationNode(void) const
{
  return mut;
}

bool COneNode :: checkroot (void) const
{
  return (parent1 == NULL && parent2 == NULL);
}


bool COneNode :: checkleaf (void) const
{
  return leaf;
}



void COneNode :: RemoveNode(void)
{
  COneNode * child = child1;
  if (child == NULL) child = child2;
  if ((child1 == NULL && child2 == NULL) || (child1 != NULL && child2 != NULL))
    {
      printf("error in Removeitselft\n");
      exit(1);
    }

  child->SetParents(parent1, parent2);

  if (parent1 != NULL) parent1->Replace(this, child);
  if (parent2 != NULL) parent2->Replace(this, child);

  child1 = NULL;    child2 = NULL;
  parent1 = NULL;   parent2 = NULL;
}

void COneNode :: Break(int siteid, COneNode* &e1, COneNode* &e2, bool matchLeft)
{
  if (siteid > right || siteid < left)
    {
      fprintf(stderr, "error in break: %d\n", siteid);
      View();
      exit(1);
    }

    if(matchLeft == true)
    {
        e2 = new COneNode (this, NULL, siteid+1, right,  &(data[siteid+1-left]), COPY);
        e1 = new COneNode (this, NULL, left,     siteid, data, NOCOPY); //data of node1
    }
    else
    {
        e1 = new COneNode (this, NULL, siteid, right,  &(data[siteid-left]), COPY);
        e2 = new COneNode (this, NULL, left,     siteid-1, data, NOCOPY); //data of node1
    }


  parent1 = e1;
  parent2 = e2;
  data  = NULL;
  mark = 1;
}

int COneNode :: GetDataLength(void) const
{
  return (int) strlen(data);
}

bool COneNode :: overlap(int il, int ir) const
{
  if (right < il) return false;
  if (ir < left) return false;
  return true;
}

void COneNode :: SetMut(void)
{
  mut = true;
}


void COneNode :: GetParents(const COneNode *&p1, const COneNode *&p2) const
{
  p1 = parent1;
  p2 = parent2;
}
