
#ifndef __COneNode__
#define  __COneNode__
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <string.h>
#include "define.h"


#define SWAP01(x) ((x) == '0' ? '1' : ((x) == '1' ? '0':(x)))

class COneNode
{
 private:

  int id; //id in vector nodes


  char *data;
  int left, right; //the starting point in snp index
  //data[0] .. data[len] = original[left]..orignial[right] : len = right - left

  COneNode *child1,  *child2;
  COneNode *parent1, *parent2;

  bool mark, mut, leaf;

 public:

  ~COneNode();
  COneNode(COneNode * child1, COneNode * child2, int left, int right,
	   char *data, bool copy);


  COneNode *clone(void);

  void SetID(int iid);
  void Deactivate(void);

  bool Mutation(int site, int id, COneNode *&e);
  void Break(int siteid, COneNode* &e1, COneNode*& e2, bool matchLeft); //break and deactive

  int GetLeft(void) const;
  int GetRight(void) const;
  int GetLength(void) const;
  int GetID(void) const;

  const char* GetData(void) const;
  char * GetDataOut(void);  //get out data and set data = null

  bool GetMark(void) const;

  void View(void) const;
  int  GetDataLength(void) const;
  bool ExistSite(int k) const;
  char GetSite(int k) const;
  bool overlap(int il, int ir) const;
  void GetLR(int &l, int &r) const;

  void SetParents(COneNode *p1, COneNode *p2);

  void SetMut(void);

  void Replace(COneNode *c1, COneNode *c2);
  void RemoveNode(void);
  bool MutationNode(void) const;

  bool checkroot(void) const;
  bool checkleaf(void) const;

  void SetLeaf(void);

  void GetChild(const COneNode* &c1, const COneNode *&c2) const;
  void GetParents(const COneNode *&p1, const COneNode *&p2) const;
  bool Contain(int k) const;
  void ClearData(void);
  void ResetData(void);



};




typedef std::vector<COneNode*> vecpCOneNode;

#endif
