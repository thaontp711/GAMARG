#include "CARG.h"

CARG :: CARG(const CData *idata)
{
  numsequences = idata->GetnHap();
  nummarkers = idata->GetnSite();
}

void CARG :: getMarginalTrees(FILE *ftrees,int noArg)
{
       COneNode *root = NULL;
	for (int i = 0; i != (int) nodes.size(); ++i)
		if (nodes[i]->checkroot())
		{
			root = nodes[i];
			break;
		}
	for (int siteid = 0; siteid < nummarkers; siteid++)
	{
		getMarginalTreeForMarker(ftrees,root,siteid, noArg);
	}
      
}

void CARG :: getMarginalTreeForMarker(FILE *ftrees,COneNode *root, int siteid, int noArg)
{
	std::vector<const COneNode*> q;
	q.push_back(root);
	int marginaltrees[3][numsequences-1];
	fprintf(ftrees, "marker %d\n", siteid);
	int node = numsequences-2;
	int id = 0;
	int i = 0;
	while (i < (int) q.size())
	{
		const COneNode *e = q[i++], *p = NULL;
		const COneNode *c1 = NULL, *c2 = NULL;
		id = e->GetID();
		while (1==1)
		{
			e->GetChild(c1, c2);
			int nc = 0;
			if (c1 != NULL && c1->Contain(siteid))  { nc ++;  p = c1; }
			if (c2 != NULL && c2->Contain(siteid))  { nc ++;  p = c2; }
			if (nc == 1)
			{
				if(p->checkleaf()) {
					for(int x = numsequences-2; x >= node; x--)
					{
						if(marginaltrees[0][x] == id) {
							marginaltrees[0][x] = p->GetID();
							break; }
						else if(marginaltrees[1][x] == id) {
							marginaltrees[1][x] = p->GetID();
							break; }
					}
					
				}
				else {
					e = p;
					continue;
				}
		
			}

			if(nc ==2)
			{
				marginaltrees[0][node] = c1->GetID();
				marginaltrees[1][node] = c2->GetID();
				marginaltrees[2][node] = id;
				node--;
				
			}
			break;
		}
		if(e->checkleaf())
			continue;
		if (c1 != NULL && c1->Contain(siteid) && (!c1->checkleaf()))  q.push_back(c1);
		if (c2 != NULL && c2->Contain(siteid) && (!c2->checkleaf()))  q.push_back(c2);

	}

	vecint a,b;
	int n = numsequences-1;
	int array[2][numsequences-2];
	for(int i = 0; i<2;i++)
	      for(int j = 0; j < numsequences-1;j++)	
		array[i][j] = -1;
	
	for(int i = 0;i<numsequences-1;i++)
	{
		n++;
		if (i!=numsequences-2){
		    a.push_back(marginaltrees[2][i]); 
		}
		marginaltrees[2][i] = n;
		b.push_back(n);
	}
	
	for(int i = 0; i<(int)a.size();i++)
	{
		int x = a[i];

		for(int j = numsequences-2;j>=0;j--)
		{
			if((marginaltrees[0][j] == x) && (array[0][j] == -1)){
				marginaltrees[0][j] = b[i];
				array[0][j] = 1;
				break; }
			else if((marginaltrees[1][j] == x) && (array[1][j] = -1)) {
				marginaltrees[1][j] = b[i];
				array[1][j] = 1;
				break; }
		}
	}
	
	for(int i = 0; i< numsequences-1;i++)
		fprintf(ftrees,"%d %d %d\n",marginaltrees[0][i],marginaltrees[1][i],marginaltrees[2][i]);
}


void CARG :: Load(vecpCOneNode *nodes, const vecint &poss)
{

  




}

