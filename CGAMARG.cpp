#include "CGAMARG.h"





CGAMARG::    CGAMARG(int nth, int ndist)
{
  
  nthread = nth;
  distgam = ndist;
  fprintf(stderr, "fm->nthread: %d; distgam = %d\n", nthread, distgam);
  nRec = 0;
  priqueue = (prioqueue_onepair*) new prioqueue_onepair;
  mu = true;
}


CGAMARG::    ~CGAMARG()
{
  for (int i = 0; i != (int) nodes.size(); ++i) delete nodes[i];
  
  if (priqueue != NULL) free(priqueue);
}


bool operator < (const onepair &x, const onepair &y)
{
  return x.match < y.match;
}

bool operator < (const gemsite &x, const gemsite &y)
{
	return x.sright - x.sleft > y.sright - y.sleft;
		//return x.matchlength < y.matchlength;
	//else
	//	return x.gametecount < y.gametecount;
}

void CGAMARG :: FindPath(int id) const
{
  /*
  std::vector<const char*> haps;
  const vecstring *hs = data->GetHaps();

  const char* comhap = NULL;

  for (int i = 0; i != hs->size(); ++i)
    if (i != id) haps.push_back(hs->at(i).c_str());
    else comhap = hs->at(i).c_str();

  vecint list1, list2;
  vecint *listid, *listidn;
  bool one = true;
  int nsite = strlen(comhap), nrc = 0, nhap = haps.size(), siteid = 0;

  while (siteid < nsite)
    {

      //      printf("siteid=%d; listid1 = %d; list2 = %d one = %d\n", siteid, list1.size(), list2.size(), one?1:0);
      if (one) { listid = &list1; listidn = &list2;}
      else {listid = &list2; listidn = &list1;}
      one = !one;



      if (listid->size() == 0)
	{
	  list1.clear();
	  list2.clear();
	  for (int i = 0; i != haps.size(); ++i)    list1.push_back(i);
	  one = true;
	  nrc ++;
	  continue;
	}


      listidn->clear();
      for (int i = 0; i != listid->size(); ++i)
	{
	  int id = listid->at(i);
	  if (haps[id][siteid] == comhap[siteid]) listidn->push_back(id);
	}

      if (listidn->size() > 0) siteid ++;
      else
	{
	  bool mut = true;
	  char chr = '.';
	  for (int i = 0; i != nhap; ++i)
	    if (haps[i][siteid] == '0' || haps[i][siteid] == '1')
	      {
		if (chr == '.') chr = haps[i][siteid];
		if (chr != haps[i][siteid]) { mut = false; break;}
	      }
	  if (mut) { one = !one; siteid++;}

	}
    }
  printf("id = %d, nrc = %d\n", id, nrc);
  */
}


void CGAMARG :: CheckCount(void) const
{
  int nsite = data->GetnSite();
  for (int i = 0; i != nsite; ++i)
    {
      std::string s = "";
      for (int j = 0; j != (int) activenodes.size(); ++j) {
	int id  = activenodes[j];
	if (nodes[id]->ExistSite(i)) s = s + nodes[id]->GetSite(i);
      }

      int c0 =0, c1 = 0, c2 = 0;
      for (int k = 0; k != (int) s.size(); ++k)
	{
	  if (s[k] == '0') c0 ++;
	  if (s[k] == '1') c1 ++;
	  if (!NOTMISSING(s[k])) { c2 ++; fprintf(stderr, "%d: %d\n", k, s[k]);}
	}
      if (c0 != count0[i] || c1 != count1[i])
	{
	  fprintf(stderr, "error: %d: %d\t%d\t%s [%d,%d, %d]\n", i, count0[i], count1[i], s.c_str(), c0,c1, c2);
	  exit(1);
	}
    }
}


void CGAMARG :: ReadData(const CData *idata)
{

  data = idata;

  nodes.clear();
  activenodes.clear();
  coallist.clear();
  prigemlist = std::priority_queue <gemsite>();
    vecop.clear();
   
    fprintf(stderr, "prigemlist size: %d\n",(int)prigemlist.size());

  count0.clear();
  count1.clear();
    count00.clear();
    count01.clear();
    count10.clear();
    count11.clear();
    sitemarked.clear();

  int nsite = data->GetnSite();

  count0.resize(nsite, 0);
  count1.resize(nsite, 0);
    count00.resize((nsite-1)*distgam, 0);
    count01.resize(distgam*(nsite-1), 0);
    count10.resize(distgam*(nsite-1), 0);
    count11.resize(distgam*(nsite-1), 0);
    sitemarked.resize(distgam*(nsite-1),false);
    fprintf(stderr, "resize: %d\n",(int)prigemlist.size());

  const vecstring *haps = data->GetHaps();

  for (int i = 0; i != (int) haps->size(); ++i)
    {
      const char *p = haps->at(i).c_str();
      for (int k = 0; k != nsite; ++k)
	{
	  if (p[k] == '0') count0[k] ++;
	  if (p[k] == '1') count1[k] ++;
    }
        for (int k = 0; k != nsite-1; ++k)
        {
            
                for (int j = 1; j != distgam+1; j++)
                {
                    if (j+k > nsite-1)
                        break;
                    //fprintf(stderr,"k = %d, k+j = %d, p[k] = %c, p[k+j] = %c\n", k, k+j, p[k], p[k+j]);
                    if ((p[k] == '0')&&(p[k] == p[k+j])) count00[k*distgam+j]++;
                    if ((p[k] == '0')&&(p[k] != p[k+j])) count01[k*distgam+j]++;
                    if ((p[k] == '1')&&(p[k] == p[k+j])) count11[k*distgam+j]++;
                    if ((p[k] == '1')&&(p[k] != p[k+j])) count10[k*distgam+j]++;
                }
                
            }
            
        }
        
    

  char* tempc = (char*)calloc(M20, sizeof(char));

  for (int i = 0; i != (int) haps->size(); ++i)
    {
      strcpy(tempc, haps->at(i).c_str());
      COneNode *e = new COneNode(NULL, NULL, 0, nsite-1, tempc, COPY);
      e->SetLeaf();
      Add(e);
    }
  free(tempc);
    for (int k = 0; k != nsite; ++k)
        if (POSSIBLEMUT(count0[k], count1[k])) mutlist.push_back(k); // mutlist

    for (int k = 0; k != nsite-1; ++k)
        for (int j = 1; j != distgam+1; j++)
        {
            if (j+k > nsite-1)
                break;

            if((count11[k*distgam+j] != 0)&&(count00[k*distgam+j] != 0)&&(count10[k*distgam+j] != 0)&&(count01[k*distgam+j] != 0)) {
              if ((count11[k*distgam+j] == 1)||(count10[k*distgam+j] == 1)||(count00[k*distgam+j] ==1)||(count01[k*distgam+j] == 1)) {
                gemsite gs;
                if (count11[k*distgam+j] == 1) {
                    gs.gamete = 3;
                    gs.sleft = k;
                    gs.sright = k+j;
                    prigemlist.push(gs);
                    
                }
                if (count10[k*distgam+j] == 1) {
                    gs.gamete = 2;
                    gs.sleft = k;
                    gs.sright = k+j;
                    prigemlist.push(gs);
                }
                if (count01[k*distgam+j] == 1) {
                    gs.gamete = 1;
                    gs.sleft = k;
                    gs.sright = k+j;
                    prigemlist.push(gs);
                }
                if (count00[k*distgam+j] == 1) {
                    gs.gamete = 0;
                    gs.sleft = k;
                    gs.sright = k+j;
                    prigemlist.push(gs);
                }
                //gs.site = k;
                //gemlist.push_back(gs);
                sitemarked[k*distgam+j] = true;
                  //fprintf(stderr, "gs = %d, sitemarked = %d\n", gs.gamete, (int)sitemarked[k*distgam+j]);
            }
                
        }
    }
	fprintf(stderr, "prigemlist size = %d, sitemarked = %d\n", (int)prigemlist.size(), (int)sitemarked.size());
}

void CGAMARG :: AddNodetoThread()
{
  std::vector<int> iret; //thread id
  iret.resize(nthread, 0);

  std::vector<pthread_t> threads; //pthreads
  threads.resize(nthread, 0);

  std::vector<onethread> ths; //data for one thread
  
  onethread e;
  e.nodes = &nodes;
  //e.coallist = &coallist;
  //e.priqueue = &priqueue;
  //e.onenode = onenode;
  e.activenodes.clear();
  ths.resize(nthread, e);

  //fprintf(stderr, "nodetothread size: %d\n",(int)nodetothread.size());
  for (int i = 0; i != (int) ths.size(); i++)
    {
	if((int)nodetothread.size() > i)
	{	
		COneNode *onenode = nodetothread[i];
		//fprintf(stderr, "-------- Add node: %d\n",i);
		//onenode->View();
		ths[i].onenode = onenode;
		//ths[i].onenode->View();
		for (int j = 0; j != (int) activenodes.size(); ++j)
		    {
			int id = activenodes[j];
			if (onenode->overlap( nodes[id]->GetLeft(), nodes[id]->GetRight()))
    			  ths[i].activenodes.push_back(id);
			//fprintf(stderr, "AAAA %d\n",(int)activenodes.size());
		    }
	}
	else
	    break;
	//fprintf(stderr, "%d------nodetothread size: %d\n",i,(int)nodetothread.size());
    }
  for (int i = 0; i != (int) ths.size(); ++i)
    {
      int nr = (int) ths[i].activenodes.size();
      ths[i].match.resize(nr, 0);
      ths[i].mismatchpoint.resize(nr, 0);
      ths[i].isLeft.resize(nr, 0);
      //ths[i].difcharn1.resize(nr, 0);
    }


  for (int i = 0; i != nthread; ++i)
    iret[i] = pthread_create(&threads[i], NULL, LeftMatch, (void*) &ths[i]);

  for (int i = 0; i != nthread; ++i) pthread_join(threads[i], NULL);


  for (int i = 0; i != (int) ths.size(); ++i)
    for (int j = 0; j != (int) ths[i].activenodes.size(); ++j)
      {
      	 if((int)nodetothread.size() > i)
		 {
			int id = (int) ths[i].activenodes[j];
			onepair e;
			//e.node1 = onenode;
			e.node1 = nodetothread[i];
			e.node2 = nodes[id];
			e.match = ths[i].match[j];
			e.mismatchpoint = ths[i].mismatchpoint[j];
			e.matchLeft = ths[i].isLeft[j];
			//e.node1chardif = ths[i].difcharn1[j];
			
			if (e.mismatchpoint == -1 && e.match > 0) 
			{ 
				coallist.push_back(e);}
			if (e.match > 0) priqueue->push(e);
		 }
      }
  
  for (int i = 0; i != (int) nodetothread.size(); ++i)
  {
	  COneNode *onenode = nodetothread[i];
	  onenode->SetID((int) nodes.size());
	  nodes.push_back(onenode);
	  activenodes.push_back(onenode->GetID());
  }
 
  nodetothread.clear();
 // CheckCount();
}

void CGAMARG :: Add(COneNode* &onenode)
{

  std::vector<int> iret; //thread id
  iret.resize(nthread, 0);

  std::vector<pthread_t> threads; //pthreads
  threads.resize(nthread, 0);

  std::vector<onethread> ths; //data for one thread

  onethread e;
  e.nodes = &nodes;
  e.onenode = onenode;
  e.activenodes.clear();
  ths.resize(nthread, e);


  //add nodes into
  int nc = 0;
  for (int i = 0; i != (int) activenodes.size(); ++i)
    {
      int id = activenodes[i];
      if (onenode->overlap( nodes[id]->GetLeft(), nodes[id]->GetRight()))
	ths[nc ++ % nthread].activenodes.push_back(id);
    }


  for (int i = 0; i != (int) ths.size(); ++i)
    {
      int nr = (int) ths[i].activenodes.size();
      ths[i].match.resize(nr, 0);
      ths[i].mismatchpoint.resize(nr, 0);
      ths[i].isLeft.resize(nr, 0);
      //ths[i].difcharn1.resize(nr, 0);
    }


  for (int i = 0; i != nthread; ++i)
    iret[i] = pthread_create(&threads[i], NULL, LeftMatch, (void*) &ths[i]);

  for (int i = 0; i != nthread; ++i) pthread_join(threads[i], NULL);


  for (int i = 0; i != (int) ths.size(); ++i)
    for (int j = 0; j != (int) ths[i].activenodes.size(); ++j)
      {
	int id = (int) ths[i].activenodes[j];
	onepair e;
	e.node1 = onenode;
	e.node2 = nodes[id];
	e.match = ths[i].match[j];
	e.mismatchpoint = ths[i].mismatchpoint[j];
	e.matchLeft = ths[i].isLeft[j];
	//e.node1chardif = ths[i].difcharn1[j];
    /*if(e.mismatchpoint != -1 && e.match > 0) {
        bool gamete = CheckGamete(e.mismatchpoint, e.matchLeft, id, e.gametetypenode1, e.gametetypenode2, e.countgametenode1, e.countgametenode2);
        e.gamete = gamete;
    }
    else
        e.gamete = false;
     */
	if (e.mismatchpoint == -1 && e.match > 0)
	{ 
		/*fprintf(stderr, "====== Missmatchpoint: %d ---- %d ========\n", e.mismatchpoint, e.matchLeft);
		e.node1->View();
		e.node2->View();	*/	
		coallist.push_back(e);}
	if (e.match > 0) priqueue->push(e);
      }


  onenode->SetID((int) nodes.size());
  /*fprintf(stderr, "\n-------- View all active nodes --------\n");
  for (int i = 0; i != (int) activenodes.size(); ++i)
  {
      int id = activenodes[i];
	nodes[id]->View();
  }*/

//fprintf(stderr, "-------- Add node: --------\n");
//onenode->View();
  nodes.push_back(onenode);
  activenodes.push_back(onenode->GetID());

  /*fprintf(stderr, "\n-------- View all active nodes --------\n");
  for (int i = 0; i != (int) activenodes.size(); ++i)
  {
    int id = activenodes[i];
    nodes[id]->View();
  }*/

}



void CGAMARG :: Mutation(int siteid, int id)
{
  bool mut = false;
  for (int i = 0; i != (int) activenodes.size(); ++i)
    {
      int nid = activenodes[i];
//fprintf(stderr, "+++ node %d, siteid: %d, id: %d \n",nodes[nid]->GetID(), siteid, id);

      COneNode *e = NULL;
      if (nodes[nid]->Mutation(siteid, id, e))
	{
	  //fprintf(stderr, "Mutation %d ---- %d !!!!!!!!\n", nodes[nid]->GetID(), siteid);

        const char *p = e->GetData();
        //e->View();
        
        int l, r;
        e->GetLR(l,r);
        //const char *p = nodes[nid]->GetData();
        //nodes[nid]->View();
        
        //fprintf(stderr, "\nleft: %d, right: %d, p[siteid-l] = %c \n", l, r,p[siteid-l]);
        //fprintf(stderr, "00 = %d, 01 = %d, 10 = %d, 11 = %d\n", count00[siteid], count01[siteid], count10[siteid], count11[siteid]);
        if(nodes[nid]->GetLength() > 1) {
            if (siteid == l) { /// compare siteid to siteid+1, siteid+2, ....siteid+distgam
                int end = 0;
                if (r-siteid > distgam)
                    end = distgam;
                else
                    end = r-siteid;
                
                //int siteidj = siteid+1;
                for (int j = 1; j <= end; j++)
                    CheckGameteMut1st(p[siteid-l], p[siteid+j-l],siteid,siteid+j);
                    //CheckGameteMut1st(p[siteid-l], p[siteid+j-l],siteid*distgam+j);
                //CheckGameteMut1st(p[siteid-l],p[siteid+1-l],siteid);
            }
            
            else if (siteid == r) { /// compare siteid to siteid-1, siteid-2, ....siteid-distgam
                int start = 0;
                if (siteid -l > distgam)
                    start = siteid-distgam;
                else
                    start = l;
                for (int j = start; j < siteid; j++)
                    CheckGameteMut2nd(p[j-l],p[siteid-l],j, siteid);
                    //CheckGameteMut2nd(p[j-l],p[siteid-l],j*distgam+siteid-j);
            }
            
            else {  // compare to both left and right site of siteid
                // from siteid to the right
                int end = 0;
                if (r-siteid > distgam)
                    end = distgam;
                else
                    end = r-siteid;
                
                for (int j = 1; j <= end; j++)
                    CheckGameteMut1st(p[siteid-l], p[siteid+j-l],siteid,siteid+j);
                    //CheckGameteMut1st(p[siteid-l], p[siteid+j-l],siteid*distgam+j);
                
                //from siteid back to the left
                int start = 0;
                if (siteid -l > distgam)
                    start = siteid-distgam;
                else
                    start = l;
                for (int j = start; j < siteid; j++)
                    CheckGameteMut2nd(p[j-l],p[siteid-l],j,siteid);
                    //CheckGameteMut2nd(p[j-l],p[siteid-l],j*distgam+siteid-j);
                
                //CheckGameteMut1st(p[siteid-l],p[siteid+1-l],siteid);
                //CheckGameteMut2nd(p[siteid-1-l],p[siteid-l],siteid-1);
            }
        }
//nodes[nid]->View();
	  //fprintf(stderr, "Mutation %d ---- %d !!!!!!!!\n", nodes[nid]->GetID(), siteid);
	  RemoveActiveID( nodes[nid]->GetID()); //must be done before add e as nodes[nid] has no data
        nodes[nid]->Deactivate();
	  Add(e);
	  mut = true;
	  break;
	}
    }

  if (!mut && 1==1)
    {
      fprintf(stderr, "cannot find mutation\n");
      exit(1);
    }
  if (id == (int)'0') { count0[siteid]  --; count1[siteid] ++;}
  else {count1[siteid]  --; count0[siteid] ++;}
  //for (int i = 0; i < (int)gemlist.size(); i++)
	//fprintf(stderr, "gem: type: %d, site: %d, mark = %d\n", gemlist[i].gamete, gemlist[i].site, (int)sitemarked[gemlist[i].site]);
  //fprintf(stderr, "id = %d, count00 = %d, count01 = %d, count10 = %d, count11 = %d, mark = %d\n",siteid, count00[siteid],count01[siteid],count10[siteid],count11[siteid], (int)sitemarked[siteid]);
}


bool CGAMARG :: PossibleMutation(void)
{

  bool mut = false;
  std::random_shuffle(mutlist.begin(), mutlist.end());
  for (int ii = 0; ii != (int) mutlist.size(); ++ii)
    {
      int i = mutlist[ii];
      if (POSSIBLEMUT( count0[i], count1[i]))
	{
	  if (count1[i] == 1)	    Mutation(i, '1');
	  else Mutation(i, '0');
	  mut = true;
	}
    }
  mutlist.clear();
  return mut;
}




void CGAMARG :: RemoveActiveID(int id)
{
  for (int i = 0; i != (int) activenodes.size(); ++i)
    if (id == activenodes[i])
      {
	activenodes[i] = activenodes[ activenodes.size()-1 ];
	activenodes.resize(activenodes.size()-1);
	return  ;
      }
}


void CGAMARG :: Coal(COneNode *node1, COneNode* node2)
{

  int l1, r1, l2, r2;
  node1->GetLR(l1, r1);
  node2->GetLR(l2, r2);
  if (l1 < l2 || (l1 == l2 && r1 > r2))
    {
      COneNode *t = node1;
      node1 = node2;
      node2 = t;
    }

  //l1 > l2 || (l1 == l2 && r1 <= r2)

  node1->GetLR(l1, r1);
  node2->GetLR(l2, r2);
//fprintf(stderr, "Coal!!!!!! \n");
//View();
//fprintf(stderr, "+++++++++++++++++++++++++++++++++++++++++++ \n");
/*
fprintf(stderr, "node1: ");
fprintf(stderr, "l1: %d, r1: %d, data: \n", l1, r1);
node1->View();
fprintf(stderr, "node2: ");
fprintf(stderr, "l2: %d, r2: %d, data: \n", l2, r2);
node2->View();
*/
  char *p1 = node1->GetDataOut(); //get data out of node 1
  char *p  = node2->GetDataOut(); //get data out of node 2

  if (r2 < r1)
    p = (char*) realloc (p, (r1 - l2 + 3) * sizeof(char)); //extend p1

  p[ MAX(r1, r2) - l2 + 1] = '\0'; //

  int id1 = 0, id = l1 - l2;
  for (int i = l1; i <= r1; ++i)
    {
      if ((id <= (r2-l2)) && NOTMISSING(p[id]) && NOTMISSING(p1[id1]))
	{
	  if (p[id] == '0') count0[i] --;
	  else count1[i] --;

	  if (POSSIBLEMUT( count0[i], count1[i])) 
	  { 	//fprintf(stderr, "+++++++++++++++++++++ Add mutation: %d, 1: %d, 0: %d \n",i, count1[i], count0[i]);
		mutlist.push_back(i);
	  }
        if ((i < r1) &&(id < (r2-l2))&&(node1->GetLength() > 1)&&(node2->GetLength() > 1)) {
            int end = 0;
            if (r2 < r1)
                end = r2;
            else
                end = r1;
            if (end-i > distgam)
                end = distgam;
            else
                end = end -i;
            
            //int idj = id +1;
            
            for (int j = 1; j <= end; j++) {
		//fprintf(stderr, "i = %d, id = %d, end = %d, p[id] = %c, p[id+j] = %c\n",i, id, end, p[id], p[id+j]);
                CheckGamete(p[id], p[id+j],i,i+j);
                //CheckGamete(p[id], p[id+j],i*distgam+j);
                //idj++;
            }
		//fprintf(stderr, "id = %d, count00 = %d, count01 = %d, count10 = %d, count11 = %d, mark = %d\n",id, count00[i],count01[i],count10[i],count11[i], (int)sitemarked[i]);
	 }
	
	}
        else if ((l1 > l2) && (id > (r2-l2)) && NOTMISSING(p1[id1]))
        {
            
            if ((i-l1) < distgam) {
                int start = 0;
                if ((i-l2) > distgam)
                    start = i-distgam;
                else
                    start = l2;
                //fprintf(stderr, "i = %d, start = %d, id = %d, id1 = %d\n", i, start, id, id1);
                for (int x = start; x < l1; x++)
                    CheckGameteCoal(p[x-l2],p1[id1],x,i);
                
            }
        }

      //if (!NOTMISSING(p[id]))
      p[id] = p1[id1];
      id1 ++;
      id  ++;
    }

  free(p1);

  /*for (int i = 0; i < (int)gemlist.size(); i++)
	fprintf(stderr, "gem: type: %d, site: %d, mark = %d\n", gemlist[i].gamete, gemlist[i].site, (int)sitemarked[gemlist[i].site]);
  */
    //fprintf(stderr,"done Coal +++\n");

  COneNode *e = new COneNode(node1, node2, l2, MAX(r1, r2), p, NOCOPY);
//int r;
//r = e->GetRight();
  //fprintf(stderr,"Coalescent node: \n");
//fprintf(stderr, "l: %d, r: %d, data: \n", l2, r);
  //e->View();
  RemoveActiveID( node1->GetID() );
  RemoveActiveID( node2->GetID() );


  node1->SetParents(e, NULL);
  node2->SetParents(e, NULL);

  node1->Deactivate();
  node1->ResetData();

  node2->Deactivate();
  node2->ResetData();

   Add(e);
}

void CGAMARG :: Recom(COneNode *node1,  COneNode* node2, int matchpoint, bool matchLeft)
{
    //fprintf(stderr, "Recom!!!!!! \n");

    /*if(matchLeft == true)
    {
        printf("Left %d\n", matchpoint);
        if (node1->GetRight() > node2->GetRight() || (node1->GetRight() == node2->GetRight() && node1->GetLeft() < node2->GetLeft())) //break the smaller node
        {
            COneNode *t = node1; node1 = node2; node2 = t;
        }
    }
    else
    {
        printf("Right %d\n", matchpoint);
        if (node1->GetLeft() < node2->GetLeft()||(node1->GetLeft() == node2->GetLeft() && node1->GetRight() > node2->GetRight())) //break the smaller node
        {
          COneNode *t = node1; node1 = node2; node2 = t;
        }
    }
*/
  int i = rand()%2;
  if ((node1->GetLength() > node2->GetLength()) || ((node1->GetLength() == node2->GetLength())&&(i==0))) //break the smaller node
    {
      COneNode *t = node1; node1 = node2; node2 = t;
    }
    
    const char* p = node1->GetData();
    int l,r;
	node1->GetLR(l,r);
    //node1->View();
    //node2->View();
    int point = 0;
    int countpoint = 0;
    int endpoint = 0;
    
    if (matchLeft == true)
        point = matchpoint;
    else
        point = matchpoint - 1;
    //fprintf(stderr, "matchpoint = %d\n",point);
    
    if (point-l >= distgam)
        countpoint = point - distgam +1;
    else
        countpoint = l;
    
    if (r-point > distgam)
        endpoint = point + distgam;
    else
        endpoint = r;
    
    for (int i = countpoint; i <= point; i++)
        for (int j = point+1; j <= endpoint; j++)
        {
            if (j-i > distgam)
                break;
            CheckGamete(p[i-l],p[j-l],i,j);
        }
    
   
    
  COneNode *e1 = NULL, *e2 = NULL;
/*fprintf(stderr,"---------------------RECOM ------------------ \n");

  fprintf(stderr,"Node1: \n");
  node1->View();
  fprintf(stderr,"Node2: \n");
  node2->View();
*/
  node1->Break(matchpoint, e1, e2, matchLeft);
  RemoveActiveID(node1->GetID());
    node1->Deactivate();

  //fprintf(stderr,"Node e1: \n");
  //e1->View();

  Add(e1);
  //fprintf(stderr,"Node e2 to coalesent: \n");
  //e2->View();
  e2->SetID(nodes.size()); nodes.push_back(e2);
  Coal(e2, node2);
}

void CGAMARG :: RecomGemnode(COneNode *node,  int matchpoint, bool matchLeft)
{
    //fprintf(stderr, "Recom Gem!!!!!! \n");
	COneNode *e1 = NULL, *e2 = NULL;
  /*fprintf(stderr,"---------------------RECOM ------------------ \n");

  fprintf(stderr,"Node1: \n");
  node->View();
  //fprintf(stderr,"Node2: \n");
  //node2->View();
*/
    //node->View();
  node->Break(matchpoint, e1, e2, matchLeft);
  RemoveActiveID(node->GetID());
  node->Deactivate();

  Add(e1);
  Add(e2);

  //fprintf(stderr,"Node e1: \n");
  //e1->View();
  //fprintf(stderr,"Node e2: \n");
  //e2->View();
 
}

bool CGAMARG :: PossibleCoal(void)
{
  //  printf("\n\n--- check coal : %d --\n", (int) coallist.size());
    while (!coallist.empty())
    {
      std::random_shuffle(coallist.begin(), coallist.end());
      onepair e = coallist.back();
      coallist.pop_back();
      COneNode *node1 = e.node1, *node2 = e.node2;
      if (node1->GetMark() || node2->GetMark()) continue;
      Coal(node1, node2);
      return true;
    }
  return false;

}

bool CGAMARG :: PossibleGem(void)
{
    int lr = 0;
	//int countvecempty = 0;
    std::vector<gemsite> vecgs;
    gemsite temp;
    temp.gamete = -1;
    /*fprintf(stderr, "\n-------- View all active nodes --------\n");
    for (int i = 0; i != (int) activenodes.size(); ++i)
    {
        int id = activenodes[i];
        nodes[id]->View();
    }*/
    
    while (!prigemlist.empty())
    {
        gemsite gs = prigemlist.top();
        prigemlist.pop();
        
	//for (int i = 0; i < (int)gemlist.size(); i++)
	   //fprintf(stderr, "gem: type: %d, left: %d, right = %d, mark = %d\n", gs.gamete, gs.sleft, gs.sright, (int)sitemarked[gs.sleft*distgam+gs.sright-gs.sleft]);

        //std::random_shuffle(gemlist.begin(), gemlist.end());
	//for (int i = 0; i < (int)gemlist.size(); i++)
	    //fprintf(stderr, "After shuffle gem: type: %d, site: %d, mark = %d\n", gemlist[i].gamete, gemlist[i].site, (int)sitemarked[gemlist[i].site]);

        //gemsite gs = gemlist.back();
        //gemlist.pop_back();
	 
        int dist = gs.sright - gs.sleft;
	
        if (sitemarked[gs.sleft*distgam+dist] == false)
            continue;

        //fprintf(stderr, "sleft: %d, sright: %d, count00 = %d, count01 = %d, count10 = %d, count11 = %d, mark = %d\n",gs.sleft, gs.sright, count00[gs.sleft*distgam+dist],count01[gs.sleft*distgam+dist],count10[gs.sleft*distgam+dist],count11[gs.sleft*distgam+dist], (int)sitemarked[gs.sleft*distgam+dist]);

        if (lr == 0) lr = dist;
        if (dist > lr) { temp = gs; break; }
        if (dist == lr) vecgs.push_back(gs);
    }
    
    if (vecgs.size() == 0) return false;
      
      if ((temp.gamete != -1) && (sitemarked[temp.sleft*distgam+temp.sright-temp.sleft] == true))  prigemlist.push(temp);
      
      
      vecint vecchecked; // 0: not checked yet, 1: checked and added to vecop; 2: checked but marked
      vecop.clear();
	//vec_option vecoplr; // for nodes that mismatchpoint is not between sleft and sright
	//vecoplr.clear();
      vec_onepair vecapair;
      for (int i = 0; i < (int)vecgs.size(); i++)
          vecchecked.push_back(0);
      
      //int len = 0;
	int countcheck = 0;
	
      while(!priqueue->empty())
      {
        onepair e = priqueue->top();
        priqueue->pop();
	 if ((e.match == 0) || (e.mismatchpoint == -1)) continue;

          COneNode *node1 = e.node1, *node2 = e.node2;
	  if (node1->GetMark() || node2->GetMark()) continue;

          int l1,r1,l2,r2;
          node1->GetLR(l1,r1);
          node2->GetLR(l2,r2);
		
          for(int i = 0; i < (int)vecgs.size(); i++)
          {
              if(vecchecked[i] != 0)
                  continue;

              COneNode *nodegem;
              if(IsNodegamete(node1,vecgs[i]))
                  nodegem = node1;
              else if(IsNodegamete(node2,vecgs[i]))
                  nodegem = node2;
              else
                  continue;
                  		
              oneoption op;
              op.id = nodegem->GetID();
              op.idvecgs = i;
              if (e.mismatchpoint > vecgs[i].sright || e.mismatchpoint < vecgs[i].sleft) {
                  //fprintf(stderr,"Not Match mismatchpoint!!!\n");
                  int k = rand() % 2;
                  if (k == 0)
        			op.brpointgs = vecgs[i].sleft;
                  else
        			op.brpointgs = vecgs[i].sright-1;
                  //vecop.push_back(op);
              }
              else {
                  //fprintf(stderr,"Match mismatchpoint!!!\n");
                  if (e.mismatchpoint == vecgs[i].sleft)
                    op.brpointgs = vecgs[i].sleft;
                  else if (e.mismatchpoint == vecgs[i].sright)
                    op.brpointgs = vecgs[i].sright-1;
                  else {
                    
                    if (e.matchLeft == true)
                        op.brpointgs = e.mismatchpoint;
                    else
                        op.brpointgs = e.mismatchpoint-1;
                  }
              //vecop.push_back(op);
	
            }
              vecop.push_back(op);
            //fprintf(stderr,"op.brpointgs = %d\n", op.brpointgs);
              
            vecchecked[i] = 2;
            countcheck++;
		
              //nodegem->View();
              
          }
        //}
        if (!e.node1->GetMark() && !e.node2->GetMark())  vecapair.push_back(e);
        if (countcheck == (int)vecgs.size())
          break;
	
      }
      
      if (vecop.empty()) {
          fprintf(stderr,"vecgs size = %d, vecop empty!!!!\n", (int)vecgs.size());
		//countvecempty++;
		//if (countvecempty == 100)
			exit(1);
          /*while (!vecapair.empty()) {
              onepair e = vecapair.back();
              vecapair.pop_back();
              if (e.node1->GetMark() || e.node2->GetMark()) continue;
              priqueue->push(e);
          }
          continue;*/
      }
      //if (!vecop.empty())
          //fprintf(stderr,"vecop size = %d!!!!\n", (int)vecop.size());
      
      COneNode *node;
      oneoption op;
	   std::random_shuffle(vecop.begin(), vecop.end());
          op = vecop.back();
          vecop.pop_back();
	
          node = nodes[op.id];
      //fprintf(stderr,"node to recomGem!!\n");
          //node->View();

          const char* p = node->GetData();
          int l,r;
          node->GetLR(l,r);
          
          int brpoint = op.brpointgs;
          int countpoint = 0;
          int endpoint = 0;
          
          //fprintf(stderr, "gs.left = %d, gs.right = %d, brpoint = %d\n", vecgs[op.idvecgs].sleft, vecgs[op.idvecgs].sright, brpoint);
          
          if (brpoint-l >= distgam)
              countpoint = brpoint - distgam +1;
          else
              countpoint = l;
          
          if (r-brpoint > distgam)
              endpoint = brpoint + distgam;
          else
              endpoint = r;
          
          for (int i = countpoint; i <= brpoint; i++)
              for (int j = brpoint+1; j <= endpoint; j++)
              {
                  if (j-i > distgam)
                      break;
                  CheckGamete(p[i-l],p[j-l],i,j);
                  //CheckGamete(p[i-l],p[j-l],i*distgam+j-i);
              }
          
          
          RecomGemnode(node,brpoint,true);
          //fprintf(stderr, "done RecomGem---\n");
          sitemarked[vecgs[op.idvecgs].sleft*distgam+vecgs[op.idvecgs].sright-vecgs[op.idvecgs].sleft] = false;
    
        //for (int i = 0; i < (int)gemlist.size(); i++)
		//fprintf(stderr, "After RecomGem: type: %d, site: %d, mark = %d\n", gemlist[i].gamete, gemlist[i].site, (int)sitemarked[gemlist[i].site]);
	 //fprintf(stderr, "count00 = %d, count01 = %d, count10 = %d, count11 = %d, mark = %d\n",count00[gs.site],count01[gs.site],count10[gs.site],count11[gs.site], (int)sitemarked[gs.site]);
      while (!vecapair.empty()) {
          onepair e = vecapair.back();
          vecapair.pop_back();
          if (e.node1->GetMark() || e.node2->GetMark()) continue;
          priqueue->push(e);
      }
      
      //fprintf(stderr, "done push back priqueue---\n");
    for (int i = 0; i != (int) vecgs.size(); ++i)
        if ((sitemarked[vecgs[i].sleft*distgam+vecgs[i].sright-vecgs[i].sleft] == true))  prigemlist.push(vecgs[i]);

      //fprintf(stderr, "temp: type: %d, sleft: %d,, sright: %d, mark = %d\n", temp.gamete, temp.sleft, temp.sright, (int)sitemarked[temp.sleft*distgam+temp.sright-temp.sleft]);
      
    //if ((temp.gamete != -1) && (sitemarked[temp.sleft*distgam+temp.sright-temp.sleft] == true))  prigemlist.push(temp);
      //fprintf(stderr, "done insert to temp ---\n");
        return true;
      
}

bool CGAMARG :: IsNodegamete(COneNode *node, gemsite gs)
{
	int l,r;
	node->GetLR(l,r);
	if ((gs.sleft < l) || (gs.sright > r))
		return false;

    const char* p = node->GetData();
	//int l = node->GetLeft();
		//if(node->GetID() == 7860)
    		//	fprintf(stderr, "gs.gamete = %d, sleft = %c, sright = %c\n",gs.gamete, p[gs.sleft-l], p[gs.sright-l]);      			
    
    if ((gs.gamete == 0)&&(p[gs.sleft-l] == p[gs.sright-l])&&(p[gs.sleft-l] == '0')) {
        return true;
    }
    else if ((gs.gamete == 1)&&(p[gs.sleft-l] == '0')&&(p[gs.sright-l] == '1')) {
        return true;
    }
    else if ((gs.gamete == 2)&&(p[gs.sleft-l] == '1')&&(p[gs.sright-l] == '0')) {
        return true;
    }
    else if ((gs.gamete == 3)&&(p[gs.sleft-l] == p[gs.sright-l])&&(p[gs.sleft-l] == '1')) {
        return true;
    }
    else
        return false;
    //else {
    //    fprintf(stderr, "error in gamete: %d\t%d\n", gs.gamete, gs.site);
    //    exit(1);
    //}
}

void CGAMARG :: CheckGameteMut1st(char p, char p1, int idleft, int idright)
{
	/*const char* p = node->GetData();
     int point = 0;
     if (matchLeft == true)
     point = siteid;
     else
     point = siteid - 1;
     */
   
    int distlr = idright-idleft;
    int siteid = idleft*distgam+distlr;

    if ((p == '0')&&(p1 == '0')) {
        count00[siteid]++;
        count10[siteid]--;
        
        //if ((count00[siteid] == 0)||(count10[siteid] == 0)||(count01[siteid] == 0)||(count11[siteid] == 0))
            //sitemarked[siteid] = false;
        if ((count10[siteid] == 1)&&(count00[siteid] != 0)&&(count01[siteid] != 0)&&(count11[siteid] != 0)) {
            gemsite gs;
            gs.gamete = 2;
            gs.sleft = idleft;
            gs.sright = idright;
            //gs.site = siteid;
            prigemlist.push(gs);
            sitemarked[siteid] = true;
        }
    }
    else if ((p == '0')&&(p1 == '1')) {
        count01[siteid]++;
        count11[siteid]--;
        
        //if ((count00[siteid] == 0)||(count10[siteid] == 0)||(count01[siteid] == 0)||(count11[siteid] == 0))
            //sitemarked[siteid] = false;
        if ((count11[siteid] == 1)&&(count10[siteid] != 0)&&(count00[siteid] != 0)&&(count01[siteid] != 0)) {
            gemsite gs;
            gs.gamete = 3;
            gs.sleft = idleft;
            gs.sright = idright;
            //gs.site = siteid;
            prigemlist.push(gs);
            sitemarked[siteid] = true;
        }
    }
    else if ((p == '1')&&(p1 == '0')) {
        count10[siteid]++;
        count00[siteid]--;
        
        //if ((count00[siteid] == 0)||(count10[siteid] == 0)||(count01[siteid] == 0)||(count11[siteid] == 0))
            //sitemarked[siteid] = false;
        if ((count00[siteid] == 1)&&(count10[siteid] != 0)&&(count01[siteid] != 0)&&(count11[siteid] != 0)) {
            gemsite gs;
            gs.gamete = 0;
            gs.sleft = idleft;
            gs.sright = idright;
            //gs.site = siteid;
            prigemlist.push(gs);
            sitemarked[siteid] = true;
        }
    }
    else if ((p == '1')&&(p1 == '1')) {
        count11[siteid]++;
        count01[siteid]--;
        
        //if ((count00[siteid] == 0)||(count10[siteid] == 0)||(count01[siteid] == 0)||(count11[siteid] == 0))
            //sitemarked[siteid] = false;
        if ((count01[siteid] == 1)&&(count10[siteid] != 0)&&(count11[siteid] != 0)&&(count00[siteid] != 0)) {
            gemsite gs;
            gs.gamete = 1;
            gs.sleft = idleft;
            gs.sright = idright;
            //gs.site = siteid;
            prigemlist.push(gs);
            sitemarked[siteid] = true;
        }
    }
    else {
        fprintf(stderr,"Gem error in Mut1st!!! \n");
        exit(1);
    }
    if ((count00[siteid] == 0)||(count10[siteid] == 0)||(count01[siteid] == 0)||(count11[siteid] == 0))
        sitemarked[siteid] = false;
}

void CGAMARG :: CheckGameteMut2nd(char p, char p1, int idleft, int idright)
{
	/*const char* p = node->GetData();
     int point = 0;
     if (matchLeft == true)
     point = siteid;
     else
     point = siteid - 1;
     */
    //fprintf(stderr, "%c, %c \n", p, p1);
    //fprintf(stderr, "%c, %c, %d, %d \n", p, p1, idleft, idright);
    int distlr = idright-idleft;
    int siteid = idleft*distgam+distlr;

    if ((p == '0')&&(p1 == '0')) {
        count00[siteid]++;
        count01[siteid]--;
        
        if ((count01[siteid] == 1)&&(count10[siteid] != 0)&&(count00[siteid] != 0)&&(count11[siteid] != 0)) {
            gemsite gs;
            gs.gamete = 1;
            gs.sleft = idleft;
            gs.sright = idright;
            //gs.site = siteid;
            prigemlist.push(gs);
            sitemarked[siteid] = true;
        }
    }
    else if ((p == '0')&&(p1 == '1')) {
        count01[siteid]++;
        count00[siteid]--;
        
        if ((count00[siteid] == 1)&&(count10[siteid] != 0)&&(count01[siteid] != 0)&&(count11[siteid] != 0)) {
            gemsite gs;
            gs.gamete = 0;
            gs.sleft = idleft;
            gs.sright = idright;
            //gs.site = siteid;
            prigemlist.push(gs);
            sitemarked[siteid] = true;
        }
    }
    else if ((p == '1')&&(p1 == '0')) {
        count10[siteid]++;
        count11[siteid]--;
        
        if ((count11[siteid] == 1)&&(count00[siteid] != 0)&&(count01[siteid] != 0)&&(count10[siteid] != 0)) {
            gemsite gs;
            gs.gamete = 3;
            gs.sleft = idleft;
            gs.sright = idright;
            //gs.site = siteid;
            prigemlist.push(gs);
            sitemarked[siteid] = true;
        }
    }
    else if ((p == '1')&&(p1 == '1')) {
        count11[siteid]++;
        count10[siteid]--;
        
        if ((count10[siteid] == 1)&&(count11[siteid] != 0)&&(count01[siteid] != 0)&&(count00[siteid] != 0)) {
            gemsite gs;
            gs.gamete = 2;
            gs.sleft = idleft;
            gs.sright = idright;
            //gs.site = siteid;
            prigemlist.push(gs);
            sitemarked[siteid] = true;
        }
    }
    else {
        fprintf(stderr,"Gem error in Mut2nd!!! \n");
        exit(1);
    }
    if ((count00[siteid] == 0)||(count10[siteid] == 0)||(count01[siteid] == 0)||(count11[siteid] == 0))
        sitemarked[siteid] = false;

}


void CGAMARG :: CheckGamete(char p, char p1, int idleft, int idright)
{
	/*const char* p = node->GetData();
    int point = 0;
    if (matchLeft == true)
        point = siteid;
    else
        point = siteid - 1;
    */
    //fprintf(stderr, "%c, %c \n", p, p1);
    //fprintf(stderr, "Check gamete: %c, %c, %d, %d \n", p, p1, idleft, idright);
    int distlr = idright-idleft;
    
    int siteid = idleft*distgam+distlr;

    if ((p == '0')&&(p1 == '0')) {
        count00[siteid]--;
        
        //if ((count00[siteid] == 0)||(count10[siteid] == 0)||(count01[siteid] == 0)||(count11[siteid] == 0))
            //sitemarked[siteid] = false;
        if ((count00[siteid] == 1)&&(count10[siteid] != 0)&&(count01[siteid] != 0)&&(count11[siteid] != 0)) {
            gemsite gs;
            gs.gamete = 0;
            gs.sleft = idleft;
            gs.sright = idright;
            //gs.site = siteid;
            prigemlist.push(gs);
            sitemarked[siteid] = true;
        }
    }
    else if ((p == '0')&&(p1 == '1')) {
        count01[siteid]--;
        
        //if ((count00[siteid] == 0)||(count10[siteid] == 0)||(count01[siteid] == 0)||(count11[siteid] == 0))
            //sitemarked[siteid] = false;
        if ((count01[siteid] == 1)&&(count10[siteid] != 0)&&(count00[siteid] != 0)&&(count11[siteid] != 0)) {
            gemsite gs;
            gs.gamete = 1;
            gs.sleft = idleft;
            gs.sright = idright;
            //gs.site = siteid;
            prigemlist.push(gs);
            sitemarked[siteid] = true;
        }
    }
    else if ((p == '1')&&(p1 == '0')) {
        count10[siteid]--;
        
        //if ((count00[siteid] == 0)||(count10[siteid] == 0)||(count01[siteid] == 0)||(count11[siteid] == 0))
            //sitemarked[siteid] = false;
        if ((count10[siteid] == 1)&&(count00[siteid] != 0)&&(count01[siteid] != 0)&&(count11[siteid] != 0)) {
            gemsite gs;
            gs.gamete = 2;
            gs.sleft = idleft;
            gs.sright = idright;
            //gs.site = siteid;
            prigemlist.push(gs);
            sitemarked[siteid] = true;
        }
    }
    else if ((p == '1')&&(p1 == '1')) {
        count11[siteid]--;
        
        //if ((count00[siteid] == 0)||(count10[siteid] == 0)||(count01[siteid] == 0)||(count11[siteid] == 0))
            //sitemarked[siteid] = false;
        if ((count11[siteid] == 1)&&(count10[siteid] != 0)&&(count01[siteid] != 0)&&(count00[siteid] != 0)) {
            gemsite gs;
            gs.gamete = 3;
            gs.sleft = idleft;
            gs.sright = idright;
            //gs.site = siteid;
            prigemlist.push(gs);
            sitemarked[siteid] = true;
        }
    }
    else {
        fprintf(stderr,"Gem error!!! \n");
        exit(1);
    }
    if ((count00[siteid] == 0)||(count10[siteid] == 0)||(count01[siteid] == 0)||(count11[siteid] == 0))
        sitemarked[siteid] = false;
}

void CGAMARG :: CheckGameteCoal(char p, char p1, int idleft, int idright)
{
	/*const char* p = node->GetData();
     int point = 0;
     if (matchLeft == true)
     point = siteid;
     else
     point = siteid - 1;
     */
    //fprintf(stderr, "%c, %c \n", p, p1);
    //fprintf(stderr, "Check gamete: %c, %c, %d, %d \n", p, p1, idleft, idright);
    int distlr = idright-idleft;
    
    int siteid = idleft*distgam+distlr;

    if ((p == '0')&&(p1 == '0')) {
        count00[siteid]++;
        
        //if ((count00[siteid] == 0)||(count10[siteid] == 0)||(count01[siteid] == 0)||(count11[siteid] == 0))
        //sitemarked[siteid] = false;
        if ((sitemarked[siteid] == false)&&(count00[siteid] == 1)&&(count10[siteid] != 0)&&(count01[siteid] != 0)&&(count11[siteid] != 0)) {
            gemsite gs;
            gs.gamete = 0;
            gs.sleft = idleft;
            gs.sright = idright;
            //gs.site = siteid;
            prigemlist.push(gs);
            sitemarked[siteid] = true;
        }
    }
    else if ((p == '0')&&(p1 == '1')) {
        count01[siteid]++;
        
        //if ((count00[siteid] == 0)||(count10[siteid] == 0)||(count01[siteid] == 0)||(count11[siteid] == 0))
        //sitemarked[siteid] = false;
        if ((sitemarked[siteid] == false)&&(count01[siteid] == 1)&&(count10[siteid] != 0)&&(count00[siteid] != 0)&&(count11[siteid] != 0)) {
            gemsite gs;
            gs.gamete = 1;
            gs.sleft = idleft;
            gs.sright = idright;
            //gs.site = siteid;
            prigemlist.push(gs);
            sitemarked[siteid] = true;
        }
    }
    else if ((p == '1')&&(p1 == '0')) {
        count10[siteid]++;
        
        //if ((count00[siteid] == 0)||(count10[siteid] == 0)||(count01[siteid] == 0)||(count11[siteid] == 0))
        //sitemarked[siteid] = false;
        if ((sitemarked[siteid] == false)&&(count10[siteid] == 1)&&(count00[siteid] != 0)&&(count01[siteid] != 0)&&(count11[siteid] != 0)) {
            gemsite gs;
            gs.gamete = 2;
            gs.sleft = idleft;
            gs.sright = idright;
            //gs.site = siteid;
            prigemlist.push(gs);
            sitemarked[siteid] = true;
        }
    }
    else if ((p == '1')&&(p1 == '1')) {
        count11[siteid]++;
        
        //if ((count00[siteid] == 0)||(count10[siteid] == 0)||(count01[siteid] == 0)||(count11[siteid] == 0))
        //sitemarked[siteid] = false;
        if ((sitemarked[siteid] == false)&&(count11[siteid] == 1)&&(count10[siteid] != 0)&&(count01[siteid] != 0)&&(count00[siteid] != 0)) {
            gemsite gs;
            gs.gamete = 3;
            gs.sleft = idleft;
            gs.sright = idright;
            //gs.site = siteid;
            prigemlist.push(gs);
            sitemarked[siteid] = true;
        }
    }
    else {
        fprintf(stderr,"Gem error!!! \n");
        exit(1);
    }
    if ((count00[siteid] == 0)||(count10[siteid] == 0)||(count01[siteid] == 0)||(count11[siteid] == 0))
        sitemarked[siteid] = false;
}


bool CGAMARG :: Recom(void)
{


  if (priqueue->size() > 50000000  && priqueue->size() > activenodes.size() * activenodes.size() * 10)
    {
      prioqueue_onepair *temp = (prioqueue_onepair*) new prioqueue_onepair;
      while (!priqueue->empty())
	{
	  onepair e = priqueue->top();
	  priqueue->pop();
	  if (!e.node1->GetMark() && !e.node2->GetMark())  temp->push(e);
	}

      delete priqueue;
      priqueue = temp;
    }
    
    int len = 0;
    std::vector<onepair> qe;
    onepair temp;
    temp.match = 0;
    
    while (!priqueue->empty())
    {
        onepair e = priqueue->top();
        priqueue->pop();
        
        COneNode *node1 = e.node1, *node2 = e.node2;
        if (e.match == 0) continue;
        if (node1->GetMark() || node2->GetMark()) continue; // find a used node
        
        if (node1->GetLength() > node2->GetLength()) //break the smaller node
    	{
      		COneNode *t = node1; node1 = node2; node2 = t;
    	}
    	//const char *p = node1->GetData();
    	int countsite = 0;
    	int l,r;
    	node1->GetLR(l,r);
    	int start = 0;
    	int end = 0;
    	if(e.matchLeft == true) {
    		start = e.mismatchpoint+1;
    		end = r;
    	}
    	else {
    		start = l;
    		end = e.mismatchpoint-1;
    	}
    	for (int i = start; i<=end;i++) {
    			if(count0[i] == 0 || count1[i] == 0)
    				countsite++;
    			else
    				break;
    	}
    	//fprintf(stderr,"countsite = %d, e.match = %d \n", countsite, e.match);
    	if (countsite == e.match)
    		continue;
    		
        if (len == 0) len = e.match;
        if (e.match < len) { temp = e; break; }
        if (e.match == len) qe.push_back(e);
    }
    
    if (qe.size() == 0) return false;

    int i = rand()% qe.size();
    onepair e = qe[i];
    //    bool idnode = optionvec[i].isnode1;
    //fprintf(stderr,"qe size = %d, %d \n", (int)qe.size(), (int)idnode);
      //  if(idnode == false)
    	//	Recomlargernode(e.node2, e.mismatchpoint, e.matchLeft);
    	//else
	//if (e.node1->GetLength() > e.node2->GetLength()) //break the smaller node
    //	{
      //		COneNode *t = e.node1; e.node1 = e.node2; e.node2 = t;
    	//}

    
    //CheckGamete(p[point-l], p[point+1-l], point);
    Recom(e.node1, e.node2, e.mismatchpoint, e.matchLeft);
    	 
  nRec++;
  for (int i = 0; i != (int) qe.size(); ++i)
    if (!qe[i].node1->GetMark() && !qe[i].node2->GetMark())  priqueue->push(qe[i]);
  //for (int i = 0; i != (int) qe_temp.size(); ++i)
    //if (!qe_temp[i].node1->GetMark() && !qe_temp[i].node2->GetMark())  priqueue->push(qe_temp[i]);
  if ((temp.match!= 0) && !temp.node1->GetMark() && !temp.node2->GetMark())  priqueue->push(temp);
  return true;
}



void CGAMARG :: BuildARG(void)
{
  int oldsize = 0;
  while (true)
    {
	/*if((int)nodetothread.size()==0) {
            fprintf(stderr,"check!!!\n");
            CheckCount(); }*/
            //CheckCount();
      if ((int) nodes.size()   > oldsize  + 10000)
	{
	  fprintf(stderr, "size = %d, active size = %d, nsite = %d, mul = %d, priqueue = %d, coal = %d\n", (int) nodes.size(), (int) activenodes.size(), (int)  data->GetnSite(), (int) mutlist.size(), (int) priqueue->size(), (int)  coallist.size());
	  oldsize = (int) nodes.size();
      }

      /*
      for (int i = 0; i != (int) nodes.size(); ++i)
	if (nodes[i]->GetDataLength() > 0 && nodes[i]->GetLength() != nodes[i]->GetDataLength())
	  printf("%d\t", nodes[i]->GetLength(), nodes[i]->GetDataLength());
      */
       
      if (PossibleCoal()) continue;
      if (PossibleMutation()) continue;
         //fprintf(stderr, "size = %d, active size = %d, nsite = %d, mul = %d, priqueue = %d, coal = %d, gemlist = %d\n", (int) nodes.size(), (int) activenodes.size(), (int)  data->GetnSite(), (int) mutlist.size(), (int) priqueue->size(), (int)  coallist.size(), (int) prigemlist.size());
        if (PossibleGem()) continue;
      if (!Recom()) break;
    }
    //Check();
    	while((int)activenodes.size()>1) {
    		int id1 = activenodes[0];
		int r,l;
		nodes[id1]->GetLR(l,r);
		for (int i = 1; i < (int)activenodes.size(); i++) {
			int id2 = activenodes[i];
			int l2,r2;
			nodes[id2]->GetLR(l2,r2);
			if(l == r2+1 || r == l2-1) { 
			   Coal(nodes[id1],nodes[id2]);
			   break;
			}
		}
    		
		//nodes[id1]->View();
		//nodes[id2]->View();
      	}
    //gemlist.clear();
    //coallist.clear();
  //Check();
    
  fprintf(stderr, "\n**************************finish buiding ARG*************\n");
  fprintf(stderr, "Rec num %d\n", nRec);
    
    View();
}


void CGAMARG :: View(void) const
{
  for (int i = 0; i != (int) activenodes.size(); ++i)
    nodes[ activenodes[i] ]->View();
}


void CGAMARG :: Check(const COneNode *root, int siteid) const
{
  std::vector<const COneNode*> q;
  q.push_back(root);

  std::vector<int> getids;
  int i = 0;
  while (i < (int) q.size())
    {
      const COneNode *e = q[i++], *p = NULL;
      const COneNode *c1 = NULL, *c2 = NULL;

      while (1==1)
	{
	  e->GetChild(c1, c2);

	  int nc = 0;
	  if (c1 != NULL && c1->Contain(siteid))  { nc ++;  p = c1;}
	  if (c2 != NULL && c2->Contain(siteid))  { nc ++;  p = c2;}
	  if (nc == 1) e = p;
	  else break;
	}

      if (e->checkleaf()) {
	getids.push_back(e->GetID());
	continue;
      }
      if (c1 != NULL && c1->Contain(siteid))  q.push_back(c1);
      if (c2 != NULL && c2->Contain(siteid))  q.push_back(c2);
    }

  sort(getids.begin(), getids.end());

  bool error = false;
  int nhap = data->GetnHap();
  for (int i = 1; i < (int) getids.size(); ++i)
    if (getids[i] == getids[i-1]) error = true;

  if ((int)getids.size() != nhap || error)
    {
      fprintf(stderr, "getid size: %d\t nhap: %d \n", (int)getids.size(),nhap);
      fprintf(stderr, "error at site: %d\n", siteid);
      for (int i = 0; i != (int) getids.size(); ++i)
	fprintf(stderr, "%d\t", getids[i]);
      fprintf(stderr, "\n");
      exit(1);
    }
}

void CGAMARG :: Check(void) const
{
  COneNode *root = NULL;
  for (int i = 0; i != (int) nodes.size(); ++i)
    if (nodes[i]->checkroot())
      {
	root = nodes[i];
	break;
      }

  fprintf(stderr, "checking \n");
  for (int i = 0; i != data->GetnSite(); ++i) {
    if (i % 1000 == 0) printf("i = %d\n", i);
    Check(root, i);
  }
  fprintf(stderr, "ok in check\n");
}

void CGAMARG :: Write2ARG(FILE *fout)
{

  int nrm = 0;
  std::vector<const COneNode *> pp;
  pp.resize(4);
  for (int i = 0; i != (int) nodes.size(); ++i)
    if (nodes[i]->MutationNode())
      {

	fprintf(stderr, "\n>>>>>>> Mutation node %d \n", nodes[i]->GetID());
	nodes[i]->GetParents(pp[0], pp[1]);
      	nodes[i]->GetChild(pp[2], pp[3]);
      for (int j = 0; j != 4; ++j)
	if (pp[j] == NULL) fprintf(stderr, "\t-1");
	else fprintf(stderr, "\t%d", pp[j]->GetID());
	fprintf(stderr, "\n");

	nodes[i]->RemoveNode();
	delete nodes[i];
	nodes[i] = NULL;
	nrm ++;
      }

  fprintf(stderr, "rm node = %d --> %d\n", (int) nodes.size(), nrm);

  vecint poss;
  data->GetPos(poss);

  int l = 0;
  for (int i = 0; i != (int) nodes.size(); ++i)
    if (nodes[i] != NULL) nodes[l++] = nodes[i];
  nodes.resize(l);
  for (int i = 0; i != (int) nodes.size(); ++i)    nodes[i]->SetID(i);

  const vecstring *names = data->GetSampleNames();

  fprintf(fout, "sample %d\n", (int) names->size());
  fprintf(fout, "%s", names->at(0).c_str());
  for (int i = 1; i != (int) names->size(); ++i)
    fprintf(fout, "\t%s", names->at(i).c_str());




  fprintf(fout, "\n");
  std::vector<const COneNode *> p;
  p.resize(4);
  fprintf(fout, "nodes %d\n", (int) nodes.size());
  //fprintf(fout, "Rec num %d\n", nRec);
  //fprintf(fout, "\n -------- ARG --------\n");
  for (int i = 0;  i != (int) nodes.size(); ++i)
    {

      fprintf(fout, "%d\t%d\t%d", nodes[i]->GetID(), poss[nodes[i]->GetLeft()], poss[nodes[i]->GetRight()]);
      //      fprintf(fout, "%d\t%d\t%d", nodes[i]->GetID(), nodes[i]->GetLeft(), nodes[i]->GetRight());

      nodes[i]->GetParents(p[0], p[1]);
      nodes[i]->GetChild(p[2], p[3]);

      for (int j = 0; j != 4; ++j)
	if (p[j] == NULL) fprintf(fout, "\t-1");
	else fprintf(fout, "\t%d", p[j]->GetID());
      fprintf(fout, "\n");
    }

}
