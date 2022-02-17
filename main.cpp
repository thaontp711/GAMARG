#include <stdio.h>
#include "CGAMARG.h"
#include "CARG.h"



void help(void)
{
  printf("./gamarg [options]\n");
  printf("  -excludedsnp <file of excluded snps>\n");
  printf("  -includeedsnp <file of included snps>\n");
  printf("  -excludedsample <file of excluded samples>\n");
  printf("  -includedsample <file of excluded samples>\n");
  printf("  -out <output ARG file>\n");
  printf("  -nthread <number of threads>\n");
  printf("  -narg <number of args>\n");
  printf("  -mapping <number of permutations>\n");
  printf("  -data <haps> <sample info>\n");
    printf("  -ndistgam <distance to check gamete types>\n");
  printf("  -trees <data info> <output marginal trees> \n");

}


int main(int ArgI, char* ArgC[])
{
  if (ArgI == 1) { help(); return 0;}
  CData *data = new CData();
  //CARG *arg = new CARG();
  
  std::string argfolder;
  std::string prefix;
  
  int nthread = 0;
  int ndistgam = 0;
  int narg = 0;
  bool printArg = false;
  bool printTrees = false;
  for (int i = 0; i != ArgI; ++i)
    {
      if (strcmp(ArgC[i], "-excludedsnp") == 0)  data->ReadrmPoss(ArgC[i+1]);
      if (strcmp(ArgC[i], "-includeedsnp") == 0) data->ReadselectedPoss(ArgC[i+1]);
      if (strcmp(ArgC[i], "-excludedsample")  == 0) data->Readrmsamples(ArgC[i+1]);
      if (strcmp(ArgC[i], "-includedsample")  == 0) data->Readselectedsamples(ArgC[i+1]);
      if (strcmp(ArgC[i], "-out") == 0) { printArg = true; argfolder = ArgC[i+1]; }
      if (strcmp(ArgC[i], "-narg") == 0) narg = atoi(ArgC[i+1]);
      if (strcmp(ArgC[i], "-nthread") == 0) nthread = atoi(ArgC[i+1]);
	  if (strcmp(ArgC[i], "-ndistgam") == 0) ndistgam = atoi(ArgC[i+1]);
	  if (strcmp(ArgC[i], "-trees") == 0) {printTrees = true; prefix = ArgC[i+1];}
    }
  
  for (int i = 0; i != ArgI; ++i)
    if (strcmp(ArgC[i], "-data") == 0)
      data->ReadData(ArgC[i+1],  ArgC[i+2]);
  
  for (int k = 0; k != narg; ++k)
    {
      fprintf(stderr, "ARG %d\n",k);
      CGAMARG *fm = new CGAMARG(nthread, ndistgam);
      fm->ReadData(data);
      fm->BuildARG();
      
      if(printArg == true) // print ARG
      {
		  char convertStrarg[10];
		  sprintf(convertStrarg,"%d",k);
		  char aa[260];
		  strcpy(aa,argfolder.c_str());
		  strcat(aa,convertStrarg);
	         FILE *farg = stdout;
		  farg = fopen(aa,"w");
	         fm->Write2ARG(farg);
		  fclose(farg);
      }
      
	  if(printTrees == true) // print marginal trees
	  {
	  	CARG *arg = new CARG(data);
		char convertStr[10];
		sprintf(convertStr, "%d", k);
		char a[260];
		strcpy(a, prefix.c_str());
		strcat(a, convertStr);
		FILE *ftrees = stdout;
		ftrees = fopen(a,"w");
		arg->nodes = fm->nodes;
		arg->getMarginalTrees(ftrees,k);
		fclose(ftrees);
		delete arg;
	  }
      delete fm;
    }
  delete data;
  //delete arg;
  return 0;
}

