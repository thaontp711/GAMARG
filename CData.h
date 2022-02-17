
#ifndef __CDATA__
#define __CDATA__

#include <vector>
#include <string>
#include <algorithm>
#include "gzstream.h"
#include "utilities.h"
#include "define.h"


struct snp
{
  std::string rsid1, rsid2;
  int pos;
  char al1, al2;
};

class CData
{

 private:
  
  vecint rmposs, selectedposs;
  vecstring rmsamples, selectedsamples, allsamples, samples;
  /*
    rmsamples: list of  samples that must be removed
    selectedsamples: list of samples that includes
    allsamples: all samples read from files
    samples: list of samples after filtering
   */

  vecstring haps;

  void ReadPos(const char* fn, std::vector<int> &poss);
  void ReadNames(const char* fn, std::vector<std::string> &names);
  void FilterSamples(void);

  int nsp; //num sample

  std::vector<snp> snps;

 public:
  int numcases;
  int numcontrols;

  CData();

  void ReadrmPoss(const char*fn);
  void ReadselectedPoss(const char* fn);

  void Readrmsamples(const char*fn);
  void Readselectedsamples(const char*fn);

  void ReadSampleFile(const char *fn);
  void ReadData(const char* fnhap, const char* fnsample);

  const vecstring* GetHaps(void) const;

  int GetnSite(void) const;
  int GetnSample(void) const;
  int GetnHap(void) const;

  void GetPos(vecint &poss) const;

  const vecstring * GetSampleNames(void) const;

};


#endif
