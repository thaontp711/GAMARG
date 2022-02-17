#include "CData.h"


CData :: CData()
{

}


void CData :: ReadPos(const char *fn,   std::vector<int> &poss)
{
  igzstream inz (fn);
  int pos;
  std::string line;
  poss.clear();
  while (getline(inz, line))
    {
      sscanf(line.c_str(), "%d", &pos);
      poss.push_back(pos);
    }
  inz.close();
  sort(poss.begin(), poss.end());
}

void CData :: ReadrmPoss(const char* fn)
{
  ReadPos(fn, rmposs);
}

void CData :: ReadselectedPoss(const char* fn)
{
  ReadPos(fn, selectedposs);
}



void CData :: ReadNames(const char* fn, std::vector<std::string> &names)
{
  igzstream inz (fn);
  std::string line;
  names.clear();
  char tempc[1024];

  while (getline(inz, line)) {
    sscanf(line.c_str(), "%s", tempc);
    names.push_back(tempc);
  }
  inz.close();
  sort(names.begin(), names.end());
  //  for (int i = 0; i != names.size(); ++i) fprintf(stderr, "%s\n", names[i].c_str());
}


void CData :: Readrmsamples(const char* fn)
{
  ReadNames(fn, rmsamples);
}

void CData :: Readselectedsamples(const char* fn)
{
  ReadNames(fn, selectedsamples);
  fprintf(stderr, "selected samples = %d\n", (int) selectedsamples.size());
}

void CData :: ReadSampleFile(const char *fn)
{
  igzstream inz (fn);
  std::string line;

  //getline(inz, line);
  //getline(inz, line);
  //getline(inz, line);
  //sscanf(line.c_str(), "%d %d", &numcases, &numcontrols);
  //fprintf(stderr, "AAAAAAAAAAAAAAAAAAa %d %d\n", numcases, numcontrols);

  char tempc[1024];

  while (getline(inz, line))
    {
      sscanf(line.c_str(), "%s", tempc);
      samples.push_back(tempc);
    }
  inz.close();
}


void CData :: FilterSamples(void)
{

  if (rmsamples.size() == 0 && selectedsamples.size() == 0) return;
  nsp = 0;

  for (int i = 0; i != (int) samples.size(); ++i)
    {

      bool sel = true;
      bool rm = false;
      if (selectedsamples.size() > 0)
	sel = binary_search(selectedsamples.begin(), selectedsamples.end(), samples[i]);

      if (rmsamples.size() > 0)
	rm = binary_search(rmsamples.begin(), rmsamples.end(), samples[i]);

      if (sel && rm)
	{
	  fprintf(stderr, "sample %s exists in both the list of removed samples and selected samples\nCheck data!!!", samples[i].c_str());
	  exit(1);
	}

      if (sel && !rm)
	{
	  samples[nsp]  = samples[i];
	  haps[nsp*2]   = haps[i*2];
	  haps[nsp*2+1] = haps[i*2+1];
	  nsp ++;
	}
    }

  fprintf(stderr, "filter to %d samples and %d snps\n", nsp, (int) haps[0].size());

  samples.resize(nsp);
  haps.resize(2*nsp);

}






void CData :: ReadData(const char *fnhap, const char* fnsample)
{
  fprintf(stderr, "read data\n");

  ReadSampleFile(fnsample);
  nsp = samples.size();
  fprintf(stderr, "num of sample: %d\n", nsp);


  char tc1[1024], tc2[1024], ref[10], alt[10];
  int nhom = 0;

  igzstream inz (fnhap);
  std::string line;


  haps.resize(2*nsp, "");
  while (getline(inz, line))
    {
      snp e;
      const char *p = line.c_str();
      sscanf(p, "%s %s %d %s %s", tc1, tc2, &e.pos, ref, alt);
      //fprintf(stderr, "%s %s %d %s %s\n", tc1, tc2, &e.pos, ref, alt);
      e.rsid1 = tc1;
      e.rsid2 = tc2;
      e.al1 = ref[0];
      e.al2 = alt[0];

      if (rmposs.size()       != 0 &&  binary_search(rmposs.begin(), rmposs.end(), e.pos)) continue; //skip removed snps
      if (selectedposs.size() != 0 && !binary_search(selectedposs.begin(), selectedposs.end(), e.pos)) continue; //skip snps outside selectedposs

      Jump(p, "\t ", 5);

      int nc0 = 0, nc1 = 1;
      for (int i = 0; i != (int) strlen(p); ++i)
	{
	  if (p[i] == '0') nc0 ++;
	  if (p[i] == '1') nc1 ++;
	}

      if (nc1 == 0 || nc1 == 0)
	{
	  nhom ++;
	  continue;
	}

      snps.push_back(e);

      int sid = 0;
      while (p[0] != '\0')
	{
	  if (p[0] != ' ' && p[0] != '\t' && p[0] != '\n' && p[0] != '\0')
	    {
	      if (sid < 2*nsp)  haps[sid] += p[0];
	      sid++;
	    }
	  p++;
	}

      if (sid != 2*nsp)
	{
	  fprintf(stderr, "error in number of samples: %d\t%d\n", nsp, sid);
	  exit(1);
	}
    }
  inz.close();
  FilterSamples();
}


const vecstring* CData :: GetHaps(void) const
{
  return &(haps);
}


int CData :: GetnSite(void) const
{
  return haps[0].size();
}

int CData :: GetnHap(void) const
{
  return haps.size();
}



int CData :: GetnSample(void) const
{
  return haps.size()/2;
}


void CData :: GetPos(vecint &poss) const
{
  for (int i = 0; i !=(int) snps.size(); ++i) poss.push_back(snps[i].pos);
}

const vecstring * CData :: GetSampleNames(void) const
{
  return &samples;


}
