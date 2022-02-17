#include "utilities.h"
void Jump(const char* &p, std::string s, int n)
{
  for (int i = 0; i != n; ++i)
    {
      while (p[0] != '\0' && strchr(s.c_str(), p[0]) != NULL) p++;
      while (p[0] != '\0' && strchr(s.c_str(), p[0]) == NULL) p++;
    }

  while (p[0] != '\0' && strchr(s.c_str(), p[0]) != NULL) p++;
}
