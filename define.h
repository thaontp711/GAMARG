#ifndef __DEFINE__
#define __DEFINE__



#define MIN(x, y)  ((x) > (y) ? (y) : (x))
#define MAX(x, y)  ((x) < (y) ? (y) : (x))
#define NOTMISSING(x) (((x) == '0' || (x) == '1') ? true:false)
#define COPY true
#define NOCOPY false

#define M20 20000000

#define POSSIBLEMUT(x, y)  ( ((x) + (y) > 1 && ((x) == 1 || (y) == 1)) ? true:false)
typedef std::vector<std::string> vecstring;
typedef std::vector<int> vecint;
typedef std::vector<bool> vecbool;
#endif
