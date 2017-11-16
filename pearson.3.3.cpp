// Pearson.3.3
// Charles Phillips
// Department of Electrical Engineering and Computer Science
// University of Tennessee
// 4-1-2012

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <sys/time.h>
#include <sys/stat.h>
#include <cstdlib>
#include <limits>
#include <string.h>

using namespace std;

// branch prediction macros
#define likely(x)	__builtin_expect(!!(x), 1)
#define unlikely(x)	__builtin_expect(!!(x), 0)

struct CorrVector
{
  string label;
  int size;
  int numMissing;
  double sum;
  double sumOfSquares;
  double squaredSum;
  double sumDivN;
  double stddev;
  double* values;
};

const string VERSION = "pearson3 v3.3";
vector <CorrVector> CORRVECTORS;
bool HISTOGRAM=false;
bool WEL=false;
bool EL=false;
bool NUMBERS=false;
bool POSITIVEONLY=false;
bool NEGATIVEONLY=false;
bool PVALUE=false;
bool SPEARMAN=false;
bool CORRELATION_MATRIX=false;
double THRESHOLD=999.0;
char* FILENAME=NULL;
int LastSize;
double* RVALUES;
string FirstRow;
int ZeroVarianceCount = 0;
CorrVector CV1, CV2;               //for use with missing values

void usage()
{
  cerr << "usage: pearson3 filename -flag1 -flag2 . . .\n\n";
  cerr << "\tfilename - Tab-separated file of expression values. The first row\n";
  cerr << "\t           should contain column labels and the first column should\n";
  cerr << "\t           contain labels for each row. The upper-left cell is\n";
  cerr << "\t           treated as a dummy label and can be blank. Any non-numeric\n";
  cerr << "\t           entry, including an empty string, is interpreted as a \n";
  cerr << "\t           missing value. Blank lines or lines with fewer than 2 values\n";
  cerr << "\t           are ignored\n";
  cerr << "\n\t  ----- flags -----\n";
  cerr << "\t-h - output a frequency histogram\n";
  cerr << "\t-w - output a weighted edge list\n";
  cerr << "\t-e - output an unweighted edge list\n";
  cerr << "\t-m - output a correlation matrix\n";
  cerr << "\t-s - use Spearman's correlation\n";
  cerr << "\t-n - output numbers for each vertex instead of the original labels\n";
  cerr << "\t-pos - only output edges with positive correlations\n";
  cerr << "\t-neg - only output edges with negative correlations\n";
  cerr << "\t-pvalue - use correlation p-values instead of correlations\n";
  cerr << "\t          (suitable for data with missing values)\n";
  cerr << "\n\tOne and only one of -h, -w, -e, or -m must be specified. If -w or -e is\n";
  cerr << "\tspecified, then an additional threshold parameter between 0 and 1 must\n";
  cerr << "\talso be specified. Any edges at or above the threshold will be retained\n";
  cerr << "\tin the output graph. If -pvalue is specified, then the two-tailed\n";
  cerr << "\tcorrelation p-value is used instead of the correlation; values at or\n";
  cerr << "\tbelow the threshold will be retained. In output p-values, 0 means <0.0001.\n";
  cerr << "\tRows with ultra-low-variance (< .000001) are assumed to have 0 correlation\n";
  cerr << "\twith all other rows.\n";
  exit(1);
}

//----------------------------------------------------------
// Returns the number of seconds since it was last called.
//----------------------------------------------------------
double elapsed()
{
  static double t1;
  timeval t;

  gettimeofday(&t, NULL);

  double t2 = t.tv_sec + (t.tv_usec/1000000.0);

  double elapsed = t2 - t1;
  t1 = t2;

  return elapsed;
}

//========================================
// Returns whether or not a file exists
//========================================
bool fileExists(const char* filename)
{
  struct stat buffer;
  if (stat(filename, &buffer))
    return false;
  return true;
}

//==================================================
// Removes any trailing whitespace from a string.
//==================================================
void chompBack(string& str)
{
  int i = str.length() - 1;

  while (i>=0 && (str[i] == '\n' || str[i] == '\t' || str[i] == '\r' || str[i] == ' ') )
    i--;

  str = str.substr(0, i+1);
}

//=======================================================================
// splits a null terminated character array s into a vector of strings
// c is the character to split on
//=======================================================================
void split(vector<string>& v, const char* s, char c)
{
  v.clear();
  while (true)
  {
    const char* begin = s;

    while (*s != c && *s) { ++s; }

    v.push_back(string(begin, s));

    if (!*s) break;

    if (!*++s) break;
  }

  if (*--s == c)                            //if last character is c, append another empty string
    v.push_back("");
}

//===================================================================
// alternate version of the above function
// does not append an empty string if the last character matches c
//===================================================================
//void split2(vector<string>& v, const char* s, char c)
//{
//  v.clear();
//  while (true)
//  {
//    const char* begin = s;

//    while (*s != c && *s) { ++s; }

//    v.push_back(string(begin, s));

//    if (!*s) break;

//    if (!*++s) break;
//  }
//}

//====================
// used by t_to_p()
//====================
double zip(double q, int i, int j)
{
  double zz = 1.0;
  double  z = 1.0;
  while(i <= j)
  {
      zz = zz * q * (i * 1.0) / ((i * 1.0) + 1);
      z = z + zz;
      i = i + 2;
  }
  return z;
}

//===================================
// converts t statistic to p-value
//===================================
double t_to_p(double t, int df) 
{
  t = fabs(t);
  static const double pj2 = M_PI / 2;

  double rt = t / sqrt(df * 1.0);
  double fk = atan(rt);

  if(df == 1)
    return 1 - fk / pj2;

  double dk = cos(fk);

  if(df % 2 == 1)
    return 1 - (fk + sin(fk) * dk * zip(dk*dk, 2, df-3)) / pj2;
  else
    return 1 - sin(fk) * zip(dk*dk, 1, df-3);
}

//==================================
// Fisher's r to z transformation
//==================================
inline double rtoz(double r) 
{
  return 0.5 * log((1.0 + r) / (1.0 - r));
}

//============================================
// normal cumulative probability function
// x is the value, mu is the average, sigma
// is the standard deviation
//============================================
inline double normal_cdf(double x, double mu, double sigma) 
{
  return (1.0 - 0.5 * (1.0 + erf((x - mu) / (sigma * sqrt(2.0)))));
}

// The following two functions were found online with the following note:
//
// File: CDF.cpp Christer Karlsson This code implements a function that calculates the standard
// normal CDF (x), using an approximation from Abromowitz and Stegun Handbook of Mathematical 
// Functions. http://www.math.sfu.ca/~cbm/aands/page_932.htm

//===================================================
// The Gaussian p.d.f with mean = 0 and stddev = 1 
//===================================================
double Z(const double x) 
{ 
  return (1/sqrt(2*M_PI))*exp(-x*x/2.0 ); 
} 

//=====================================================
// standard normal cumulative probability function
// about twice as fast as the normal_cdf() function,
// so use when the mean is 0 and the stddev is 1
//=====================================================
double CDF(const double x) 
{ 
  const double b1 = 0.319381530; 
  const double b2 = -0.356563782; 
  const double b3 = 1.781477937; 
  const double b4 = -1.821255978; 
  const double b5 = 1.330274429; 
  const double p = 0.2316419; 

  if(x >= 0.0) 
  { 
    double t = 1.0 / (1.0 + p*x); 
    return (1.0 - Z(x)*t* (t*(t*(t*(t*b5 + b4) + b3) + b2) + b1)); 
  }
  else 
  { 
    double t = 1.0 / ( 1.0 - p*x ); 
    return ( Z(x)*t* (t*(t*(t*(t*b5 + b4) + b3) + b2) + b1));
  }
} 


//==============================================
// Returns the pvalue for a given correlation
// and number of data points
//==============================================
double r_to_p(double r, int n)
{
  if (n < 4)                        // only calculate for 4 or more data points
    return 1.0;

  if (fabs(r) < .0000001)
    return 1.0;

  double t = r * sqrt((n-2) / (1 - r * r));        // t statistic
  double p = t_to_p(t, n-2);
  //cout << "r: " << r << "\tn: " << n << "\tt: " << t << "\tp: " << p << endl; 
  return p;

//old method
// double y, z, zt;
// z = rtoz(fabs(r));
// zt = z / sqrt(1.0 / (n-3));            //divide by the standard error
// y = normal_cdf(zt, 0, 1);
// return 2*y;                  //multiply for two-tailed test
}

//===============================================
// Returns the correlation for a given p-value
// and number of data points
//===============================================
double p_to_r(double p, int n)
{
  double p2;                
  double r = .5;
  double adder = .25; 

  while (true)
  {
    p2 = r_to_p(r, n);

    if (fabs(p2 - p) < .00000000000001)
      break;

    if (p2 < p)
      r -= adder;
    else
      r += adder;

    adder /= 2.0;        
  }      
  return r;
}

//=========================================================================
// returns the index of the correlation value; for calculating histogram
//=========================================================================
int hash(double x)
{
  x *= 100;

  if (x < 0)
    x += 100.5;
  else
    x += 100.5000000001;

  return (int)x;
}

//================================================================
// Check the command-line arguments and set up the global flags
//================================================================
void parseArgs(int argc, char** argv)
{
  if (argc < 3 || argc > 6) usage();

  for (int i=1; i<argc; i++)
  {
    //must be one and only one file name
    if ( argv[i][0] != '-'
         && (argv[i][0] != '.' || (argv[i][0] == '.' && (argv[i][1] == '.' || argv[i][1] == '/')))
         && ( (argv[i][0] != '0' && argv[i][0] != '1') || FILENAME == NULL))
    {
      
      if (FILENAME == NULL)
        FILENAME = argv[i];
      else
      {
        cerr << "can only specify one file name or file name must come before threshold \n";
        exit(1);
      }
    }
    else                                           //parse flags
    {
      if (!strcmp(argv[i], "-m"))
      {
        CORRELATION_MATRIX = true;
      }
      else if (!strcmp(argv[i], "-s"))
      {
        SPEARMAN = true;
      }    
      else if (!strcmp(argv[i], "-h"))
      {
        if (THRESHOLD != 999)
        {
          cerr << "cannot specify threshold with -h flag\n";
          exit(1);
        }
        if (POSITIVEONLY)
        {
          cerr << "cannot specify both -h and -pos flags\n";
          exit(1);
        }
        if (NEGATIVEONLY)
        {
          cerr << "cannot specify both -h and -neg flags\n";
          exit(1);
        }

        if (!HISTOGRAM && !WEL && !EL)
          HISTOGRAM = true;
        else
        {
          cerr << "Can only specify one of -h, -w, or -e\n";
          exit(1);
        }
      }
      else if (!strcmp(argv[i], "-w"))
      {
        if (!HISTOGRAM && !WEL && !EL)
          WEL = true;
        else
        {
          cerr << "Can only specify one of -h, -w, or -e\n";
          exit(1);
        }
      }
      else if (!strcmp(argv[i], "-e"))
      {
        if (!HISTOGRAM && !WEL && !EL)
          EL = true;
        else
        {
          cerr << "Can only specify one of -h, -w, or -e\n";
          exit(1);
        }
      }
      else if (!strcmp(argv[i], "-n"))
      {
        if (!NUMBERS)
          NUMBERS = true;
        else
          usage();
      }
      else if (!strcmp(argv[i], "-pos"))
      {
        if (HISTOGRAM)
        {
          cerr << "cannot specify both -h and -pos flags\n";
          exit(1);
        }
        if (NEGATIVEONLY)
        {
          cerr << "cannot specify both -pos and -neg flags\n";
          exit(1);
        }
        if (POSITIVEONLY)
        {
          cerr << "cannot specify -pos flag more than once\n";
          exit(1);
        }
        POSITIVEONLY = true;
      }
      else if (!strcmp(argv[i], "-neg"))
      {
        if (HISTOGRAM)
        {
          cerr << "cannot specify both -h and -neg flags\n";
          exit(1);
        }
        if (POSITIVEONLY)
        {
          cerr << "cannot specify both -pos and -neg flags\n";
          exit(1);
        }
        if (NEGATIVEONLY)
        {
          cerr << "cannot specify -neg flag more than once\n";
          exit(1);
        }
        NEGATIVEONLY = true;
      }
      else if ((argv[i][0] == '.' && argv[i][1] != '/' && argv[i][1] != '.')
               || argv[i][0] == '0' || argv[i][0] == '1')
      {
        if (THRESHOLD != 999)
        {
          cerr << "can only specify one threshold\n";
          exit(1);
        }
        if (!HISTOGRAM)
        {
          THRESHOLD = atof(argv[i]);
          if (THRESHOLD > 1.0)
          {
            cerr << "threshold must be between 0 and 1\n";
            exit(1);
          }
        }
        else
        {
          cerr << "cannot specify threshold with -h flag\n";
          exit(1);
        }
      }
      else if (!strcmp(argv[i], "-pvalue"))
      {
        if (PVALUE)
        {
          cerr << "-pvalue flag can only be specified once\n";
          exit(1);
        }
        PVALUE = true;
      }
      else if (argv[i][0] == '-')
      {
        cerr << "unrecognized flag: " << argv[i] << endl;
        exit(1);
      }
      else
        usage();
    }
  }

  if (FILENAME == NULL)
  {
    cerr << "file name not given\n";
    exit(1);
  }

  if (WEL && THRESHOLD == 999)
  {
    cerr << "threshold needed for -w flag\n";
    exit(1);
  }
  if (EL && THRESHOLD == 999)
  {
    cerr << "threshold needed for -e flag\n";
    exit(1);
  }
  if (!HISTOGRAM && !WEL && !EL && !CORRELATION_MATRIX)
  {
    cerr << "must specify one of -h, -e, -w, or -m flags\n";
    exit(1);
  }
  if (PVALUE && SPEARMAN)
  {
    cerr << "cannot specify both -p and -s\n";
    exit(1);
  }

//print the flags
//  cerr << "filename: " << FILENAME << endl;
//  cerr << "HISTOGRAM: " << HISTOGRAM << endl;
//  cerr << "WEL: " << WEL << endl;
//  cerr << "EL: " << EL << endl;
//  cerr << "NUMBERS: " << NUMBERS << endl;
//  cerr << "POSITIVEONLY: " << POSITIVEONLY << endl;
//  cerr << "NEGATIVEONLY: " << NEGATIVEONLY << endl;
//  cerr << "THRESHOLD: " << THRESHOLD << endl;
}

//===========================================================
// Reads the input file and sets up the global CORRVECTORS
//===========================================================
void readInputFile(char* FILENAME)
{
  int i, previoussize, numMissing=0;
  ifstream f;
  string str;
  vector <string> columnLabels;
  vector <string> line;
  CorrVector corrvec;
  bool unequalWarning = false;

  f.open(FILENAME);
  if (!f) { cerr << "Error opening input file '" << FILENAME << "'\n"; exit(1); }

  getline(f, str);   	                         //read first row
  FirstRow = str;

  chompBack(str);
  split(columnLabels, str.c_str(), '\t');        //find the number of columns
  cerr << columnLabels.size() << " columns, ";

  int rows=0;
  while (getline(f, str))                        //read in each row
  {
    split(line, str.c_str(), '\t');

    if (line.size() < 3 || line[0] == "")        //ignore blank lines, lines with
      continue;                                  //fewer than 2 values, or blank labels

    if (rows>1 && line.size() != previoussize)
      unequalWarning = true;
    previoussize = line.size();

    corrvec.label = line[0];                     //set up the correlation vector
    corrvec.size = line.size()-1;
    corrvec.sum = corrvec.sumOfSquares = corrvec.numMissing = 0;
    corrvec.values = new double[corrvec.size];

    bool hasMissingValue = false;
    for (i=1; i<line.size(); i++)                //parse each value
    {
      corrvec.values[i-1] = atof(line[i].c_str());

      //check for missing value
      if (corrvec.values[i-1] == 0.0 && line[i] != "0" && line[i] != "0.0")
      {
        corrvec.values[i-1] = numeric_limits<double>::max();
        corrvec.numMissing++;
        hasMissingValue = true;
      }
    }

    if (hasMissingValue)
      numMissing++;

    //cerr << "line has " << corrvec.size << " elements\n";

    CORRVECTORS.push_back(corrvec);
    rows++;
  }

  f.close();

  cerr << rows << " rows\n";
  cerr << numMissing << " rows contain missing values\n";

  if (unequalWarning)
    cerr << "NOTE: not all rows contain the same number of values\n";
}

//==============================================================
// convert a CorrVector to ranks in preparation for Spearman's
//==============================================================
void convert_to_ranks(CorrVector& X)
{
  int i, j, numidentical;
  double rank;
  vector <double> vec;
  vector <double> vecsave;

//  for (i=0; i<X.size; i++)
//    cout << X.values[i] << "\t";
//  cout << endl;

  for (i=0; i<X.size; i++)
  {
    vec.push_back(X.values[i]);
    vecsave.push_back(X.values[i]);
  }

  sort(vec.begin(), vec.end());
  sort(vecsave.begin(), vecsave.end());
  
//  for (i=0; i<vec.size(); i++)
//    cout << vec[i] << "\t";
//  cout << endl;

  for (i=0; i<vec.size(); i++)
  {
    numidentical = 1;
    while (i < vec.size() && vec[i+1] == vec[i])
    {
      i++;
      numidentical++;
    }

    rank = 0.0;
    for (j=i-numidentical+1; j<=i; j++)
      rank += j+1;
    rank /= numidentical;

    for (j=i-numidentical+1; j<=i; j++)
      vec[j] = rank;
  }

//  for (i=0; i<vec.size(); i++)
//    cout << vec[i] << "\t";
//  cout << endl;

  for (i=0; i<X.size; i++)
    for (j=0; j<vecsave.size(); j++)
      if (X.values[i] == vecsave[j])
      {
        X.values[i] = vec[j];
        break;
      }

//  for (i=0; i<X.size; i++)
//    cout << X.values[i] << "\t";
//  cout << endl << endl;

}


//===============================================
// precompute values in the global CORRVECTORS
//===============================================
void precomputeValues()
{
  int i, j, largest=0;

  if (SPEARMAN)                    //convert to ranks for Spearman's rank correlation
    for (i=0; i<CORRVECTORS.size(); i++)
      if (!CORRVECTORS[i].numMissing > 0)
        convert_to_ranks(CORRVECTORS[i]);

  for (i=0; i<CORRVECTORS.size(); i++)
  {
    if (CORRVECTORS[i].size > largest)
      largest = CORRVECTORS[i].size;

    for (j=0; j<CORRVECTORS[i].size; j++)
    {
      CORRVECTORS[i].sum += CORRVECTORS[i].values[j];
      CORRVECTORS[i].sumOfSquares += CORRVECTORS[i].values[j] * CORRVECTORS[i].values[j];
    }
    CORRVECTORS[i].squaredSum = CORRVECTORS[i].sum * CORRVECTORS[i].sum;
    CORRVECTORS[i].sumDivN = CORRVECTORS[i].sum / CORRVECTORS[i].size;

    CORRVECTORS[i].stddev = CORRVECTORS[i].sumOfSquares - CORRVECTORS[i].squaredSum / CORRVECTORS[i].size;

    if (CORRVECTORS[i].stddev <= 0.000001)           //deal with floating-point rounding error
      CORRVECTORS[i].stddev = 0.0;                   //due to ultra-low-variance in some probes
    else
      CORRVECTORS[i].stddev = sqrt(CORRVECTORS[i].stddev);
  }

  CV1.values = (double*)malloc(largest * sizeof(double));            //allocate space for the two global CorrVectors
  CV2.values = (double*)malloc(largest * sizeof(double));

  if (PVALUE)
  {
    int largest = 0;                           //find largest vector
    for (i=0; i<CORRVECTORS.size(); i++)
      if (CORRVECTORS[i].size > largest)
        largest = CORRVECTORS[i].size;

    RVALUES = new double[largest+1];

    for (i=4; i<=largest; i++)
      RVALUES[i] = p_to_r(THRESHOLD, i);
  }

}

//================================
// print the global CORRVECTORS
//================================
void printCORRVECTORS()
{
  for (int i=0; i<CORRVECTORS.size(); i++)
  {
    for (int j=0; j<CORRVECTORS[i].size; j++)
    {
      if (CORRVECTORS[i].values[j] == numeric_limits<double>::max())
        cerr << ".\t";
      else
        cerr << CORRVECTORS[i].values[j] << "\t";
    }
    cerr << endl;
  }
}

//=============================
// print a single CorrVector
//=============================
void printCorrVector(const CorrVector& corrvector)
{
  for (int i=0; i<corrvector.size; i++)
    cerr << corrvector.values[i] << "\t";
  cerr << endl;
}

//=============================================================
// Computes the Pearson correlation between two CorrVectors
// when the vectors have the same size and no missing values
//=============================================================
inline double correlationNoMissing(const CorrVector& X, const CorrVector& Y)
{
  int i;
  double r, numerator=0;

//  for (i=0; i<X.size; i++)
//    numerator += X.values[i] * Y.values[i];

  int size4 = X.size - X.size % 4;                //unroll the two-line loop
                                                  //commented out above
  for (i=0; i<size4; i+=4)
  {
    numerator += X.values[i] * Y.values[i];
    numerator += X.values[i+1] * Y.values[i+1];
    numerator += X.values[i+2] * Y.values[i+2];
    numerator += X.values[i+3] * Y.values[i+3];
  }

  for (i=size4; i<X.size; i++)
    numerator += X.values[i] * Y.values[i];       //end loop unroll

  numerator -= X.sum * Y.sumDivN;

  double denominator = X.stddev * Y.stddev;

  if (denominator != 0.0)                         //if either vector has zero variance
    r = numerator / denominator;                  //return zero for the correlation
  else
    r = 0.0;

  if (r > 1.0 && r < 1.1)
    r = 1.0;
  if (r < -1.0 && r > -1.1)
    r = -1.0;

  if (r < -1.0 || r > 1.0) {
    cerr << "INTERNAL ERROR: correlation out of bounds: " << r << endl;

    for (i=0; i<X.size; i++)
      cerr << X.values[i] << "\t";
    cerr << endl;
    cerr << "stddev = " << X.stddev << endl;

    for (i=0; i<Y.size; i++)
      cerr << Y.values[i] << "\t";
    cerr << endl;    
    cerr << "stddev = " << Y.stddev << endl;

    cerr << "stddev X: " << X.stddev << endl;
    cerr << "stddev Y: " << Y.stddev << endl;
    cerr << "numerator: " << numerator << endl;
    cerr << "denominator: " << denominator << endl;
    
//    double devnumerX = X.sumOfSquares - X.squaredSum/X.size;
//    double devnumerY = Y.sumOfSquares - Y.squaredSum/Y.size;

//    CORRVECTORS[i].stddev = sqrt(CORRVECTORS[i].sumOfSquares - CORRVECTORS[i].squaredSum / CORRVECTORS[i].size);

//    cerr << "devnumerX = " << devnumerX << endl;
//    cerr << "devnumerY = " << devnumerY << endl;

    exit(1);
  }

  LastSize = X.size;                              //in case PVALUE flag is set
  return r;
}



//============================================================
// Computes the Pearson correlation between two CorrVectors
// when the vectors have different sizes or missing values
//============================================================
double correlationMissing(const CorrVector X, const CorrVector Y)
{
  int i, j, size=0, smallest;
  double r, sumxy=0, numerator=0;

  X.size <= Y.size ? smallest = X.size : smallest = Y.size;    //get size of smaller vector
  CV1.label = X.label;
  CV2.label = Y.label;
  CV1.sum = CV2.sum = 0;
  CV1.sumOfSquares = CV2.sumOfSquares = 0;

  for (i=0, j=0; i<smallest; i++)       //set up the two vectors
  { 
    if (X.values[i] != numeric_limits<double>::max() && Y.values[i] != numeric_limits<double>::max())
    {
      CV1.values[j] = X.values[i];
      CV2.values[j] = Y.values[i];
      j++;
    }
  }
  CV1.size = CV2.size = j;

  if (SPEARMAN)
  {
    convert_to_ranks(CV1);
    convert_to_ranks(CV2);

//    cerr << CV1.label << "\t";
//    for (i=0; i<CV1.size; i++)
//      cerr << CV1.values[i] << "\t";
//    cerr << endl;
//    cerr << CV2.label << "\t";
//    for (i=0; i<CV2.size; i++)
//      cerr << CV2.values[i] << "\t";
//    cerr << endl << endl;
  }

  for (i=0; i<CV1.size; i++)       //set up the two vectors
  {
      CV1.sum += CV1.values[i];
      CV2.sum += CV2.values[i];

      CV1.sumOfSquares += CV1.values[i] * CV1.values[i];
      CV2.sumOfSquares += CV2.values[i] * CV2.values[i];

      numerator += CV1.values[i] * CV2.values[i];

      size++;
  }

  CV1.stddev = sqrt(CV1.sumOfSquares - (CV1.sum * CV1.sum) / size);
  CV2.stddev = sqrt(CV2.sumOfSquares - (CV2.sum * CV2.sum) / size);

  numerator -= CV1.sum * CV2.sum / size;

  double denominator = CV1.stddev * CV2.stddev;

  if (denominator != 0 && size > 2)               //if either vector has zero variance or
    r = numerator / denominator;                  //the vectors have fewer than 3 points in 
  else                                            //common return zero for the correlation
    r = 0;

//  cerr << "sum: " << A.sum << "\t" << B.sum << endl;
//  cerr << "sumOfSquares: " << A.sumOfSquares << "\t" << B.sumOfSquares << endl;
//  cerr << "squaredSum: " << A.squaredSum << "\t" << B.squaredSum << endl;
//  cerr << "sumDivN: " << A.sumDivN << "\t" << B.sumDivN << endl;
//  cerr << "stddev: " << A.stddev << "\t" << B.stddev << endl;
//  cerr << "numerator: " << numerator << endl;
//  cerr << "denominator: " << denominator << endl;
//  cerr << "correlation: " << r << endl;

  if (r > 1.0 && r < 1.1)                        //deal with floating-point rounding errors
    r = 1.0;
  if (r < -1.0 && r > -1.1)
    r = -1.0;

  if (r < -1.0 || r > 1.0) { 
    cerr << "INTERNAL ERROR: correlation out of bounds: " << r << endl;
    exit(1);
  }

  LastSize = size;                               //in case PVALUE flag is set
  return r;
}


//===================================================
// returns the correlation between two CorrVectors
//===================================================
double correlation(CorrVector& X, const CorrVector& Y) 
{
  double r;

  if (likely(X.numMissing == 0 && Y.numMissing == 0 && X.size == Y.size))
    r = correlationNoMissing(X, Y);
  else
    r = correlationMissing(X, Y);

  if (r == r)             // check if correlation is NaN 
    return r;
  else
  {
    ZeroVarianceCount++;
    return 0.0;
  }
}

//==========================================
// output a histogram of the correlations
//==========================================
void OutputHistogram()
{
  int i, j;
  cerr << "outputting histogram\n";

  unsigned long long* counts = (unsigned long long*)calloc(201, sizeof(unsigned long long));   //64 bit counter

  for (i=0; i<CORRVECTORS.size(); i++)
  {
    for (j=i+1; j<CORRVECTORS.size(); j++)
    {
      double r = correlation(CORRVECTORS[i], CORRVECTORS[j]);
      if (r != r) 
      {
        cerr << "Error in OutputHistogram(), r is NAN!!!\n";
        exit(1);
      }
      if (r < -1.0 || r > 1.0) { cerr << "ERROR: r is " << r << endl; exit(1); }
      counts[hash(r)]++;
      //cerr << hash(r) << endl;
    }
  }

  for (i=0; i<201; i++)
    cout << ((double)i-100.0)/100.0 << "\t" << counts[i] << endl;

  cerr << ZeroVarianceCount << " correlations were set to 0.0 because one correlate had zero variance\n";

  free(counts);
}

//================================
// output a correlation matrix
//================================
void OutputCorrelationMatrix()
{
  int i, j;
  double r;
  cerr << "outputting correlation matrix\n";

  for (i=0; i<CORRVECTORS.size(); i++)
    if (NUMBERS)
      cout << "\t" << i;
    else
      cout << "\t" << CORRVECTORS[i].label;
  cout << endl;

  for (i=0; i<CORRVECTORS.size(); i++)
  {
    if (NUMBERS)
      cout << i;
    else
      cout << CORRVECTORS[i].label;

    for (j=0; j<CORRVECTORS.size(); j++)
    {
      if (j != i)
        r = correlation(CORRVECTORS[i], CORRVECTORS[j]);
      else
        r = 1.0;

      r = (int((r + .0005) * 1000.0)) / 1000.0;   //round to 3 decimal places

      cout << "\t" << r;
    }
    cout << endl;
  }
}


//=========================================================
// finds the first available file name for the temp file
//=========================================================
string GetTempFilename()
{
  int i;
  string filename;
  char fnum[2];

  for (i=1; i<=99; i++)            //find first free temp file name
  {
    sprintf(fnum, "%d", i);
    filename = "el";       
    filename += fnum;
    filename += ".temp";
    if (!fileExists(filename.c_str()))
       break;
  }
  if (i >= 99) { cerr << "more than 99  el?.temp files exist\n"; exit(1); }
  return filename;
}

//=================================================================
// Writes the contents of a file to STDOUT then deletes the file
//=================================================================
void WriteTempFileToSTDOUT(string filename)
{
  ifstream f;
  string str;

  f.open(filename.c_str(), ios::in);
  if (!f) { cerr << "Error opening file '" << filename << "' for reading\n"; exit(1); }

  while(getline(f, str))      
    cout << str << endl;

  f.close(); 
  str = "/bin/rm " + filename;               
  int xyz = system(str.c_str());
}

//===================================================
// Output either a weighted or an unweighted graph
//===================================================
void OutputGraph()
{
  int i, j, numvertices=0;
  unsigned long long numedges=0;         //64 bit
  double p, x;
  unsigned long long numcorrelations=0;
  string tempfilename;
  char fnum[2];
  ofstream temp;
  bool* vertexpresent = new bool[CORRVECTORS.size()];

  for (i=0; i<CORRVECTORS.size(); i++)
    vertexpresent[i] = false;

  tempfilename = GetTempFilename();          //find first free temp file name

  temp.open(tempfilename.c_str(), ios::out);        //open the temp file
  if (!temp) { cerr << "Error opening file '" << tempfilename << "' for writing"; exit(1); }
 
  for (i=0; i<CORRVECTORS.size(); i++)          //write each correlation >= threshold
    for (j=i+1; j<CORRVECTORS.size(); j++)      //to the temp file
    {
      numcorrelations++;
      double r = correlation(CORRVECTORS[i], CORRVECTORS[j]);

      if ((!PVALUE && fabs(r) >= THRESHOLD) || (PVALUE && fabs(r) >= RVALUES[LastSize]))
        if ( (!POSITIVEONLY && !NEGATIVEONLY) || (POSITIVEONLY && r>=0) || (NEGATIVEONLY && r<=0) )
        {
          if (EL)
            temp << CORRVECTORS[i].label << "\t" << CORRVECTORS[j].label << endl;
          else if (WEL)
          {
            if (PVALUE)
            {
              p = r_to_p(r, LastSize);
              x = (int((p + .00005) * 10000.0)) / 10000.0;   //round to 4 decimal places
            }
            else
              x = (int(r * 10000.0)) / 10000.0;              //truncate to 4 decimal places

            temp << CORRVECTORS[i].label << "\t" << CORRVECTORS[j].label << "\t" << x << endl;
          }
          else
          { 
             cerr << "INTERNAL ERROR: neither EL nor WEL specified in OutputGraph() function\n";
             exit(1);
          }
          vertexpresent[i] = vertexpresent[j] = true;
          numedges++;
        }
    }

  temp.close();

  for (i=0; i<CORRVECTORS.size(); i++)    //count how many vertices are actually in the graph
    if (vertexpresent[i])
      numvertices++;

  delete [] vertexpresent;

  cerr << numedges << " of " << numcorrelations << " correlations (" << (numedges / (numcorrelations * 1.0)) * 100.0 <<  "%) were at or above the threshold\n";
  cerr << "output graph has " << numvertices << " vertices, " << numedges << " edges\n"; 

  cout << numvertices << "\t" << numedges << endl;
  WriteTempFileToSTDOUT(tempfilename);
}

int main(int argc, char** argv)
{
  cerr << VERSION << endl;;

  elapsed();

  parseArgs(argc, argv);       //parse the arguments and set globals

  readInputFile(FILENAME);     //read the file and set up the CORRVECTOR
  precomputeValues();          //perform all possible precomputations

//  printCORRVECTORS();

  if (HISTOGRAM)
    OutputHistogram();
  else if (WEL || EL)
    OutputGraph();
  else if (CORRELATION_MATRIX)
    OutputCorrelationMatrix();
  else
  {
    cerr << "Internal error: None of HISTOGRAM, WEL, or EL was specified\n"; 
    exit(1);
  }

  cerr << "elapsed time: " << elapsed() << " sec.\n";
}
