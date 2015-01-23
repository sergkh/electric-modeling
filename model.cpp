#include <iostream>
#include <math.h>
#include <sstream>
#include <fstream>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <iomanip>

using namespace std;

bool debug = false;

struct InParams {
  double p0;
  double Vc;
  double E;
  double a;
  double n;

  InParams() {
    p0 = -1;
    Vc = 0.01;
    E = 1;
    a = 1;
    n = 2;
  }
};

struct Results {
  double t;
  double* A;
  double* p;
  double* delta;
};

bool hasOption(int argc, char* argv[], const char* option) {
  for (int i = 1; i < argc; i++) 
  {
    if (!strcmp(argv[i], option)) 
    {
      return true;
    }
  }

  return false;
}

char* getOption(int argc, char* argv[], const char* option) {
  for (int i = 1; i < argc; i++) 
  {
    if (!strcmp(argv[i], option)) 
    {
      return argv[i + 1];
    }
  }

  return NULL;  
}

int readParams(char *fileName, InParams **params) {
  if(debug) clog << "Reading params from file: " << fileName << endl;

  int lines = 0;
  ifstream infileLines(fileName);

  if(!infileLines.is_open()) 
  {
    if(debug) clog << "File is empty" << endl;		
    return 0;
  }

  string line;
  while (std::getline(infileLines, line)) if(line.length() > 0) lines++;
  infileLines.close();

  if(lines == 0) return 0;

  (*params) = new InParams[lines];
  InParams *parArray = (*params);
  int curIdx = 0;

  ifstream infile(fileName);
  while (getline(infile, line)) 
  {
    if(line.length() == 0) continue;

    // find params on the line
    istringstream iss(line);
    
    while(iss) 
    {
      string st;
      iss >> st;

      if(st.length() == 0) continue;

      if(!st.compare(0, 3, "p0=")) 
      {
        parArray[curIdx].p0 = atof(st.substr(3, st.length() - 2).c_str());
      }

      if(!st.compare(0, 3, "Vc=")) 
      {
        parArray[curIdx].Vc = atof(st.substr(3, st.length() - 2).c_str());
      }

      if(!st.compare(0, 2, "E=")) 
      {
        parArray[curIdx].E = atof(st.substr(2, st.length() - 1).c_str());
      }

      if(!st.compare(0, 2, "a=")) 
      {
        parArray[curIdx].a = atof(st.substr(2, st.length() - 1).c_str());
      }

      if(!st.compare(0, 2, "n=")) 
      {
        parArray[curIdx].n = atof(st.substr(2, st.length() - 1).c_str());
      }
    }

    if(debug) 
    {
      clog << "Input[" << curIdx + 1 << "]: p0=" << parArray[curIdx].p0 
           << ", Vc=" << parArray[curIdx].Vc
           << ", E=" << parArray[curIdx].E
           << ", a=" << parArray[curIdx].a
           << ", n=" << parArray[curIdx].n
           << endl;
    }

    curIdx++;
  }

  infile.close();

  return lines;
}

void printResults(Results* results, int resultsSize, int paramsCount, char* rowSep, ostream& out) 
{
  // header line
  out << "Time";
  for(int parIdx = 1; parIdx <= paramsCount; parIdx++) out << ",p" << parIdx << ",Delta" << parIdx;
  out << endl;

  // results
  for(int row = 0; row < resultsSize; row++) 
  {
    out << results[row].t;

    for(int col = 0; col < paramsCount; col++) 
    {
      out << "," << setprecision(10) << results[row].p[col] << (rowSep? rowSep : ";") << results[row].delta[col];
    }

    out << endl;
  }
}

void printUsage() 
{
  cout << "usage: model -in input.txt -out output.csv [-tmin min_time] [-tmax max_time] [-tdelta time_delta] [-debug]" << endl;
  cout << "usage: model -help" << endl;
  cout << "K-process modulation with trapeze method utility." << endl;
  cout << "The following options are available:" << endl;
  cout << " -help\tShows this help." << endl;
  cout << " -tmin\tSpecify starting point of time interval. Defaults to 0." << endl;
  cout << " -tmax\tSpecify end of time interval. Defaults to 100." << endl;
  cout << " -tdelta\tSpecify increase of time interval. Defaults to 1." << endl;
  cout << " -in file_name\tSpecifies input file name. File should conatin at least one line in the following format:" << endl;
  cout << "\t\ta=1 n=1 Vc=0.01 p0=-1 E=3" << endl;
  cout << "\t\tEach line will create one result sequence in the output csv." << endl;
  cout << " -out file_name\tSpecifies output file, which will contain time row and rows for each input line. If not specified default output is used." << endl;
  cout << " -sep symbol\tRow separator symbol or string for CSV file. Defaults to ';'." << endl; 
  cout << " -debug\tShows debug output during estimation." << endl;
  cout << "The following options are available:" << endl;
  cout << "Output: F, E as csv file." << endl;
}

int main(int argc, char* argv[])
{
  if (hasOption(argc, argv, "-help")) 
  {
    printUsage();
    return 0;
  }

  debug = hasOption(argc, argv, "-debug") || hasOption(argc, argv, "/?");

  int tmin = 0;
  int tmax = 100;
  int tdelta = 1;

  if(hasOption(argc, argv, "-tmin")) tmin = atoi(getOption(argc, argv, "-tmin"));
  if(hasOption(argc, argv, "-tmax")) tmax = atoi(getOption(argc, argv, "-tmax"));
  if(hasOption(argc, argv, "-tdelta")) tdelta = atoi(getOption(argc, argv, "-tdelta"));

  if(debug) clog << "Time interval: " << tmin << " – " << tmax << ", step: " << tdelta << endl;

  char* paramsFile = getOption(argc, argv, "-in");
  
  if(!paramsFile) {
    cerr << "No input file specified" << endl;
    printUsage();
    return 1;
  }

  InParams* inParams = NULL;
  int paramsCount = readParams(paramsFile, &inParams);

  if(paramsCount == 0) {
    cout << "File is empty. Exiting." << endl;
    return 1;
  }

  int resultsSize = ((tmax - tmin) / tdelta) + 1;

  int idxA = 0;  

  Results* results = new Results[resultsSize];
  for(int i = 0; i < resultsSize; i++) 
  {
    results[i].A = new double[paramsCount];
    results[i].p = new double[paramsCount];
    results[i].delta = new double[paramsCount];

  }

  // Calculations section

  for(int parIdx = 0; parIdx < paramsCount; parIdx++) 
  {
    // ci is holding all input parameters
    InParams& ci = inParams[parIdx];

    // precompute some values
    double currCoef = (ci.Vc / (ci.n + 1)) * exp( (1 - ci.a) / ci.E );  
    if(debug) clog << "Params[" << parIdx + 1 << "] Vc/(n+1) * exp((1-a)/E)=" << currCoef << endl;

    double timeCoef = exp(-1/ci.E);
    if(debug) clog << "Params[" << parIdx + 1 << "] exp(-1/E)=" << timeCoef << endl;

    double deltaCoef = (1 - ci.p0) * ci.Vc * exp(-ci.a/ci.E);
    if(debug) clog << "Params[" << parIdx + 1 << "] 4*(1 - p0) * Vc * exp(-a/E) = " << deltaCoef << endl;

    double eCoef = ci.Vc  * exp((1 - ci.a) / ci.E);
    if(debug) clog << "Params[" << parIdx + 1 << "] Vc * exp((1-a)/E) = " << eCoef << endl;

    int serIdx = 0;

    for (int t = tmin; t <= tmax; t += tdelta) 
    {
      double a = currCoef * (pow(1 + t * timeCoef, ci.n + 1) - 1);
      double p = 1 - (1 - ci.p0) * exp(-a);

      double expPow = eCoef * (pow(1 + t*timeCoef, ci.n + 1) - 1) / (ci.n + 1);

      double delta = deltaCoef * pow(1 + t*timeCoef, ci.n) * exp(-expPow);

      results[serIdx].t = t;
      results[serIdx].A[parIdx] = a;
      results[serIdx].p[parIdx] = p;
      results[serIdx].delta[parIdx] = delta;

      if(debug) clog << "Params[" << parIdx + 1 << setprecision(10)
                     << "] A(" << t << ")=" << a 
                     << "\tP(" << t << ")=" << p 
                     << "\tdelta(" << t << ")=" << delta 
                     << endl;

      serIdx++;
    }

  }
  
  char* rowSep = getOption(argc, argv, "-sep");
  char* outFile = getOption(argc, argv, "-out");

  if(outFile) 
  {
    ofstream ofs(outFile, ofstream::out | ofstream::trunc);
    printResults(results, resultsSize, paramsCount, rowSep, ofs);
  } 
  else 
  {
    printResults(results, resultsSize, paramsCount, rowSep, cout);
  }

  // free all dynamic data

  for(int i = 0; i < resultsSize; i++) {
    delete results[i].A;
    delete results[i].p;
    delete results[i].delta;
  }
  
  delete results;
  delete inParams;

  return 0;
}