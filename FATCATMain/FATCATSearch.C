// FATCATSearch.pl in C++
// In order to Use OpenMP to do database searching in parallel
// Perl does not support OpenMP
// Running FATCAT for a query against a list of structures
// Compile: g++ -fopenmp -o FATCATSearch FATCATSearch.C
// Example usage: ./FATCATSearch 1dzrA.pdb Jan14-15_pdb_90.descript -i1 ./ -i2 /data/www/fatcat-data/pdb_pdb -o [/path/]mult.aln
// Zhanwen Li @ Mar. 2015


#include <cstdlib> //getenv
#include "omp.h"
#include "basic.h"

using namespace ARRAY;

// global vars
char* fatcatpath;
char prog[200];


int main (int argc, char* argv[]) {

  time_t timebegin = time(NULL);

  // set up env vars and global vars
  fatcatpath = getenv("FATCAT");
  if (fatcatpath==NULL) { printf("Stop: FATCAT environment variable is not found\n"); exit(1); }
  sprintf(prog, "%s/FATCATMain/FATCAT", fatcatpath);

  // check and validate inputs
  if (argc<5) { 
    printf("FATCATSearch  <query> <target-list> <-o output-file> FATCAT-parameters(refer FATCAT usage: -o outfile must be given)\n");
    printf(" when only part of the items from the target-list are used (e.g divide job in cluster), use target-list:beg:end\n");
    printf(" for instance, /path/to/targetlist/scop175_90.descript, i.e search query agains targets in file scop175_90.description\n\n");
    exit(1);
  }
  char query[200]; strcpy(query, argv[1]);
  char targetListF[200]; strcpy(targetListF, argv[2]);
  char fatcatOutFile[2000]; fatcatOutFile[0]='\0';
  int fatcatI[100]; 
  int j=-1;
  for (int i=3; i<argc; i++) {
    if (!strcmp(argv[i], "-o") && argc>i+1) { strcpy(fatcatOutFile, argv[++i]); }
    else if (!strcmp(argv[i], "-q") || !strcmp(argv[i], "-m")) { continue; } // skip -q, -m
    else { fatcatI[++j]=i;};
  }
  if (fatcatOutFile[0]=='\0') { printf("Please use \"-o outfile\" for FATCAT program output\n"); exit(1);}

  // output file and output file dir
  char outDir[500] = {'.', '\0'}; // default is current dir
  char outFile[500]; strcpy(outFile, fatcatOutFile);
  int p=0; 
  for (int i=0; i<strlen(fatcatOutFile); i++) {
    if (fatcatOutFile[i] == '/') { p=i; }
  }
  if (p != 0) {   strcpy(outFile, fatcatOutFile+p+1);  strncpy(outDir, fatcatOutFile, p); } 

  //printf("Dir: %s; File: %s\n", outDir, outFile);

  char fatcatPars[200]; fatcatPars[0]='\0';
  for (int k=0; k<=j; k++) { strcat(fatcatPars, argv[fatcatI[k]]); strcat(fatcatPars, " ");}
  //printf("FATCAT pars: %s\n", fatcatPars);
  //printf("FATCAT output file: %s\n", fatcatOutFile);

  // read target list into an array of strings
  // max of 100M targets, each line with max len of 100 chars for PDB codes
  ifstream in(targetListF);
  if (!in) { printf("Open file %s error\n", targetListF); exit(1); }
  int max=1000000;
  char **targetList = NewMatrix <char> (max, 100); 
  int n=0;
  char str[5000];
  while(!in.eof()) {
    str[0]='\0'; in.getline(str,5000); 
    if (str[0]=='#' || str[0]=='\n' || str[0]=='\0') continue;
    sscanf(str, "%s", targetList[n]);
    //printf("%s\n", str);
    n++;
  }
  in.close();

  //printf("Total number of target: %d\n", n);

  // check whether file reading is correct
  //for (int i=0; i<n; i++) {    printf("%s\n", targetList[i]);  }
  
  // OpenMp: parallelly compute pairwise alignments 
  int threads=0;
  char rmComd[200]; sprintf(rmComd, "rm -f %s/_%s.????_", outDir, outFile); system(rmComd);
  #pragma omp parallel for
  for (int i=0; i<n; i++) {
    threads = omp_get_num_threads(); // OpenMP can auto-detect by the number of cores
    int id = omp_get_thread_num();
    char command[500];
    sprintf(command, "%s -p1 %s -p2 %s.pdb  %s -q >> %s/_%s.%04d_", prog, query, targetList[i], fatcatPars, outDir, outFile, id);
    //printf("%s\n", command);
    system(command);
  }

  //printf("Number of threads: %d\n", threads);
  // free memory
  if (targetList!=NULL) { DelMatrix <char> (targetList, max); }
  
  // Combine the results and delete the result files from each threads
  char comd[200]; sprintf(comd, "cat %s/_%s.????_ > %s", outDir, outFile, fatcatOutFile);
  system(comd);
  sprintf(comd, "rm -f %s/_%s.????_", outDir, outFile); system(comd);
  
  time_t timeend = time(NULL);
  double diff = difftime(timeend, timebegin);
  printf("Total time %.1f\n", diff);

}

