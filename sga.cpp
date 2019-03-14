/**************************************************
  SGA - Simple Genetic Algorithm
  ------------------------------
  Sin Bowl minimization 2012
    - recent bug fixes:   09/2012 by John Johnson
    - slight mod to rand: 11/2017 by Scott Gordon
    - modified to work with 3 variable functions and added disaster to
      improve randomization by Andrei Kuzmiankov 11/2017
  ------------------------------
  SGA Adapted from David Goldberg:
    "Genetic Algorithms in Search, Optimization, and Machine Learning"
     Addison-Wesley, 1989.
  ------------------------------
  Unix version
  Compile with:  g++ sga.c
***************************************************/

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <time.h>
using namespace std;

#define POPULATION_SIZE 20    // population size - number of strings
#define CHROM_LENGTH    32    // binary string length of each individual
#define PMUT            0.08  // probability of flipping each bit
#define MAX_GEN         10000 // GA stops after this many generations
#define GEN_REP         1     // report is generated at these intervals
#define ELITE           1     // 1=elitism,  0=no elitism
#define MAXMIN          -1    // -1=minimize, 1=maximize

/***************************************************************
****  random fraction between 0.0 and 1.0                  *****
****************************************************************/
#define fracrand() ((double)rand()/RAND_MAX)

void   initialize_population();
void   crossover(int parent1, int parent2, int child1, int child2);
void   mutation();
void   tselection();
int* decode(int index);
void   getpreviousbest();
double evaluate(int rawX, int rawY);
double* convRange(int rawX, int rawY);
int    flip(double prob);
void   statistics();
void   elite();
void   finalreport();

struct individual
{
  int xV,yV;
  unsigned char string[CHROM_LENGTH];
  double fitness;
};

struct individual pool[POPULATION_SIZE];
struct individual new_pool[POPULATION_SIZE];
struct individual beststring, verybest;

int selected[POPULATION_SIZE];
int generations;

/*********************************************************/
void disaster(){
  int* p;
  for (int i=rand() % 1; i<POPULATION_SIZE-1; )
  {
      for (int j=0; j<CHROM_LENGTH; j++ )
        pool[i].string[j] = flip(0.5);
      p= decode(i);
      pool[i].xV = p[0];
      pool[i].yV = p[1];
      pool[i].fitness = evaluate(pool[i].xV,pool[i].yV );
      i= i+2;
    }
  }


int main()
{
  cout.setf(ios::fixed); cout.setf(ios::showpoint); cout.precision(4);
  int i;
  int *p;
  generations = 0;
  if (MAXMIN==-1) verybest.fitness = 999999; else verybest.fitness=-999999;

  srand(time(NULL));
  initialize_population();
  generations = 1;

  do
  {
    getpreviousbest();

    if((generations % 8) == 0 ){
      disaster();
    }

    /*** SELECTION ***/
    tselection();

    /*** CROSSOVER ***/
    for (i=0; i<POPULATION_SIZE; i=i+2)
      crossover(selected[i],selected[i+1],i,i+1);

    /*** MUTATION ***/
    mutation();

    /*** EVALUATE ***/
    for (i=0; i<POPULATION_SIZE; i++)
    {
      p = decode(i);
      pool[i].xV = p[0];
      pool[i].yV = p[1];
      pool[i].fitness = evaluate(pool[i].xV,pool[i].yV );
    }

    if (ELITE==1)
      elite();

    if (generations % GEN_REP == 0)
      statistics();

  } while (++generations < MAX_GEN);

  finalreport();
  return(0);
}

/*********************************************************
  3-2 Tournament Selection 
**********************************************************/
void tselection()
{ int i;
  for (i=0; i<POPULATION_SIZE; i+=2)
  {
    int r = (int) (fracrand()*POPULATION_SIZE);
    int s = (int) (fracrand()*POPULATION_SIZE);
    int t = (int) (fracrand()*POPULATION_SIZE);

    if ( ((MAXMIN*pool[r].fitness) >= (MAXMIN*pool[s].fitness))
      || ((MAXMIN*pool[r].fitness) >= (MAXMIN*pool[t].fitness)))
    {
      if ((MAXMIN*pool[s].fitness) > (MAXMIN*pool[t].fitness))
        { selected[i] = r; selected[i+1] = s; }
      else
        { selected[i] = r; selected[i+1] = t; }
    }
    else
    {
      if ( ((MAXMIN*pool[s].fitness) >= (MAXMIN*pool[r].fitness))
        || ((MAXMIN*pool[s].fitness) >= (MAXMIN*pool[t].fitness)))
      {
        if ((MAXMIN*pool[r].fitness) > (MAXMIN*pool[t].fitness))
        { selected[i] = s; selected[i+1] = r; }
        else
        { selected[i] = s; selected[i+1] = t; }
      }
      else
      {
        if ( ((MAXMIN*pool[t].fitness) >= (MAXMIN*pool[r].fitness))
          || ((MAXMIN*pool[t].fitness) >= (MAXMIN*pool[s].fitness)))
        {
          if ((MAXMIN*pool[r].fitness) > (MAXMIN*pool[s].fitness))
          { selected[i] = t; selected[i+1] = r; }
          else
          { selected[i] = t; selected[i+1] = s;}
} } } } }

/*********************************************************
  Elitism - copy best string to 0th position of new pool
**********************************************************/
void elite()
{
  int* p =decode(0);
  if ((MAXMIN*beststring.fitness) > (MAXMIN*evaluate(p[0],p[1])))
  {
    pool[0].fitness = beststring.fitness;
    pool[0].xV = beststring.xV;
    pool[0].yV = beststring.yV;
    for (int j=0; j<CHROM_LENGTH; j++)
      pool[0].string[j] = beststring.string[j];
  }
}

/*********************************************************
    Initialize pool to random binary values
**********************************************************/
void initialize_population()
{
  int* p;
  for (int i=0; i<POPULATION_SIZE; i++)
  {
    for (int j=0; j<CHROM_LENGTH; j++ )
      pool[i].string[j] = flip(0.5);
    p= decode(i);
    pool[i].xV = p[0];
    pool[i].yV = p[1];
    pool[i].fitness = evaluate(pool[i].xV,pool[i].yV );
  }
  statistics();
}

/*************************************************************
  - Determine and display best string from previous generation.
  - Maintain very best string from all runs.
**************************************************************/
void getpreviousbest()
{
  cout.setf(ios::fixed); cout.setf(ios::showpoint); cout.precision(4);
  double* p;
  if (MAXMIN==-1) beststring.fitness=999999; else beststring.fitness=-999999;

  for (int i=0; i<POPULATION_SIZE; i++)
  {
    if ((MAXMIN*pool[i].fitness) > (MAXMIN*beststring.fitness))
    {
      beststring.xV = pool[i].xV;
      beststring.xV = pool[i].xV;
      for (int j=0; j<CHROM_LENGTH; j++)
        beststring.string[j] = pool[i].string[j];
      beststring.fitness = pool[i].fitness;
    }
  }   

  if (generations % GEN_REP == 0)
  {
    cout << endl << "   Best string: ";
    for (int j=0;j<CHROM_LENGTH;j++)
      cout << (int) beststring.string[j];
    p = convRange(beststring.xV,beststring.yV);
    cout << " values: " << "X " << p[0]<<"Y "<<p[1]<<endl;
    cout << " fitness: " << beststring.fitness << endl; 
  }

  if ((MAXMIN*beststring.fitness) > (MAXMIN*verybest.fitness))
  {
    verybest.xV = beststring.xV;
      verybest.yV = beststring.yV;
    for (int j=0; j<CHROM_LENGTH; j++)
      verybest.string[j] = beststring.string[j];
    verybest.fitness = beststring.fitness;
  }
}

/*********************************************************
      one-point crossover
**********************************************************/
void crossover (int parent1, int parent2, int child1, int child2)
{
  int i, site;
  site = (int) (fracrand()*CHROM_LENGTH);
  for (i=0; i<CHROM_LENGTH; i++)
  {
    if ((i<=site) || (site==0))
    {
      new_pool[child1].string[i] = pool[parent1].string[i];
      new_pool[child2].string[i] = pool[parent2].string[i];
    }
    else
    {
      new_pool[child1].string[i] = pool[parent2].string[i];
      new_pool[child2].string[i] = pool[parent1].string[i];
    }
  }
}

/*********************************************************
    Bitwise mutation  - also transfers strings to pool
**********************************************************/
void mutation()
{
  int i,j;
  for (i=0; i<POPULATION_SIZE; i++)
  {
    for (j=0; j<CHROM_LENGTH; j++)
      if (flip(PMUT)==1)
        pool[i].string[j] = ~new_pool[i].string[j] & 0x01;
      else
        pool[i].string[j] = new_pool[i].string[j] & 0x01;
  }
}

/*********************************************************
    Convert bitstring to positive integer
**********************************************************/
int* decode(int index)
{
  int v1 = 0, v2= 0;
  int* p = new int[2];
  for (int i=0; i<CHROM_LENGTH/2; i++)
    v1 += (int) pow(2.0,(double)i) * pool[index].string[(CHROM_LENGTH/2)-1-i];
  for (int i=0; i<CHROM_LENGTH/2; i++)
    v2 += (int) pow(2.0,(double)i) * pool[index].string[CHROM_LENGTH-1-i];
 p[0] = v1;
 p[1] =v2;
  return p;
}

/*********************************************************
   F(X) = .1absX - sinX
*********************************************************/
double evaluate(int vx, int vy)
{

  double* p= convRange(vx, vy);
  double g = (double) ((2.0*fabs(3.0*p[1]))-12.0*sin(p[0]))*2.0*p[1]*cos(2.0*p[0]);
  return(g);
}

/*********************************************************
 Convert positive integer to desired floating point range.
 Problem-specific - change for different string lengths
**********************************************************/
double* convRange(int rawX, int rawY)
{
  double* p = new double[2];
  double outval = ((((double)rawX)/65535.0)*120.0)-60.0;
  p[0]=outval;
  double outval1 = ((((double)rawY)/65535.0)*120.0)-60.0;
  p[1] = outval1;
  return p;
}

/*********************************************************
    Do a biased coin toss based on a probability
**********************************************************/
int flip(double prob)
{
  return((fracrand()<prob)?1:0);
}

/*********************************************************
    Report printed at each generation
**********************************************************/
void statistics()
{
  int i,j;
  int* p;
  double* pp;
  cout.setf(ios::fixed); cout.setf(ios::showpoint); cout.precision(4);

  cout << "\nGENERATION: " << generations << endl << "Selected Strings: ";
  for (i=0; i<POPULATION_SIZE; i++) cout << selected[i] << " ";
  cout << endl << "\n\tX\tY\tf(x)\t\tnew_str\t\tX\tY";
  for (i=0; i<POPULATION_SIZE; i++)
  {
    cout << endl << "   ";
    pp = convRange(pool[i].xV,pool[i].yV);
    cout << pp[0] << " "<< pp[1]<<"\t"<<pool[i].fitness << "\t";
    for (j=0; j<CHROM_LENGTH; j++)
      cout << (int) pool[i].string[j];
    p =decode(i);
    pp=convRange(p[0],p[1]);
    cout << "\t" << pp[0]<<" "<<pp[1];
  }
  cout << endl;
}

/*********************************************************
    Report printed at the very end of execution
**********************************************************/
void finalreport()
{
  cout << "=======================================================" << endl;
  cout << "Best result over all generations:" << endl;
  for (int j=0; j<CHROM_LENGTH; j++)
    cout << (int) verybest.string[j];
  cout << endl;
  double* p = convRange(verybest.xV,verybest.yV);
  cout << "Decoded values = " << "X " << p[0]<<"Y "<<p[1]<<endl;
  cout << "  Fitness = " << verybest.fitness << endl;
}

