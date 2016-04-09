#include <iostream>
#include <stdlib.h>
#include <climits>
#include <float.h>
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>
#include <boost/foreach.hpp>
#include <vector>
#include <iterator>
#include <cstdio>
#include <map>
#include <queue>
#include <sys/time.h>
#include <cstring>
#include <list>
#include <boost/lexical_cast.hpp>
#include <sstream>
#include <string>
#include <boost/date_time.hpp>
#include <boost/date_time/date_facet.hpp>
#include <boost/date_time/time_facet.hpp>
#include "tokenizer.h"

#define GRAPH_FORMAT 1
#define MAX(a,b)  ((a)>(b)?(a):(b))
#define MIN(a,b)  ((a)<(b)?(a):(b))
#define ABS(x)    ((x)>0?(x):-(x))

//#define REAL 1
//#define DEBUG 1
//#define DEBUG2 1
//#define WEEKDAYS 0

using namespace std;
namespace bt = boost::posix_time;

typedef boost::unordered_map< boost::dynamic_bitset<>, boost::dynamic_bitset<> > biclset; //used for storing a bicluster
typedef boost::unordered_map< boost::dynamic_bitset<>, int > subset; //used for storing rows/columns of a bicluster
typedef boost::unordered_map< boost::dynamic_bitset<>, int >::iterator subsit; //iterator on rows/columns of a bicluster


subset *subsets; //list of all subset from a given pattern

long   noRows; //noRows of input file
long   noCols; //noCols of input

size_t steps; //PBBA input parameter: level of discretization
unsigned long   minNoRows;  //PBBA input parameter: minimal number of rows of bilusters
unsigned long   minNoCols;  //PBBA input parameter: minimal number of rows of bilusters

//static char *atom = NULL;
//static char delims[] = "\t\r\n"; //delimeters used in a file
int id=0; //id of bicluster


namespace boost {
    template <typename B, typename A>
     std::size_t hash_value(const boost::dynamic_bitset<B, A>& bs) {
        std::vector<B, A> v;
        boost::to_block_range(bs, std::back_inserter(v));
        return boost::hash_value(v);
    }
}


struct Prop {
  int level;
  int colID;
};


bool operator<(const Prop& p1, const Prop& p2) {
    return p1.level < p2.level;
}


ostream& operator<<(ostream& out, const boost::dynamic_bitset<> &b)  {
  for (uint j=0; j<b.size(); j++) {
    if (b[j]>0)
      out<< j << " ";
  }
  return out;
}

template <typename T>
ostream& printData(ostream& out, const boost::dynamic_bitset<> &b, vector<T> v) {
  for (uint j=0; j<b.size(); j++) {
    if (b[j]>0)
      out<< v.at(j) << " ";
  }
  return out;
}



ostream& operator<<(ostream& out, const biclset &map)  {

  BOOST_FOREACH(const biclset::value_type &b, map)  {
//    if( b.first.count()<minNoCols || b.second.count()<minNoRows && b.first.size()>0 && b.second.size()>0 )
    if( b.first.count()<minNoCols || b.second.count()<minNoRows )
      continue;
    out << "#" << id << ":"<< b.second <<"| " <<b.first;
    out<< endl;
    id++;
    #ifndef GRAPH_FORMAT
    if (id>1000) {
	cout << "Showing TOP 1000 of " << map.size() << " elements! Halting" << endl;
	exit(0);
    }
    #endif
  }
  return out;
}



void showArray(int *input, int row) {
    for (int j=0; j<noCols;++j) {
      if (input[j+row*noCols]) 
      cout << row << " ";
    else
      cout << input[j+row*noCols]<< " ";
    }
    cout <<endl;
}


/*
 * constructs the input matrix for PBBA algorithm
 */
int* construct(int *input, size_t nrows, size_t ncols) {
  int* array=new int[nrows*ncols];
  subsets = new boost::unordered_map< boost::dynamic_bitset<>, int >[nrows];
//  subset = new boost::unordered_map< int, boost::dynamic_bitset<> >)[nrows];

  for (uint j=0; j<ncols;++j) {
    for (uint a=0; a<nrows; a++) {
      array[j+a*ncols]=0;
     }
  }


  for (uint j=0; j<ncols;++j) {
    for (uint a=0; a<nrows; a++) {
      for (uint b=a+1; b<nrows; b++) {

        if (input[j+a*ncols]!=0 && input[j+b*ncols]!=0 ) {
          if ( input[j+a*ncols] == input[j+b*ncols] ) {
            if (array[j+a*ncols]!=0) {
              array[j+a*ncols]=MIN(array[j+a*ncols],(int)a+1);
            }
            else  {
              array[j+a*ncols]=a+1;
	    }
            array[j+b*ncols]=MAX(array[j+b*ncols],(int)a+1);
            break;
          }
        }  //new trans
      } // for b loop
    } // for a loop
  }

  for (uint j=0; j<ncols;++j) 
    for (uint a=0; a<nrows; a++) 
      if (input[j+a*ncols]!=0 && array[j+a*ncols]==0) {
        array[j+a*ncols]=a+1;
      }
  return array;
}


/*
 * restricts the specific column patterns in rows, so that they weren't included
 */
bool restricted (subset *s, boost::dynamic_bitset<> *common, int value) {
 subsit f = s->find(*common);
 if (f!=s->end()) {
   #ifdef DEBUG3
   cout << "<<<" << (*f).first <<endl;
   #endif
//   if ( (*f).second==value ) 
     return true;
 }
  for (subsit it = s->begin();  it != s->end();  ++it) {
    if ( (*common).is_subset_of( (*it).first) ) {
//      if ((*it).second==value) {
         #ifdef DEBUG3
         cout << "<<<" << (*it).first <<endl;
         #endif
         return true;
//      }
    }
  }
  return false;
}



/*
 * performs motif propagation for seed row. Seed row.
 */
void insert (biclset &map, int seed, std::pair< boost::dynamic_bitset<>, boost::dynamic_bitset<> > insert) {
  #ifdef DEBUG2
  cout << insert.first<<"|" << insert.second << " "<<endl;
  cout << "SIMILARITIES FOR SEED " << seed << endl;
  #endif
  queue< std::pair< boost::dynamic_bitset<>, boost::dynamic_bitset<> > > changes;
  boost::dynamic_bitset<> *dynp = &insert.first;


  if (restricted(&subsets[seed], dynp,1) || insert.first.none()) {
        #ifdef DEBUG2
      cout <<"COVERED!" << endl;
        #endif
      return;
  }

  if (insert.first.any()) {
    if (map.find(insert.first)==map.end()) {
      #ifdef DEBUG2
      cout << "S" << insert.first << "|" << insert.second << endl;
      #endif
      map[insert.first]=insert.second;
    } else {
      #ifdef DEBUG2
      cout << "R" << insert.first << "|" << map[insert.first] << endl;
      #endif
      map[insert.first]=map[insert.first]|insert.second;
    }
  }


  boost::dynamic_bitset<> common;
  BOOST_FOREACH(biclset::value_type &b, map) {
  // j-> b
  // rows -> b.second

  common = b.first & insert.first;    //coland
  dynp=&common;

  if (restricted(&subsets[seed], dynp,1)) {
      #ifdef DEBUG2
      cout <<"COV-subs!" << endl;
      #endif
      continue;
    }
    if ( common.any() && common.size()>=minNoCols ) {
//    if ( common.count() >= minNoCols ) {
      #ifdef DEBUG2
      cout << "=" <<common <<" "<<endl;
      #endif
      biclset::iterator f = map.find(common);
      if (f!=map.end()) {
        changes.push(std::make_pair(common, b.second | insert.second | f->second));
	#ifdef DEBUG2
        cout << "M" << changes.back().second << endl;
        #endif
      } else {
          changes.push(std::make_pair(common, insert.second|b.second));
	#ifdef DEBUG2
        cout << "+" << changes.back().second << endl;
	#endif
      }
    }
    common.reset();
  }  //end boost-foreach
  
  while (!changes.empty()) {
      biclset::iterator f = map.find(changes.front().first);
      if (f!=map.end())
        map[changes.front().first]|=changes.front().second;
      else
        map[changes.front().first]=changes.front().second;
    changes.pop();

  }
  #ifdef DEBUG2
  cout <<endl;
  #endif

}

/*
 * Performs propagation of the given motif (seed) 
 */

void nib(int *input, int *ranked, int seed, size_t nrows, size_t ncols) {
  biclset biclseed;
//  cout << "NIB:" <<seed << " " << nrows << " " << ncols << endl;
  boost::dynamic_bitset<> pattern(ncols);
  boost::dynamic_bitset<> row(nrows);
  priority_queue<Prop> mask;

//  int* mask=new int[ncols];

  for (uint j=0; j<ncols; j++)  {
    if (ranked[j+seed*ncols]!=0) {
      pattern.set(j);
      //Prop elem; //={seed+1,j};
      Prop elem={seed,j};
//      elem.level=seed+1;
      elem.level=ranked[j+seed*ncols];
      elem.colID=j;
      mask.push(elem);
    }
  }

  row.set(seed);
  #ifdef DEBUG2
    cout << seed << " " << pattern <<"|" << row << endl;
    cout << "BICLSEED" << endl << biclseed <<endl;
  #endif

  if (pattern.count()<minNoCols) {
    pattern.clear();
    row.clear();
    while (!mask.empty())
      mask.pop();
    return;
  }
  if (pattern.count()>=minNoCols)
    insert(biclseed, seed, std::make_pair(pattern,row));

  int nextlevel=0;
  for (int level=seed+1; level>=0; level=nextlevel,nextlevel=0) {
    boost::dynamic_bitset<> resembles(ncols);
    boost::dynamic_bitset<> rowsim(nrows);
    #ifdef DEBUG
      cout << "level=" << level << endl;
    #endif
    //find next transition
    if (mask.empty())
      break;
    nextlevel=mask.top().level;
    //    if (nextlevel==0)
    //      break;    
    #ifdef DEBUG
      cout << "nextlevel=" << nextlevel << " " << mask.top().level << endl;
    #endif

    while ( mask.top().level==nextlevel && !mask.empty()) {
      int j=mask.top().colID;
      resembles.set(j);

      if ( ranked[j+(mask.top().level-1)*ncols] != mask.top().level ) {
        Prop elem;
        elem.level=ranked[j+(mask.top().level-1)*ncols];
        elem.colID=j;
        mask.push(elem);
      }
      mask.pop();
    }

    if (resembles.count()>=minNoCols) {
      subsets[nextlevel-1].insert(boost::unordered_map< boost::dynamic_bitset<>, int >::value_type(resembles,1));
      rowsim.set(nextlevel-1);
      rowsim.set(seed);
      insert( biclseed, seed, std::make_pair(resembles,rowsim));
    }
//    resembles.reset();
//    rowsim.reset();
    if (biclseed.size()>500000) {
      cerr << "Halted prematurely >500000 biclusters" << endl;
      exit(1);
    }
  }  //end for
  cout << biclseed;
  row.clear();
  pattern.clear();
  biclseed.clear();
}








void nib(int *input, int *ranked, size_t nrows, size_t ncols) {
  for(int i=nrows-1; i>=0; i--) {
    cout << "** SEED "<<i << " **"<<endl;
    nib(input,ranked,i,nrows,ncols);
  }
}


void preprocess_array(double *raw, int steps, int *input) {

  list<double> values;
  srand(time(NULL));
  double *stepped = new double[steps];

  for (int i = 0L; i < noRows; i++) {
    for (int j = 0L; j < noCols; j++) {
      double d=raw[j+i*noCols];
      values.push_back(d);
    }
  }

  int level= (noCols*noRows %steps ?  noCols*noRows/steps+1 :  noCols*noRows/steps );
  values.sort();

  int i=0;
  for(list<double>::iterator it=values.begin(); it!=values.end();it++,i++) {
      stepped[i/level]=(*it);

  }

  //cout << "W:"<<*(values.begin()) << "-"<< stepped[0] << "|" << *values.rbegin() << "-" << stepped[steps-2] << "|" << stepped[steps-1] << endl;
  for (int i = 0L; i < noRows; i++) {
    for (int j = 0L; j < noCols; j++) {
      int x=0;
      while ( double v=stepped[x] < raw[j+i*noCols]) {
        x++;
      }
      input[j+i*noCols]= x<steps ? x+1 : steps;
    }
  }
  for (int i=0; i<noRows; i++) {
    for (int j=0; j<noCols; j++) {
      cout << input[j+i*noCols] << "x";
    }
    cout << endl;
  }
  values.clear();
  for (int i=0; i<steps; i++)
    stepped[i]=0;
  delete stepped;
}

void preprocess_graph(double *raw, int *input) {
    for (int i=0; i<noRows; i++) {
      for (int j=0; j<noCols; j++) {
        input[j+i*noCols] = raw[j+i*noCols];
        cout << input[j+i*noCols] << "\t";
      }
      cout << endl;
    }
}



int main( int argc, char *argv [])
{
  if(argc!=5) {
    perror("Usage: ./pbba <file.csv> steps minNoRows minNoCols");
    exit(1);
  }
  srand(time(NULL));
  steps=atoi(argv[2]);
  minNoRows=atoi(argv[3]);
  minNoCols=atoi(argv[4]);
  list<Data> dataset;
  std::vector<std::string> vec;
  string line;

  preverify(argv[1],noRows,noCols);
  int *input = new int[noRows*noCols];
  vector<relation> pos(noCols);
  tokenize(argv[1], dataset, input, pos);
/*
  //fclose(fp);
  #ifdef GRAPH_FORMAT
  preprocess_graph(raw, input);
  #else
  preprocess_array(raw, steps, input);
  #endif
  */
  int *out= new int[noRows*noCols];
  out=construct(input,noRows,noCols);

/*
  for (int i = 0L; i < noRows; i++) {
    for (int j = 0L; j < noCols; j++) {
    cout << input[j+noCols*i] << " ";
    }
  }
*/
  struct timeval begin_t, end_t;
  gettimeofday(&begin_t, NULL);
  nib(input, out, noRows, noCols);
  gettimeofday(&end_t, NULL);
  cout << endl;
  cout <<  "begin.tv_sec = " << begin_t.tv_sec << ", beg.tv_usec = " <<  begin_t.tv_usec << endl;
  cout <<  "end.tv_sec = " << end_t.tv_sec << ", end.tv_usec = " <<  end_t.tv_usec << endl;
  long elapsed_utime;    // processing time in microseconds
  long elapsed_seconds;  // diff between seconds counter
  long elapsed_useconds; // diff between microseconds counter
  elapsed_seconds  = end_t.tv_sec  - begin_t.tv_sec;
  elapsed_useconds = end_t.tv_usec - begin_t.tv_usec;
  elapsed_utime = (elapsed_seconds) * 1000000 + elapsed_useconds;
  cout << endl;


//  cout << "Processing time = " << elapsed_utime << " microsecond(s)" << endl;
  long microseconds = elapsed_utime % 1000000;
  long seconds = (elapsed_utime / 1000000) % 60;
  long minutes = (elapsed_utime / 1000000) / 60;
  cout  << "Processing time = " << minutes << " minute(s), " <<  seconds << " second(s), " << microseconds << " microsecond(s)" << endl;
  cout << elapsed_utime << endl;
  return 0;
} /* main */

