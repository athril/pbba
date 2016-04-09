#ifndef __TOKENIZE_H__
#define __TOKENIZE_H__

#include <iostream>
#include <sstream>
#include <boost/bimap.hpp>

using namespace std;
using namespace boost;

typedef boost::bimap< std::string, int > relation;
typedef relation::value_type rel;
const size_t SIZEX=7;

struct Data {
  string content[SIZEX];
  Data(tokenizer<escaped_list_separator<char> > &tok);

  friend std::ostream &operator<<(std::ostream &os, Data const &d) {
    for (int i=0; i<SIZEX; i++)
      os << d.content[i];
    return os;
  }

  friend std::istream &operator>>(std::istream &is, Data &d) {
    for (int i=0; i<SIZEX; i++)
      is >> d.content[i];
    return is;
  }
};


void preverify(char *filename, long &noRows, long &noCols);
int* tokenize(char *filename, list<Data> &dataset, int *input, vector<relation> &pos);

#endif


