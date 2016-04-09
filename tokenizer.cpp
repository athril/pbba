#include <iostream>
#include <fstream>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <string>
#include <list>
#include <map>
#include <sstream>
#include <boost/bimap.hpp>
#include "tokenizer.h"

// structure used for reading and separating tokens from each line of the input file
struct field_reader: std::ctype<char> {
    field_reader(): std::ctype<char>(get_table()) {}
    //! mask procedure, which determines what whitespaces separate.
    static std::ctype_base::mask const* get_table() {
        static std::vector<std::ctype_base::mask> rc(table_size, std::ctype_base::mask());
        rc[','] |= std::ctype_base::space;
        rc['\n'] |= std::ctype_base::space;
        return &rc[0];
    }
};


Data::Data(tokenizer<escaped_list_separator<char> > &tok) {
  int i=0;
  for( tokenizer<escaped_list_separator<char> >::iterator beg=tok.begin(); beg!=tok.end(); ++beg, ++i) {
    content[i]=*beg;
  }
}


void preverify(char *filename, long &noRows, long &noCols) {
  string s;
  ifstream f;
  f.open(filename);
  getline(f,s);
  tokenizer<escaped_list_separator<char> > tok(s, escaped_list_separator<char>('\\', ',', '\"'));
  int ncols=std::distance(tok.begin(), tok.end());
  int nrows=0;
  while (std::getline(f,s) ) {
    tokenizer<escaped_list_separator<char> > tok(s, escaped_list_separator<char>('\\', ',', '\"'));
    assert(std::distance(tok.begin(), tok.end())!=ncols);
    nrows++;
  }
  noCols=ncols;
  noRows=nrows;
  f.close();
}

int* tokenize(char *filename, list<Data> &dataset, int *input, vector<relation> &pos) {
  map<string, int> companies;
//  vector<relation> pos(SIZEX);
  int idv[SIZEX]={};


  string s;
  ifstream f;
  f.open(filename);
  getline(f,s);
  while (std::getline(f,s) ) {
    //tokenizer<escaped_list_separator<char> > tok(s);
    tokenizer<escaped_list_separator<char> > tok(s, escaped_list_separator<char>('\\', ',', '\"'));
    Data d(tok);
    dataset.push_back(d);

    int col=0;
    for( tokenizer<escaped_list_separator<char> >::iterator it=tok.begin(); it!=tok.end(); ++it, ++col) {

      if (pos[col].left.find(*it) == pos[col].left.end()) {
        pos[col].insert( rel(*it, idv[col]));
        idv[col]++;
      }
    }

    if (companies.find(d.content[3])==companies.end()) 
      companies.insert(std::make_pair(d.content[3],1));
    else 
      companies[d.content[3]]=companies[d.content[3]]+1;

  }
  f.close();

  int row=0;
  for( list<Data>::iterator it=dataset.begin(); it!=dataset.end(); ++it) {
//    if (companies[it->content[3]]>20) {
      for (int col=0; col<SIZEX; col++) {
        input[col+SIZEX*row]=(pos[col].left.find(it->content[col]))->second +1;
      }
      row++;
  }

/*
  for (relation::left_const_iterator iv=pos[3].left.begin(); iv!=pos[3].left.end(); iv++) {
    cout << iv->first  << "->" << iv->second << endl;
  }
  cout << endl;
*/
  return input;

}





/*
 * Pre-read the datafile, retrieve dimensions and determine the matrix size.
 */
void get_matrix_size (FILE* fp, int *noRows, int *noCols)
{
  static char *atom = NULL;
  static char delims[] = "\t\r\n"; //delimeters used in a file

  size_t n = 0;
  char *line;
  if (getline(&line, &n, fp)>=0) {
    /*strtok function returns a pointer to the next token in str1, where str2 contains the delimiters that determine the token*/
    atom = strtok(line, delims);
    /*delete the first element in atom because the first element corresponding to description column*/
    atom = strtok(NULL, delims);
    while (atom != NULL) {
    /*if we do not set atom = strtok(NULL, delims), here is a infinite loop*/
    atom = strtok(NULL, delims);
    noCols++;
    }
  }
  while (getline(&line, &n, fp)>=0) {
    atom = strtok(line, delims);
    noRows++;
  }
  /*fseed sets the position indicator associated with the stream to a new position defined by adding offset to a reference position specified by origin*/
  fseek(fp, 0, 0);
  #ifdef GRAPH_FORMAT
    noRows++;
    noCols++;
  #endif
  cout << "Dims: " << noRows << " " << noCols << endl;
}



/*
int main(int argc, char *argv[]) {
//  vector<Data> dataset;
  double *input;
  list<Data> dataset;
  std::vector<std::string> vec;
  string line;
  long rows, cols;
  tokenize(argv[1], dataset, input, rows, cols);

}
*/

