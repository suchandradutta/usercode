#ifndef __ANAUTIL__HH
#define __ANAUTIL__HH

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <iterator>

using namespace std;

class AnaUtil {

public:
  static void tokenize(const string& str, vector<string>& tokens, const string& delimiters=" ");
  static void bit_print(int value, int pos=32, ostream& os=cout);

  /* PRINT_ELEMENTS()
   * - prints optional C-string optcstr followed by
   * - all elements of the collection coll
   * - separated by spaces
   */
  template <class T>
  static inline void PRINT_ELEMENTS (const T& coll, const char* optcstr="") 
  {
    typename T::const_iterator pos;

    cout << optcstr;
    for (pos = coll.begin(); pos != coll.end(); ++pos)
      cout << *pos << endl;
    cout << endl;
  }

  /* INSERT_ELEMENTS (collection, first, last)
   * - fill values from first to last into the collection
   * - NOTE: NO half-open range
   */
  template <class T>
  static inline void INSERT_ELEMENTS (T& coll, int first, int last)  
  {
      for (int i=first; i<=last; ++i) {
          coll.insert(coll.end(),i);
      }  
  }
  /* INSERT_ELEMENTS (collection, first, last)
   * - fill values from first to last into the collection
   * - NOTE: NO half-open range
   */
  template <class T>
  static inline void COPY_ELEMENTS (const T& sourceColl, T& destColl)  
  {
    destColl.clear();
    typename T::const_iterator pos; 
    for (pos = sourceColl.begin(); pos != sourceColl.end(); ++pos)   
      destColl.push_back(*pos); 
  }
};
#endif
