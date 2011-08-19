#include <iostream>
#include <climits>

#include "TMath.h"
#include "TLorentzVector.h"
#include "AnaUtil.h"

void AnaUtil::tokenize(const string& str, vector<string>& tokens, const string& delimiters) {
  // Skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);

  // Find first "non-delimiter".
  string::size_type pos = str.find_first_of(delimiters, lastPos);

  while (string::npos != pos || string::npos != lastPos)  {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));

    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);

    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}
void AnaUtil::bit_print(int value, int pos, ostream& os) { 
  static const int INT_BIT = 4*CHAR_BIT; 
  int i, mask = 1 << 31; 
   
  if (pos > INT_BIT) i = INT_BIT; 
  for (i = 1; i <= (INT_BIT - pos); ++i) { 
    value <<= 1; 
  } 
  os.put(' ');
  for (i = 1; i <= pos; ++i) { 
    os.put(((value & mask) == 0) ? '0' : '1'); 
    value <<= 1; 
    if ((INT_BIT - pos + i) % CHAR_BIT ==0 && i != INT_BIT) os.put(' '); 
  } 
  os << endl; 
} 
