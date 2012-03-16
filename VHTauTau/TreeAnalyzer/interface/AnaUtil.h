#ifndef __ANAUTIL__HH
#define __ANAUTIL__HH

#include <climits>
#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <string>
#include <vector>
#include <map>

#include "TMath.h"
#include "TLorentzVector.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::setw;

using std::string;
using std::vector;
using std::map; 
using std::abs;
using std::sqrt;

namespace AnaUtil {

  //public:
  static inline 
  void tokenize(const string& str, vector<string>& tokens, const string& delimiters=" ") {
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
  static inline 
  void bit_print(int value, int pos=32, ostream& os=cout) {
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
  static inline
  double deltaPhi(double phia, double phib) {
    double dphi = abs(phia - phib);
    if (dphi > TMath::Pi()) dphi = 2 * TMath::Pi() - dphi;
    return dphi;
  }
  static inline 
  double deltaPhi(const TLorentzVector &a, const TLorentzVector& b) {
    return deltaPhi(a.Phi(), b.Phi());
  }
  static inline 
  double deltaR(const TLorentzVector &a, const TLorentzVector& b) {
    double dphi = deltaPhi(a,b);
    double deta = a.Eta() - b.Eta();
    return sqrt(dphi * dphi + deta * deta);
  }
  static inline double cutValue(map<string, double>& m, string cname) {
    if (m.find(cname) == m.end()) {
      cerr << ">>> key: " << cname << " not found in the map!" << endl;
      for (map<string,double>::const_iterator jt  = m.begin(); 
                                              jt != m.end(); ++jt)  
        cerr << jt->first << ": " << setw(7) << jt->second << endl;
    }
    assert(m.find(cname) != m.end());
    return m.find(cname)->second;
  }

  // ------------------------------------------------------------------------
  // Convenience routine for filling 1D histograms. We rely on root to keep 
  // track of all the histograms that are booked all over so that we do not 
  // have to use any global variables to save the histogram pointers. Instead, 
  // we use the name of the histograms and gROOT to retrieve them from the 
  // Root object pool whenever necessary. This is the closest one can go to 
  // hbook and ID based histogramming
  // -------------------------------------------------------------------------
  static inline
  TH1* getHist1D(const char* hname) {
    TObject *obj = gDirectory->GetList()->FindObject(hname); 
    if (!obj) {
      cerr << "**** Histogram for <<" << hname << ">> not found!" << endl;
      return 0;
    }
    TH1 *h = 0;
    if (obj->InheritsFrom("TH1D"))
      h = dynamic_cast<TH1D*>(obj);
    else if (obj->InheritsFrom("TH1C"))
      h = dynamic_cast<TH1C*>(obj);
    else if (obj->InheritsFrom("TH1K"))
      h = dynamic_cast<TH1K*>(obj);
    else if (obj->InheritsFrom("TH1S"))
      h = dynamic_cast<TH1S*>(obj);
    else if (obj->InheritsFrom("TH1I"))
      h = dynamic_cast<TH1I*>(obj);
    else
      h = dynamic_cast<TH1F*>(obj);

    if (!h) {
      cerr << "**** 1D Histogram <<" << hname << ">> not found" << endl;
      return 0;
    }
    return h;
  }
  static inline
  TH1* getHist1D(const string& hname) {
    return getHist1D(hname.c_str());
  }

  template <class T>
  static inline bool fillHist1D(const char* hname, T value, double w=1.0) 
  {
    TH1* h = getHist1D(hname);
    if (!h) return false;
    h->Fill(value, w);
    return true;
  }
  template <class T>
  static inline bool fillHist1D(const string& hname, T value, double w=1.0) {
    return fillHist1D(hname.c_str(), value, w);
  }
  // ---------------------------------------------
  // Convenience routine for filling 2D histograms
  // ---------------------------------------------
  template <class T1, class T2>
  static inline bool fillHist2D(const char* hname, T1 xvalue, T2 yvalue, double w=1.0) 
  {
    TObject *obj = gDirectory->GetList()->FindObject(hname); 

    TH2 *h = 0;
    if (obj->InheritsFrom("TH2D"))
      h = dynamic_cast<TH2D*>(obj);
    else if (obj->InheritsFrom("TH2C"))
      h = dynamic_cast<TH2C*>(obj);
    else if (obj->InheritsFrom("TH2S"))
      h = dynamic_cast<TH2S*>(obj);
    else if (obj->InheritsFrom("TH2I"))
      h = dynamic_cast<TH2I*>(obj);
    else
      h = dynamic_cast<TH2F*>(obj);

    if (!h) {
      cerr << "**** 2D Histogram <<" << hname << ">> not found" << endl;
      return false;
    }
    h->Fill(xvalue, yvalue, w);
    return true;
  }
  template <class T1, class T2>
  static inline bool fillHist2D(const string& hname, T1 xvalue, T2 yvalue, double w=1.0) 
  {
    return fillHist2D(hname.c_str(), xvalue, yvalue, w);
  }
  // --------------------------------------------------
  // Convenience routine for filling profile histograms
  // --------------------------------------------------
  static inline
  TProfile* getProfile(const char* hname) 
  {
    TProfile *h = dynamic_cast<TProfile*>(gDirectory->GetList()->FindObject(hname));
    if (!h) {
      cerr << "**** Profile Histogram <<" << hname << ">> not found" << endl;
      return 0;
    }
    return h;
  }
  static inline
  TProfile* getProfile(const string& hname) 
  {
    return getProfile(hname.c_str());
  }

  static inline 
  bool fillProfile(const char *hname, float xvalue, float yvalue, double w=1.0) 
  {
    TProfile *h = getProfile(hname);
    if (!h) {
      // cerr << "**** Profile Histogram <<" << hname << ">> not found" << endl;
      return kFALSE;
    }

    h->Fill(xvalue, yvalue, w);
    return true;
  }
  static inline 
  bool fillProfile(const string& hname, float xvalue, float yvalue, double w=1.0) 
  {
    return fillProfile(hname.c_str(), xvalue, yvalue, w);
  }
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
