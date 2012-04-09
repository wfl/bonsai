
#include <string>
#include <iostream>
#include <sstream>

using namespace std;

class IndelAllele {

    friend ostream& operator<<(ostream&, const IndelAllele&);
    friend bool operator==(const IndelAllele&, const IndelAllele&);
    friend bool operator!=(const IndelAllele&, const IndelAllele&);
    friend bool operator<(const IndelAllele&, const IndelAllele&);
  
    
public:
    bool insertion;
    bool deletion;
    int length;
    int position;
    string sequence;

    IndelAllele(bool i, bool d, int l, int p, string s)
        : insertion(i), deletion(d), length(l), position(p), sequence(s)
    { }

};

ostream& operator<<(ostream& out, const IndelAllele& indel) {

    string t;
    if (indel.insertion)
        t = "i";
    else if (indel.deletion)
        t= "d";

    out << t <<  ":" << indel.position << ":" << indel.sequence;
    return out;
    }

bool operator==(const IndelAllele& a, const IndelAllele& b) {
    return (a.insertion == b.insertion
            && a.deletion == b.deletion
            && a.length == b.length
            && a.position == b.position
            && a.sequence == b.sequence);
    }

bool operator!=(const IndelAllele& a, const IndelAllele& b) {
    return !(a==b);
    }

bool operator<(const IndelAllele& a, const IndelAllele& b) {
    ostringstream as, bs;
    as << a;
    bs << b;
    return as.str() < bs.str();
    }

