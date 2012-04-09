
#include "GenFunc.h"

using namespace std;


// Function to return the index of a sorted vector in ascending order
vector<int> GetIndexForSortedDoubleVector(vector<double>& inputvectortobesorted)
{
    vector<int> index;
    vector<double> vectortobesorted = inputvectortobesorted;
    vector<double> sortedvector = vectortobesorted;
    sort(sortedvector.begin(),sortedvector.end());

    for (unsigned int i = 0; i < sortedvector.size(); ++i)
    {
        int it = find(vectortobesorted.begin(), vectortobesorted.end(), sortedvector.at(i))- vectortobesorted.begin();
        index.push_back(it);
        vectortobesorted.at(it) = 999999;
    }
   return index;
}

vector<int> GetIndexForSortedFloatVector(vector<float>& inputvectortobesorted)
{
    vector<int> index;
    vector<float> vectortobesorted = inputvectortobesorted;
    vector<float> sortedvector = vectortobesorted;
    sort(sortedvector.begin(),sortedvector.end());

    for (unsigned int i = 0; i < sortedvector.size(); ++i)
    {
        int it = find(vectortobesorted.begin(), vectortobesorted.end(), sortedvector.at(i))- vectortobesorted.begin();
        index.push_back(it);
        vectortobesorted.at(it) = 999999;
    }
   return index;
}

// Function to return the index of a sorted vector in Decending order
vector<int> GetIndexForSortedDecendingDoubleVector(vector<double>& inputvectortobesorted)
{
    vector<int> index;
    vector<double> vectortobesorted = inputvectortobesorted;
    vector<double> sortedvector = vectortobesorted;
    sort(sortedvector.begin(),sortedvector.end(), greater<double>());

    for (unsigned int i = 0; i < sortedvector.size(); ++i)
    {
        int it = find(vectortobesorted.begin(), vectortobesorted.end(), sortedvector.at(i))- vectortobesorted.begin();
        index.push_back(it);
        vectortobesorted.at(it) = 999999;
    }
   return index;
}

vector<int> GetIndexForSortedDecendingFloatVector(vector<float>& inputvectortobesorted)
{
    vector<int> index;
    vector<float> vectortobesorted = inputvectortobesorted;
    vector<float> sortedvector = vectortobesorted;
    sort(sortedvector.begin(),sortedvector.end(), greater<float>());

    for (unsigned int i = 0; i < sortedvector.size(); ++i)
    {
        int it = find(vectortobesorted.begin(), vectortobesorted.end(), sortedvector.at(i))- vectortobesorted.begin();
        
        index.push_back(it);
        vectortobesorted.at(it) = 999999;
    }
   return index;
}



int StepFunctionLimitingSelectionNumber(int VectorSize)
{
    int selectnumber = 2;

    if (VectorSize == 1)
        selectnumber = 1;
    else if ( (VectorSize > 1) && (VectorSize <= 6) )
        selectnumber = 2;
    else if ( (VectorSize > 6) && (VectorSize <= 12) )
        selectnumber = 3;
    else if ( (VectorSize > 12) && (VectorSize <= 19) )
        selectnumber = 4;
    else if ( (VectorSize > 19) && (VectorSize <= 25) )
        selectnumber = 5;
    else if ( (VectorSize > 25) && (VectorSize <= 32) )
        selectnumber = 6;
    else if ( (VectorSize > 32) && (VectorSize <= 39) )
        selectnumber = 7;
    else if ( (VectorSize > 39) && (VectorSize <= 47) )
        selectnumber = 8;
    else if ( (VectorSize > 47) && (VectorSize <= 55) )
        selectnumber = 9;
    else if ( (VectorSize > 55) && (VectorSize <= 62) )
        selectnumber = 10;
    else if ( (VectorSize > 62) && (VectorSize <= 70) )
        selectnumber = 11;
    else if ( (VectorSize > 70) && (VectorSize <= 75) )
        selectnumber = 12;
    else if ( (VectorSize > 75) && (VectorSize <= 81) )
        selectnumber = 13;
    else if ( (VectorSize > 81) && (VectorSize <= 86) )
        selectnumber = 14;
    else if (VectorSize > 86)
        selectnumber = 15;

    return selectnumber;
}

// -----------------------------
// String Convertion Functions
// -----------------------------

// Convert string to a unsigned long integer (Michael's code)
unsigned int GetUnsignedInt(const string & str) {
    char * end_ptr = NULL;
    const char * s = str.c_str();
    unsigned int ConvertedStr;

    if (s == end_ptr) {
        cerr << "ERROR: Fail to convert the string " << s << " to an unsigned integer." << endl;
        exit(1);
    } else { // Go ahead and convert the string to unsigned int
        ConvertedStr = strtoul(s, &end_ptr, 10);
        //cout << "CHECKING: What is s in GetUnsignedInt function? " << s << endl;
    }

    return ConvertedStr;
}

// Convert string to a long integer
long int GetInt(const string & str) {
    char * end_ptr = NULL;
    const char * s = str.c_str();
    long int ConvertedStr;

    if (s == end_ptr) {
        cerr << "ERROR: Fail to convert the string " << s << " to a integer." << endl;
        exit(1);
    } else { // Go ahead and convert the string to a long int
        ConvertedStr = strtol(s, &end_ptr, 10);
    }

    return ConvertedStr;
}

// Convert string to a double
double GetDouble(const string & str) {
    char * end_ptr = NULL;
    const char * s = str.c_str();
    double ConvertedStr;

    if (s == end_ptr) {
        cerr << "ERROR: Fail to convert the string " << s << " to a double." << endl;
        exit(1);
    } else { // Go ahead and convert the string to a double
        ConvertedStr = strtod(s, &end_ptr);
    }

    return ConvertedStr;
}

// Function to find the average of the data in Vector type
float FindAverage(vector<unsigned int> & datanum, vector<unsigned int> & dataden) {
    float ave = 0;
    if (dataden.size()>1) {
        for (unsigned int i = 0; i < datanum.size(); i++) {
            ave += (float) datanum.at(i) / (float) dataden.at(i);
        }
    }
    else {
        for (unsigned int i = 0; i < datanum.size(); i++) {
            ave += (float) datanum.at(i) / (float) dataden.at(0);
        }
    }

    ave = ave / (float) datanum.size();

    return ave;
}

// Function to find median of rge data in Vector type
float FindMedian(vector<unsigned int> & datanum, vector<unsigned int> & dataden) {
    vector<float> vec;
    float med;

    // Get the DBSnp fractions for each elements in the data vector
    if (dataden.size()>1) {
        for (unsigned int i = 0; i < datanum.size(); i++) {
            vec.push_back((float) datanum.at(i) / (float) dataden.at(i));
        }
    }
    else {
        for (unsigned int i = 0; i < datanum.size(); i++) {
            vec.push_back((float) datanum.at(i) / (float) dataden.at(0));
        }
    }

    // Sort the elements in tempvec vector (using Selection sort algorithm)
    sort(vec.begin(), vec.end());

    // Find the median of the data
    if (vec.size()%2 == 1) {    // Odd Set
        float index = (vec.size()/2)-0.5;
        med = vec.at(index);
    }
    else {    // even set
        float index = (vec.size()/2);
        med = (vec.at(index)+vec.at(index-1))/2;
    }

    return med;
}

// Function to sort the elements in vector by Selection Sort Method
vector<float> SelectionSortMethod(vector<float> data) {
    for (unsigned int i =0; i< data.size()-1; i++) {
        unsigned int min = i;
        for (unsigned int j=i+1; j < data.size(); j++) {
            if (data.at(j) < data.at(min))
                min = j;
        }

        if (i != min) {
            float swap = data.at(i);
            data.at(i) = data.at(min);
            data.at(min) = swap;
        }
    }

    return data;
}


