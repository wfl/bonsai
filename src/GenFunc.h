
#ifndef _GENFUNC_H
#define	_GENFUNC_H

#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <algorithm>


std::vector<int> GetIndexForSortedDoubleVector(std::vector<double>& inputvectortobesorted);
std::vector<int> GetIndexForSortedDecendingDoubleVector(std::vector<double>& inputvectortobesorted);

std::vector<int> GetIndexForSortedFloatVector(std::vector<float>& inputvectortobesorted);
std::vector<int> GetIndexForSortedDecendingFloatVector(std::vector<float>& inputvectortobesorted);

int StepFunctionLimitingSelectionNumber(int VectorSize);

//Strings conversion function
unsigned int GetUnsignedInt(const std::string & str);
long int GetInt(const std::string & str);
double GetDouble(const std::string & str);

// Function to find the average of the data in Vector type
float FindAverage(std::vector<unsigned int> & datanum, std::vector<unsigned int> & dataden);
// Function to find the Median of the data in Vector type
float FindMedian(std::vector<unsigned int> & datanum, std::vector<unsigned int> & dataden);
// Function to sort the data (vector) using the Selection Sort Algorithm
std::vector<float> SelectionSortMethod(std::vector<float> & data);

#endif	/* _GENFUNC_H */

