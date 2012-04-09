// ***************************************************************************
// FastaIndex.h (c) 2010 Erik Garrison <erik.garrison@bc.edu>
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 5 February 2010 (EG)
// ---------------------------------------------------------------------------

#ifndef _FASTA_H
#define _FASTA_H

#include <map>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdint.h>
#include <stdio.h>
#include <algorithm>
#include "LargeFileSupport.h"
#include <sys/stat.h>
#include "split.h"
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>

class FastaIndexEntry
{
    friend std::ostream& operator<<(std::ostream& output, const FastaIndexEntry& e);
    public:
        FastaIndexEntry(std::string name, int length, long long offset, int line_blen, int line_len);
        FastaIndexEntry(void);
        ~FastaIndexEntry(void);
        std::string name;  // sequence name
        int length;  // length of sequence
        long long offset;  // bytes offset of sequence from start of file
        int line_blen;  // line length in bytes, sequence characters
        int line_len;  // line length including newline
        void clear(void);
};

class FastaIndex : public std::map<std::string, FastaIndexEntry>
{
    friend std::ostream& operator<<(std::ostream& output, FastaIndex& i);
    public:
        FastaIndex(void);
        ~FastaIndex(void);
        std::vector<std::string> sequenceNames;
        void indexReference(std::string refName);
        void readIndexFile(std::string fname);
        void writeIndexFile(std::string fname);
        std::ifstream indexFile;
        FastaIndexEntry entry(std::string key);
        void flushEntryToIndex(FastaIndexEntry& entry);
        std::string indexFileExtension(void);
};

class FastaReference
{
    public:
        FastaReference(std::string reffilename);
        std::string filename;
        ~FastaReference(void);
        FILE* file;
        FastaIndex* index;
        std::vector<FastaIndexEntry> findSequencesStartingWith(std::string seqnameStart);
        std::string getSequence(std::string seqname);
        // potentially useful for performance, investigate
        // void getSequence(std::string seqname, std::string& sequence);
        std::string getSubSequence(std::string seqname, int start, int length);
        std::string sequenceNameStartingWith(std::string seqnameStart);
        long unsigned int sequenceLength(std::string seqname);
};

#endif
