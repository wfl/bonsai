
#include <sstream>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include <utility>

#include "Realign.h"
#include "ModifiedSmithWatermanGotoh.h"
#include "GenFunc.h"
#include "Fasta.h"
#include "api/BamAux.h"

using namespace std;
using namespace BamTools;


// Public Functions

RealignFunctionsClass::RealignFunctionsClass()
{
    m_isLowComplexityRegion = false;
}

RealignFunctionsClass::~RealignFunctionsClass()
{ }

void RealignFunctionsClass::Clear()
{
    m_Haplotypevec.clear();
    m_ChosenBestHaplotypesIndex.clear();

    //m_CoverageByPositionInWindow.clear();
    m_Position_IndelCounts.clear();
    m_Position_MismatchCounts.clear();

    m_MismatchesLociList.clear();
    m_DeletionLociList.clear();
    m_InsertionLociList.clear();

    m_CoverageByPositionInWindow.clear();
    m_isLowComplexityRegion = false;
}

// ---------------------------
// Group Alignments function
// ---------------------------
float RealignFunctionsClass::getBaseQualityASCII2float(string BaseQualstr)
{
    char *writable = new char[BaseQualstr.size() + 1];
    copy(BaseQualstr.begin(), BaseQualstr.end(), writable);
    writable[BaseQualstr.size()] = '\0';
    float f_BaseQual = (static_cast<int> (*writable) - 33);
    delete[] writable;

    return f_BaseQual;
}

float RealignFunctionsClass::getBaseQualityASCII2float(char BaseQualstr)
{
    float f_BaseQual = (static_cast<int> (BaseQualstr) - 33);

    return f_BaseQual;
}


void RealignFunctionsClass::UpdateStartandEndPositionForHaplotype(BamAlignment& alignment, vector<int>& minmaxRefSeqPos)
{
    int flanking_leftwindow = alignment.Length;

    // Check for Alignment positions
    if (minmaxRefSeqPos.at(0) < 0)
    {
        if ( (alignment.Position-flanking_leftwindow) > 0)
            minmaxRefSeqPos.at(0) = alignment.Position - flanking_leftwindow;
        else
            minmaxRefSeqPos.at(0) = 0;

        minmaxRefSeqPos.at(1) = alignment.GetEndPosition();
    }
    else
    {
        if (minmaxRefSeqPos.at(0) > alignment.Position-flanking_leftwindow)
        {
            if ( (alignment.Position-flanking_leftwindow) > 0)
                minmaxRefSeqPos.at(0) = alignment.Position - flanking_leftwindow;
            else
                minmaxRefSeqPos.at(0) = 0;
        }
        if (minmaxRefSeqPos.at(1) < alignment.GetEndPosition())
            minmaxRefSeqPos.at(1) = alignment.GetEndPosition();
    }

    //cerr << alignment.Name << " :  Minmax " << alignment.Position << " End " << alignment.GetEndPosition() << " NewMIN " <<  minmaxRefSeqPos.at(0) << endl;
}


void RealignFunctionsClass::ParseAlignmentsAndExtractVariantsByBaseQualv(vector<SalRealignInfo>& AlGroups, SalRealignInfo& alr, vector<SalRealignInfo>::iterator& Iter,
        BamAlignment& alignment, const string& referenceSequence, vector<int>& minmaxRefSeqPos, int winstartpos, int winendpos, float MinBaseQual, bool IsupdateMap)
{
    //cerr << "#Start: Start Parse " << endl;

     // Update the minimum and maxsimum of start and end position of the alignments
    UpdateStartandEndPositionForHaplotype(alr.al, minmaxRefSeqPos);

    string str = "ZZZZZZZZZZZZZZZZZZZZZ";
    //string str = "tL2";

    int rp = 0; // read position, 0-based relative to read
    int cgp = 0; // current sequence position, 0-based relative to currentSequence
    int sp = alr.al.Position; // sequence position

    vector<CigarOp>::const_iterator cigarIter = alr.al.CigarData.begin();
    vector<CigarOp>::const_iterator cigarEnd  = alr.al.CigarData.end();

    // Parameters required for update
    float numMismatches = 0;
    float numNonIndelNucleotide = 0;
    int longestINDELlen = 0;

    int currentNumOfIDS = 0;
    int currentNumOfm = 0;
    int longestContiguousMismatch = 0;
    int longestSoftclippingLength = 0;
    //
    if (alr.al.IsDuplicate() || !alr.al.IsPrimaryAlignment())
    {
        for ( ; cigarIter != cigarEnd; ++cigarIter )
        {
            int L   = cigarIter->Length;
            char T  = cigarIter->Type;

            if (T == 'M')
            { // match or mismatch
                numNonIndelNucleotide += L;

                int ContiguousMismatch = 0;
                for (int l = 0; l< L; l++)
                {
                    char refbase = (char) referenceSequence.at(cgp + l);
                    char querybase = (char) alr.al.QueryBases.at(rp + l);

                    if ((toupper(querybase) == toupper(refbase)) ||
                            (tolower(querybase) == tolower(refbase)) ||
                            (querybase == '=')) // if querybase match
                    {
                        if (longestContiguousMismatch < ContiguousMismatch)
                            longestContiguousMismatch = ContiguousMismatch;

                        if (ContiguousMismatch > 0)
                            currentNumOfm++;

                        ContiguousMismatch = 0;
                    }
                    else
                    {
                        numMismatches++;
                        ContiguousMismatch++;
                    }
                }

                if (longestContiguousMismatch < ContiguousMismatch)
                    longestContiguousMismatch = ContiguousMismatch;

                if (ContiguousMismatch > 0)
                    currentNumOfm++;

                sp  += L;
                cgp += L;
                rp  += L;
            }
            else if (T == 'D')
            {
                currentNumOfIDS++;
                if (longestINDELlen < L)
                    longestINDELlen = L;

                sp  += L;  // update sample position
                cgp += L;
            }
            else if (T == 'I')
            {
                currentNumOfIDS++;
                if (longestINDELlen < L)
                    longestINDELlen = L;

                rp  += L;
            }
            else if (T == 'S')
            { // soft clip, clipped sequence present in the read not matching the reference
            // skip these bases in the read (because not included in the start position)
                //sp  += L;
                //cgp += L;
                currentNumOfIDS++;

                if (longestSoftclippingLength < L)
                    longestSoftclippingLength = L;


                rp  += L;
            }
            else if (T == 'H') { }// hard clip on the read, clipped sequence is not present in the read
            else if (T == 'N')
            { // skipped region in the reference not present in read, aka splice
                sp += L;
                cgp += L;
            }
        }
    }
    else // NOT-DUPLICATE/NOT-unmappedReads/ISPRIMARYALIGNMENT
    {
        // Initialize or Update the converage container
        if (m_CoverageByPositionInWindow.empty())
            m_CoverageByPositionInWindow.assign((winendpos-winstartpos+1),0);

        unsigned int startpos = 0;
        unsigned int endpos = m_CoverageByPositionInWindow.size();
        if (alr.al.Position > winstartpos)
            startpos = alr.al.Position-winstartpos;

        // NOTE from NEW BAMTOOLS: alignment.GetEndPosition() = StartPosition + readLength
        if (alr.al.GetEndPosition() < winendpos)
            endpos = alr.al.GetEndPosition() - winstartpos+1;
        else
            endpos = winendpos-winstartpos+1;

        for (unsigned i=startpos; i<endpos; i++)
            m_CoverageByPositionInWindow.at(i)++;
        
        bool stop = false;
        int cigarcount = 0;

        while (stop != true)
        {

            if ( (IsupdateMap == true) && (sp <= winendpos) )
            {
                //alr.currentReadPosition = rp;
                alr.cigarindex = cigarcount;
                alr.currentAlPosition = sp;
                //alr.currentGenomeSeqPosition = cgp;
            }
            else if ( (IsupdateMap == false) &&(sp <= winendpos) )
            {
                //(*Iter).currentReadPosition = rp;
                (*Iter).cigarindex = cigarcount;
                (*Iter).currentAlPosition = sp;
                //(*Iter).currentGenomeSeqPosition = cgp;                
            }

            int L   = alr.al.CigarData.at(cigarcount).Length;
            char T  = alr.al.CigarData.at(cigarcount).Type;

            if (alr.al.Name.find(str) != string::npos)
                cerr << alr.al.Name << "   LT " << L << T << endl;

                if (T == 'M')
                {
                    string snps;

                    bool updateLastsnp = false;
                    double MismatchpreviousPos;
                    double MismatchstartPos;
                    float accbasequality = 0;
                    int alleleLength = 0;

                    numNonIndelNucleotide += L;
                    int ContiguousMismatch = 0;

                    for (int l = 0; l < L; l++)
                    {
                        if (alr.al.Name.find(str) != string::npos)
                            cerr << referenceSequence.length() << " " << cgp << endl;
                        
                        char refbase = (char) referenceSequence.at(cgp + l);
                        char querybase = (char) alr.al.QueryBases.at(rp + l);

                        if (alr.al.Name.find(str) != string::npos)
                            cerr << "refbase " << refbase << " Q : " << querybase << endl;

                            if ((toupper(querybase) == toupper(refbase)) ||
                                    (tolower(querybase) == tolower(refbase)) ||
                                    (querybase == '=')) // if querybase match
                            {
                                // REALIGN
                                if (longestContiguousMismatch < ContiguousMismatch)
                                    longestContiguousMismatch = ContiguousMismatch;

                                if (ContiguousMismatch > 0)
                                    currentNumOfm++;

                                ContiguousMismatch = 0;
                            }
                            else
                            {
                                numMismatches++;
                                ContiguousMismatch++;
                                if ((winstartpos <= (sp + l)) && ((sp + l) <= winendpos))
                                {                                    
                                    float BaseQual = getBaseQualityASCII2float(alr.al.Qualities.at(rp + l));

                                    if (alr.al.Name.find(str) != string::npos)
                                        cerr << "Cigar=M " << winstartpos << ":" << winendpos << " ; " << (sp + l) << " " << alignment.Name << " : " << refbase << " " << querybase << " ;Q= " << BaseQual << endl;

                                    // REALIGN
                                    double newsp = sp + l;
                                    if (snps.empty())
                                    {
                                        snps = snps + querybase;
                                        MismatchpreviousPos = newsp;
                                        MismatchstartPos = newsp;
                                        accbasequality = BaseQual;
                                        alleleLength = 1;
                                    }
                                    else
                                    {
                                        if (MismatchpreviousPos == (newsp - 1))
                                        {
                                            snps = snps + querybase;

                                            MismatchpreviousPos = newsp;
                                            accbasequality += BaseQual;
                                            alleleLength++;
                                        }
                                        else
                                        {
                                            if (alignment.Name.find(str) !=string::npos )
                                                cerr << alignment.Name << " : " << MismatchstartPos << " accbasequality= " << accbasequality << " snps= " << snps << " bQual= " << alignment.Qualities.substr(rp + l, 1) << "R= " << refbase << " " << querybase <<endl;

                                            if ((accbasequality/(float)alleleLength) >= MinBaseQual)
                                            {
                                                MismatchNPosition MP(MismatchstartPos, snps);

                                                // UPdate Variant Position and COunts
                                                if (m_Position_MismatchCounts.empty())
                                                {
                                                    vector<int> tmp;
                                                    tmp.push_back(1);
                                                    tmp.push_back(snps.length());
                                                    m_Position_MismatchCounts[MP] = tmp;
                                                }
                                                else
                                                {
                                                    map<MismatchNPosition, vector<int> >::reverse_iterator rit = m_Position_MismatchCounts.rbegin();
                                                    if (MismatchstartPos > (*rit).first.POSITION)
                                                    {
                                                        vector<int> tmp;
                                                        tmp.push_back(1);
                                                        tmp.push_back(snps.length());
                                                        m_Position_MismatchCounts[MP] = tmp;
                                                    }
                                                    else
                                                    {
                                                        map<MismatchNPosition, vector<int> >::iterator it = m_Position_MismatchCounts.find(MP);
                                                        if (it != m_Position_MismatchCounts.end())
                                                        {
                                                            (*it).second.at(0)++;
                                                            (*it).second.at(1) += snps.length();
                                                        }
                                                        else
                                                        {
                                                            vector<int> tmp;
                                                            tmp.push_back(1);
                                                            tmp.push_back(snps.length());
                                                            m_Position_MismatchCounts[MP] = tmp;
                                                        }
                                                    }
                                                }
                                            }
                                            snps.clear();

                                            snps = snps + querybase;
                                            MismatchpreviousPos = newsp;
                                            MismatchstartPos = newsp;
                                            accbasequality = BaseQual;
                                            alleleLength = 1;
                                        }
                                    }
                                }
                            }
                        }

                //REALIGN : UPdate last snp entry
                // MultiMap Key: SNPs, Value=Start Position vector
                 if (alignment.Name.find(str) !=string::npos )
                     cerr << updateLastsnp << " snps = " << snps << endl;

                if ( (updateLastsnp == false) && (!snps.empty()) ) //( ((sp + l) >  winendpos) && (updateLastsnp == false) && (!snps.empty()) )
                {
                    updateLastsnp = true;

                    if (alignment.Name.find(str) !=string::npos )
                        cerr << "continue SNP " << MismatchstartPos << " accbasequality= " << accbasequality << " snps= " << snps << endl;

                    if ((accbasequality/(float)alleleLength) >= MinBaseQual)
                    {
                        MismatchNPosition MP(MismatchstartPos, snps);

                        // UPdate Variant Position and COunts
                        if (m_Position_MismatchCounts.empty())
                        {
                            vector<int> tmp;
                            tmp.push_back(1);
                            tmp.push_back(snps.length());

                            m_Position_MismatchCounts[MP] = tmp;
                        }
                        else
                        {
                            map<MismatchNPosition, vector<int> >::reverse_iterator rit = m_Position_MismatchCounts.rbegin();
                            if (MismatchstartPos > (*rit).first.POSITION)
                            {
                                vector<int> tmp;
                                tmp.push_back(1);
                                tmp.push_back(snps.length());

                                m_Position_MismatchCounts[MP] = tmp;
                            }
                            else
                            {
                                map<MismatchNPosition, vector<int> >::iterator it = m_Position_MismatchCounts.find(MP);

                                if (it != m_Position_MismatchCounts.end())
                                {
                                    (*it).second.at(0)++;
                                    (*it).second.at(1)+= snps.length();
                                }
                                else
                                {
                                    vector<int> tmp;
                                    tmp.push_back(1);
                                    tmp.push_back(snps.length());

                                    m_Position_MismatchCounts[MP] = tmp;
                                }
                            }
                        }
                    }
                    snps.clear();
                }


                if (longestContiguousMismatch < ContiguousMismatch)
                    longestContiguousMismatch = ContiguousMismatch;

                if (ContiguousMismatch > 0)
                    currentNumOfm++;

                sp += L;
                cgp += L;
                rp += L;

                if (alignment.Name.find(str) !=string::npos )
                    cerr << "DONE M " << endl;
                
                }

                else if (T == 'D')
                {
                    currentNumOfIDS++;

                    if (longestINDELlen < L)
                        longestINDELlen = L;

                    if ((winstartpos <= sp) && (sp <= winendpos))
                    {
                        // Prune based on the two left aand right bases quality
                        int offsetStart, offsetEnd;
                        if (cigarcount == 0)
                        {
                            offsetStart = 0;
                            offsetEnd = 2;
                        }
                        else if (cigarcount == (int) alr.al.CigarData.size() - 1)
                        {
                            offsetStart = -2;
                            offsetEnd = 0;
                        }
                        else
                        {
                            if ((rp + 1) <= ((int) alr.al.Qualities.length() - 1))
                                offsetEnd = 2;
                            else if (rp == ((int) alr.al.Qualities.length() - 1))
                                offsetEnd = 1;

                            if ((rp - 1) <= 0)
                                offsetStart = -1;
                            else
                                offsetStart = -2;
                        }

                        float sumLeftRightbasequal = 0;
                        for (int offset = offsetStart; offset < offsetEnd; offset++)
                        {
                            if (offset != 0)
                                sumLeftRightbasequal += getBaseQualityASCII2float(alr.al.Qualities.at(rp + offset));
                        }

                        if ((sumLeftRightbasequal / (offsetEnd - offsetStart)) >= MinBaseQual)
                        {
                            // REALIGN
                            double newsp = sp;
                            vector<string> seq;
                            seq.push_back(referenceSequence.substr(cgp, L));
                            AlleleNPosition AP(newsp, seq);

                            if (alr.al.Name.find(str) != string::npos)
                                cerr << " -Cigar=D " << alignment.Name << " : " << winendpos << " =? " << sp << endl;

                            // UPdate Variant Position and COunts
                            if (m_Position_IndelCounts.empty())
                            {
                                vector<int> tmp;
                                tmp.push_back(0);
                                tmp.push_back(1);
                                tmp.push_back(L);

                                m_Position_IndelCounts[AP] = tmp;
                            }
                            else
                            {
                                map<AlleleNPosition, vector<int> >::reverse_iterator rit = m_Position_IndelCounts.rbegin();
                                if ((AP.POSITION > (*rit).first.POSITION))//( (AP.POSITION > (*rit).first.POSITION) || ((AP.POSITION > (*rit).first.POSITION) && (AP.SEQ.compare((*rit).first.SEQ) != 0) ) )
                                {
                                    vector<int> tmp;
                                    tmp.push_back(0);
                                    tmp.push_back(1);
                                    tmp.push_back(L);

                                    m_Position_IndelCounts[AP] = tmp;
                                }
                                else
                                {
                                    map<AlleleNPosition, vector<int> >::iterator it = m_Position_IndelCounts.find(AP);

                                    if (alr.al.Name.find(str) != string::npos)
                                        cerr << "     -Cigar=D ; checking new delseq " << referenceSequence.substr(cgp, L) << endl;

                                    if (it != m_Position_IndelCounts.end())
                                    {
                                        vector<string>::const_iterator its = find((*it).first.SEQ.begin(), (*it).first.SEQ.end(), referenceSequence.substr(cgp, L));

                                        (*it).second.at(1)++;
                                        (*it).second.at(2) += L;

                                        if (its == (*it).first.SEQ.end())
                                        {
                                            vector<string> keyseq((*it).first.SEQ);
                                            //cerr << "           -Cigar=D ; (*it).first.SEQ.size()= " << (*it).first.SEQ.size() << " keyseq.size()= " << keyseq.size() << endl;
                                            keyseq.push_back(referenceSequence.substr(cgp, L));
                                            AlleleNPosition APKey(newsp, keyseq);
                                            vector<int> tmp((*it).second);

                                            m_Position_IndelCounts.erase(it);
                                            m_Position_IndelCounts[APKey] = tmp;
                                        }
                                    }
                                    else
                                    {
                                        vector<int> tmp;
                                        tmp.push_back(0);
                                        tmp.push_back(1);
                                        tmp.push_back(L);

                                        m_Position_IndelCounts[AP] = tmp;
                                    }
                                }
                            }
                        }
                    }
                    sp += L; // update sample position
                    cgp += L;

                    if (alr.al.Name.find(str) != string::npos)
                        cerr << " End D" << "  currentNumOfIDS= " << currentNumOfIDS << endl;
                }
                else if (T == 'I')
                {
                    currentNumOfIDS++;
                    if (longestINDELlen < L)
                        longestINDELlen = L;

                    if ((winstartpos <= sp) && (sp <= winendpos))
                    {
                        // Prune based on Insertion base quality (Include 2left and right bases ???)
                        int offsetStart, offsetEnd;

                        if (cigarcount == 0)
                        {
                            offsetStart = 0;
                            offsetEnd = L + 2;
                        }
                        else if (cigarcount == (int) alr.al.CigarData.size() - 1)
                        {
                            offsetStart = -2;
                            offsetEnd = L;
                        }
                        else
                        {
                            if ((rp + L + 1) <= ((int) alr.al.Qualities.length() - 1))
                                offsetEnd = L + 2;
                            else if ((rp + L) == ((int) alr.al.Qualities.length() - 1))
                                offsetEnd = L + 1;

                            if ((rp - 2) < 0)
                                offsetStart = -1;
                            else
                                offsetStart = -2;
                        }

                        float sumLeftRightbasequal = 0;
                        for (int offset = offsetStart; offset < offsetEnd; offset++)
                            sumLeftRightbasequal += getBaseQualityASCII2float(alr.al.Qualities.at(rp + offset));

                        if ((sumLeftRightbasequal / (float) (offsetEnd - offsetStart)) >= MinBaseQual)
                        {
                            // REALIGN
                            double newsp = sp + 0.5;

                            vector<string> seq;
                            seq.push_back(alr.al.QueryBases.substr(rp, L));
                            AlleleNPosition AP(newsp, seq);

                            if (alr.al.Name.find(str) != string::npos)
                                cerr << " -Cigar=I " << alignment.Name << " : " << winendpos << " =? " << sp+0.5 << endl;
                            

                            // UPdate Variant Position and COunts
                            if (m_Position_IndelCounts.empty())
                            {
                                vector<int> tmp;
                                tmp.push_back(1);
                                tmp.push_back(1);
                                tmp.push_back(L);

                                m_Position_IndelCounts[AP] = tmp;
                            }
                            else
                            {
                                map<AlleleNPosition, vector<int> >::reverse_iterator rit = m_Position_IndelCounts.rbegin();
                                if ((AP.POSITION > (*rit).first.POSITION)) //( (AP.POSITION > (*rit).first.POSITION) ) (AP.POSITION > (*rit).first.POSITION) ||  ( (AP.POSITION == (*rit).first.POSITION) && (AP.SEQ.compare((*rit).first.SEQ) != 0) ) )
                                {
                                    vector<int> tmp;
                                    tmp.push_back(1);
                                    tmp.push_back(1);
                                    tmp.push_back(L);

                                    m_Position_IndelCounts[AP] = tmp;
                                }
                                else
                                {
                                    map<AlleleNPosition, vector<int> >::iterator it = m_Position_IndelCounts.find(AP); // The idea is to seperate the indel length ???

                                    if (it != m_Position_IndelCounts.end())
                                    {
                                        vector<string>::const_iterator its = find((*it).first.SEQ.begin(), (*it).first.SEQ.end(), alr.al.QueryBases.substr(rp, L));

                                        (*it).second.at(1)++;
                                        (*it).second.at(2) += L;

                                        if (its == (*it).first.SEQ.end())
                                        {
                                            vector<string> keyseq((*it).first.SEQ);
                                            keyseq.push_back(alr.al.QueryBases.substr(rp, L));

                                            AlleleNPosition APKey(newsp, keyseq);
                                            vector<int> tmp((*it).second);

                                            m_Position_IndelCounts.erase(it);

                                            m_Position_IndelCounts[APKey] = tmp;
                                        }
                                    }
                                    else
                                    {
                                        vector<int> tmp;
                                        tmp.push_back(1);
                                        tmp.push_back(1);
                                        tmp.push_back(L);

                                        m_Position_IndelCounts[AP] = tmp;
                                    }
                                }
                            }
                        }
                    }
                    rp += L;

                    if (alr.al.Name.find(str) != string::npos)
                        cerr << " End I" << "  currentNumOfIDS= " << currentNumOfIDS << endl;
                }
                else if (T == 'S')
                { // soft clip, clipped sequence present in the read not matching the reference
                    // skip these bases in the read (because not included in the start position)

                    //REALIGN
                    currentNumOfIDS++;
                    if (longestSoftclippingLength < L)
                        longestSoftclippingLength = L;

                rp  += L;
                if (alr.al.Name.find(str) != string::npos)
                    cerr << " End S" << "  currentNumOfIDS= " << currentNumOfIDS << endl;
            }
            else if (T == 'H')
            { }// hard clip on the read, clipped sequence is not present in the read
            else if (T == 'N')
            { // skipped region in the reference not present in read, aka splice
                sp += L;
                cgp += L;
            }

            cigarcount++;

            if (cigarcount == (int) alr.al.CigarData.size())
            {
                stop = true;
            }

        } // end cigar iter loop
    }

    // Update the MAP
    if (IsupdateMap == true)
    {
        if (alr.al.Name.find(str) != string::npos)
            cerr << alignment.Name << "  numNonIndelNucleotide = " << numNonIndelNucleotide << " longestINDELlen = " << longestINDELlen << " currentNumOfIDS = " << currentNumOfIDS << " currentNumOfm = " << currentNumOfm << " longestContiguousMismatch = " << longestContiguousMismatch << endl;
        //// Clear vector (In case there are any data in it)
        alr.newCigars.clear();
        alr.mismatches.clear();
        alr.longestINDELlength.clear();
        alr.SWHaplotypePos.clear();
        alr.bestCigarIndex.clear();
        alr.costF.clear();
        alr.longestINDELlength.push_back(longestINDELlen);

        alr.NumOfIDS.push_back(currentNumOfIDS);
        alr.NumOfm.push_back(currentNumOfm);
        alr.longestContiguousMismatch.push_back(longestContiguousMismatch);
        alr.CigarSoftclippingLength = longestSoftclippingLength;

        alr.newCigars.push_back(alr.al.CigarData);
        alr.mismatches.push_back((numMismatches/numNonIndelNucleotide));

        AlGroups.push_back(alr);
    }
    else
    {
         if (alr.al.Name.find(str) != string::npos)
            cerr << " testing"  << (*Iter).al.Name << "  numNonIndelNucleotide = " << numNonIndelNucleotide << " longestINDELlen = " << longestINDELlen << " currentNumOfIDS = " << currentNumOfIDS << " currentNumOfm = " << currentNumOfm << " longestContiguousMismatch = " << longestContiguousMismatch << endl;

        // Clear vector (In case there are any data in it)
        (*Iter).newCigars.clear();
        (*Iter).mismatches.clear();
        (*Iter).longestINDELlength.clear();
        (*Iter).SWHaplotypePos.clear();
        (*Iter).bestCigarIndex.clear();
        (*Iter).costF.clear();
        (*Iter).longestINDELlength.push_back(longestINDELlen);

        (*Iter).NumOfIDS.push_back(currentNumOfIDS);
        (*Iter).NumOfm.push_back(currentNumOfm);
        (*Iter).longestContiguousMismatch.push_back(longestContiguousMismatch);
        (*Iter).CigarSoftclippingLength = longestSoftclippingLength;

        (*Iter).newCigars.push_back(alr.al.CigarData);
        (*Iter).mismatches.push_back((numMismatches/numNonIndelNucleotide));
    }

    //cerr << " DoneParse" << endl;
}


float RealignFunctionsClass::AlleleCostFunction(float allelecounts, float depth, bool isIndel) {
    float score = 0;
    if (isIndel == false) // mismatches
    {
        score = allelecounts / depth;
    } else {
        score = allelecounts / depth;
    }

    return score;
}


// ************************
// TESTING

bool RealignFunctionsClass::PruningByNaiveSelectionProcedureAndConstructHaplotypes2(int winstartpos, int winendpos, int32_t RefID, string RefIdstring, vector<int> minmaxRefSeqPos, FastaReference* reference)
{
    cerr << "#Start: Prune the Alleles & Construct Haplotypes ... " << endl;

    bool IsToRealign = false;
    m_isLowComplexityRegion = false;
    
    if (!m_Position_IndelCounts.empty())
    {
        IsToRealign = true;

        vector<MismatchNPosition> MismatchPositions;
        vector<AlleleNPosition> InsertionPositions;
        vector<AlleleNPosition> DeletionPositions;
        vector<float> MismatchFreq;
        vector<float> InsertionFreq;
        vector<float> DeletionFreq;

        // Calculate Frequency and        
        for (map<MismatchNPosition, vector<int> >::iterator it = m_Position_MismatchCounts.begin(); it != m_Position_MismatchCounts.end(); ++it)
        {
            unsigned int vecindex = (*it).first.POSITION - winstartpos;

            if (vecindex == m_CoverageByPositionInWindow.size())
                vecindex--;

            MismatchPositions.push_back((*it).first);
            float scorefreq = AlleleCostFunction((float) (*it).second.at(0), (float) m_CoverageByPositionInWindow.at(vecindex), false);
            MismatchFreq.push_back(scorefreq);

            //cerr.precision(10);
            //cerr << "  -MISMATCH (*it).first.POSITION = " << fixed  << (*it).first.POSITION << " ; vecindex= " << vecindex;
            //cerr <<" score = " << (*it).second.at(0) << " " <<  m_CoverageByPositionInWindow.at(vecindex) <<  " " << scorefreq << endl;
        }

        
        for (map<AlleleNPosition, vector<int> >::iterator it = m_Position_IndelCounts.begin(); it != m_Position_IndelCounts.end(); ++it)
        {
            unsigned int vecindex;
            if ((*it).second.at(0) == 1) // Insertion
            {
                vecindex = ((*it).first.POSITION - 0.5) - winstartpos;
                if (vecindex == m_CoverageByPositionInWindow.size())
                    vecindex--;

                InsertionPositions.push_back((*it).first);

                //cerr.precision(10);
                //cerr << "  -INSERTION floor((*it).first)= " << fixed << floor((*it).first.POSITION) << " ; vecindex= " << vecindex;
                
                float aveIndellen = (float) (*it).second.at(2) / (float) (*it).second.at(1);

                float weightedscorefreq = aveIndellen * AlleleCostFunction((float) (*it).second.at(1), (float) m_CoverageByPositionInWindow.at(vecindex), true);
                InsertionFreq.push_back(weightedscorefreq);

                
                //cerr << "aveIndellen " << aveIndellen << endl;
            }
            else if ((*it).second.at(0) == 0)
            {
                vecindex = (*it).first.POSITION - winstartpos;

                if (vecindex == m_CoverageByPositionInWindow.size())
                    vecindex--;

                DeletionPositions.push_back((*it).first);

                //cerr.precision(10);
                //cerr << "  -DELETION  (*it).first.POSITION= " << fixed << (*it).first.POSITION << " ; vecindex= " << fixed << vecindex;
                
                float aveIndellen = (float) (*it).second.at(2) / (float) (*it).second.at(1);
                float weightedscorefreq = aveIndellen * AlleleCostFunction((float) (*it).second.at(1), (float) m_CoverageByPositionInWindow.at(vecindex), true);
                DeletionFreq.push_back(weightedscorefreq);
                
                //cerr << "aveIndellen " << aveIndellen << "   (*it).first.SEQ.size() = " << (*it).first.SEQ.size() << endl;
            }
        }


        
        ////******************
        int numMismatches = 2; //3;
        int numDel = 2;
        int numIns = 2;

        vector<int> inslength, insrefpos, insreadpos;
        vector<string> InsSeq;
        vector<char> InsType;

        // Based on Original Reference Sequence
        string refseq = reference->getSubSequence(RefIdstring, minmaxRefSeqPos.at(0), 2*(minmaxRefSeqPos.at(1)-minmaxRefSeqPos.at(0)+1));
        SubseqInfo refinfo(RefID, refseq, minmaxRefSeqPos.at(0), inslength, insrefpos, insreadpos, InsSeq, InsType);
        m_Haplotypevec.push_back(refinfo);
       
        //// CHECK if WINDOW is homopolymerRun region
        string windowedrefseq = refseq.substr((winstartpos-minmaxRefSeqPos.at(0)),(winendpos-winstartpos+1));        
        m_isLowComplexityRegion = IsWindowInRepeatRegion(windowedrefseq);


        // Based on Mismatched Reference Sequence
        // Rank score and select limited Variants
        if (!MismatchFreq.empty()) {
            vector<int> MismatchPositionsIndex = GetIndexForSortedDecendingFloatVector(MismatchFreq);

            if ((int) MismatchPositionsIndex.size() < numMismatches)
                numMismatches = MismatchPositionsIndex.size();

            for (int i = 0; i < numMismatches; i++)                
            {
                string Mrefseq = refseq;
                Mrefseq.replace(MismatchPositions.at(MismatchPositionsIndex.at(i)).POSITION - minmaxRefSeqPos.at(0), MismatchPositions.at(MismatchPositionsIndex.at(i)).SEQ.length(), MismatchPositions.at(MismatchPositionsIndex.at(i)).SEQ);

                SubseqInfo refinfo(RefID, Mrefseq, minmaxRefSeqPos.at(0), inslength, insrefpos, insreadpos, InsSeq, InsType);
                m_Haplotypevec.push_back(refinfo);
            }
        }



        // Only reserve longest insertions and insertion that are different
        if (!InsertionFreq.empty()) {
            vector<int> IndelPositionsIndex = GetIndexForSortedDecendingFloatVector(InsertionFreq);

            if ((int) IndelPositionsIndex.size() < numIns)
                numIns = IndelPositionsIndex.size();

            for (int i = 0; i < numIns; i++) {
                inslength.clear();
                insrefpos.clear();
                insreadpos.clear();
                InsSeq.clear();
                InsType.clear();

                string refHI = refseq;
                string longestInsseq = *max_element(InsertionPositions.at(IndelPositionsIndex.at(i)).SEQ.begin(), InsertionPositions.at(IndelPositionsIndex.at(i)).SEQ.end(), stringLength());
                refHI.insert(InsertionPositions.at(IndelPositionsIndex.at(i)).POSITION - 0.5 - minmaxRefSeqPos.at(0), longestInsseq);

                inslength.push_back(longestInsseq.length());
                insrefpos.push_back(InsertionPositions.at(IndelPositionsIndex.at(i)).POSITION - 0.5);
                insreadpos.push_back(InsertionPositions.at(IndelPositionsIndex.at(i)).POSITION - 0.5 - minmaxRefSeqPos.at(0));
                InsSeq.push_back(longestInsseq);
                InsType.push_back('I');
                SubseqInfo hapinfoI(RefID, refHI, minmaxRefSeqPos.at(0), inslength, insrefpos, insreadpos, InsSeq, InsType);
                m_Haplotypevec.push_back(hapinfoI);

                for (int k = 0; k < (int) InsertionPositions.at(IndelPositionsIndex.at(i)).SEQ.size(); k++) {
                    inslength.clear();
                    insrefpos.clear();
                    insreadpos.clear();
                    InsSeq.clear();
                    InsType.clear();

                    if (longestInsseq.find(InsertionPositions.at(IndelPositionsIndex.at(i)).SEQ.at(k)) == string::npos) {
                        refHI = refseq;
                                              
                        string seq = InsertionPositions.at(IndelPositionsIndex.at(i)).SEQ.at(k);

                        refHI.insert(InsertionPositions.at(IndelPositionsIndex.at(i)).POSITION - 0.5 - minmaxRefSeqPos.at(0), seq);

                        //cerr.precision(10);
                        //cerr << " INSERTION: " << fixed << InsertionPositions.at(IndelPositionsIndex.at(i)).POSITION - 0.5 << " - " << minmaxRefSeqPos.at(0) << " ; seq= " << seq << endl;
                        inslength.push_back(seq.length());
                        insrefpos.push_back(InsertionPositions.at(IndelPositionsIndex.at(i)).POSITION - 0.5);
                        insreadpos.push_back(InsertionPositions.at(IndelPositionsIndex.at(i)).POSITION - 0.5 - minmaxRefSeqPos.at(0));
                        InsSeq.push_back(seq);
                        InsType.push_back('I');

                        SubseqInfo hapinfoI(RefID, refHI, minmaxRefSeqPos.at(0), inslength, insrefpos, insreadpos, InsSeq, InsType);
                        m_Haplotypevec.push_back(hapinfoI);

                    }
                }
            }
        }


        if (!DeletionFreq.empty())
        {
            vector<int> IndelPositionsIndex = GetIndexForSortedDecendingFloatVector(DeletionFreq);

            if ((int) IndelPositionsIndex.size() < numDel)
                numDel = IndelPositionsIndex.size();

            for (int i = 0; i < numDel; i++)
            {
                for (int k = 0; k < (int) DeletionPositions.at(IndelPositionsIndex.at(i)).SEQ.size(); k++)
                {
                    inslength.clear();
                    insrefpos.clear();
                    insreadpos.clear();
                    InsSeq.clear();
                    InsType.clear();

                    string refHI = refseq;

                    string seq = DeletionPositions.at(IndelPositionsIndex.at(i)).SEQ.at(k);
                    refHI.erase(DeletionPositions.at(IndelPositionsIndex.at(i)).POSITION - minmaxRefSeqPos.at(0), seq.length());
 
                    //cerr.precision(10);
                    //cerr << " -DELETION: " << fixed << DeletionPositions.at(IndelPositionsIndex.at(i)).POSITION << " - " << minmaxRefSeqPos.at(0) << " ; seq= " << seq << endl;

                    inslength.push_back(seq.length());
                    insrefpos.push_back(DeletionPositions.at(IndelPositionsIndex.at(i)).POSITION);
                    insreadpos.push_back(DeletionPositions.at(IndelPositionsIndex.at(i)).POSITION - minmaxRefSeqPos.at(0));
                    InsSeq.push_back(seq);
                    InsType.push_back('D');

                    SubseqInfo hapinfoI(RefID, refHI, minmaxRefSeqPos.at(0), inslength, insrefpos, insreadpos, InsSeq, InsType);
                    m_Haplotypevec.push_back(hapinfoI);
                }
            }
        }
    }

    cerr << "#Done: Prune the Alleles & Construct Haplotypes ... " << endl;

    return IsToRealign;
}

bool IsWindowInRepeatRegion(string sequence)
{
    bool isrepeatregion = false;

    char nprev, nnext;
    int lengthHrun = 0;
    //Check homopolymerruns
    for (int sp=0; sp<(int)sequence.length(); sp++)
    {
        if (sp==0)
        {
            nprev = (char) sequence.at(sp);
            lengthHrun++;
        }
        else
        {
            nnext = (char) sequence.at(sp);
            if  (nprev != nnext)
                lengthHrun = 1;
            else
                lengthHrun++;

            //cerr << nprev << " = " << nnext << " : "<< lengthHrun << endl;
            nprev = nnext;
        }

        if (lengthHrun > 5) //At least one nucleotide repeat 4 or more consecutively
        {
            isrepeatregion = true;
            break;
        }
    }

    //cheap n-nucleotide repeats
    int repeatlength = 2;
    int startPos = 0;
    int repeatunit = 3;
    while (isrepeatregion == false)
    {
        if ((startPos+(repeatlength-1)) < sequence.length() )
        {
            string cropSeq = sequence.substr(startPos,repeatlength);
            string dupcropSeq = cropSeq;
            for (int r=1 ;r<repeatunit;r++)
                dupcropSeq.append(cropSeq);

            // initialize
            // Optimized SLX SW Parameters
            float matchScore = 10.0;
            float mismatchScore = -9.0;
            float gapOpenPenalty = 15.0;
            float gapExtendPenalty = 6.66;

            int currentlongestINDELlen = 0;
            float currentfractionMismatches = 0;
            int currentNumOfIDS = 0;
            int currentNumOfm = 0;
            int currentlongestContiguousMismatch = 0;
            int currentlongestSlen = 0;
            int Totalreadlength = 0;

            float Bestscore = 0;
            unsigned int PosInSWPos;
            string Cigarstr;
            vector<CigarOp> newcigarinCigarOp, BamcigarinCigarOp;

            const char* haplotype = sequence.c_str();
            const char* query = dupcropSeq.c_str();
            const unsigned int haplotypeLen = strlen(haplotype);
            const unsigned int queryLen = strlen(query);

            //// create a new Smith-Waterman alignment object
            CSmithWatermanGotoh sw(matchScore, mismatchScore, gapOpenPenalty, gapExtendPenalty);
            sw.Align(PosInSWPos, Cigarstr, newcigarinCigarOp, BamcigarinCigarOp, Bestscore, currentlongestINDELlen, Totalreadlength, currentlongestSlen, currentfractionMismatches, currentNumOfIDS, currentNumOfm, currentlongestContiguousMismatch, haplotype, haplotypeLen, query, queryLen);

            //cerr << " Mi " << repeatlength << " " << dupcropSeq << " : " << Cigarstr << endl;
            
            if ( (Cigarstr.find("D") == string::npos) && (Cigarstr.find("I") == string::npos) && (Cigarstr.find("m") == string::npos) && (Cigarstr.find("S") == string::npos) )
            {
                isrepeatregion =true;
                break;
            }
            else
                startPos++;
        }
        else
        {
            startPos = 0;
            repeatlength++;     //Increment repeat length
        }

        if (repeatlength > 5)
            break;
    }

    cerr << " Repeat Region?  " << sequence << "  true? " << isrepeatregion << endl;

    return isrepeatregion;
}


// Original SMITH-WATERMAN WITH SmithWaterman Scoring: Allow more than 2 best selected haplotypes that support most of the best aligned reads with a threshold
 void RealignFunctionsClass::SelectHaplotypeCandidates_SmithWatermanBSv(vector<SalRealignInfo>& AlGroups, vector<int> minmaxRefSeqPos, bool SetLowComplexityRegionSWGapExt)
 {
    cerr << "#Start: Apply SW to realign all Alignments and Search for Best Haplotypes ..." << endl;


    // Optimized SLX SW Parameters
    float matchScore = 10.0;
    float mismatchScore = -9.0;
    float gapOpenPenalty = 15.0;
    float gapExtendPenalty = 6.66;

    if ( (m_isLowComplexityRegion==true) && (SetLowComplexityRegionSWGapExt==true) )
        gapExtendPenalty = 0;

    // ----------------------------------------------------------------
    // // realign reads with SW for unique haplotype and keep score
    // ----------------------------------------------------------------
    unsigned int readperfectmatchcount = 0;
    vector<float> HaplotypeBestMatch;
    HaplotypeBestMatch.assign(m_Haplotypevec.size(), 0);

    for (vector<SalRealignInfo>::iterator Iter = AlGroups.begin(); Iter != AlGroups.end(); ++Iter)
    {
        string str = "ZZZZZZZZZZZZZZZZ";
        //string str= "tL2_1941_2437_0:0:0_1:0:0_92";

        (*Iter).newCigars.clear();
        (*Iter).mismatches.clear();
        (*Iter).longestINDELlength.clear();
        (*Iter).NumOfIDS.clear();
        (*Iter).NumOfm.clear();
        (*Iter).longestContiguousMismatch.clear();
        (*Iter).SWHaplotypePos.clear();
        (*Iter).costF.clear();
        (*Iter).bestCigarIndex.clear();

        const char* query = (*Iter).al.QueryBases.c_str();

        int hapindex = 0;
        float prevBestscore = -1;        
        
        for (vector<SubseqInfo>::iterator haliter = m_Haplotypevec.begin(); haliter != m_Haplotypevec.end(); ++haliter)
        {
            // initialize
            int currentlongestINDELlen = 0;
            float currentfractionMismatches = 0;
            int currentNumOfIDS = 0;
            int currentNumOfm = 0;
            int currentlongestContiguousMismatch = 0;
            int currentlongestSlen = 0;
            int Totalreadlength = 0;

            float Bestscore = 0;

            unsigned int haplotypePosInSWPos;
            string newToRealignCigarstr;
            vector<CigarOp> newToRealigncigarinCigarOp, BamcigarinCigarOp;

            //cerr << "  (*Iter).costF.front() = " << (*Iter).costF.front() << "  Align " << endl;

            const char* haplotype = haliter->SEQ.c_str();

            const unsigned int haplotypeLen = strlen(haplotype);
            const unsigned int queryLen = strlen(query);

            //// create a new Smith-Waterman alignment object
            CSmithWatermanGotoh sw(matchScore, mismatchScore, gapOpenPenalty, gapExtendPenalty);
            sw.Align(haplotypePosInSWPos, newToRealignCigarstr, newToRealigncigarinCigarOp, BamcigarinCigarOp, Bestscore, Totalreadlength, currentlongestINDELlen, currentlongestSlen, currentfractionMismatches, currentNumOfIDS, currentNumOfm, currentlongestContiguousMismatch, haplotype, haplotypeLen, query, queryLen);

            if ((*Iter).al.Name.find(str) != string::npos)
            {
                cerr << "    CHECKHapSize " << m_Haplotypevec.size() << " " << hapindex << endl;
                cerr << "    CHECKINGSW: " << (*Iter).al.Name <<  " queryLen= " << queryLen << " ;swpos= " << haplotypePosInSWPos << " AlNewCigar " << newToRealignCigarstr << " size " << newToRealigncigarinCigarOp.size() << " longestindel " << currentlongestINDELlen;
                cerr << " mismatches " <<  currentfractionMismatches << endl;
            }

            // Update alignment struct
            (*Iter).newCigars.push_back(newToRealigncigarinCigarOp);
            (*Iter).mismatches.push_back(currentfractionMismatches);
            (*Iter).longestINDELlength.push_back(currentlongestINDELlen);
            (*Iter).NumOfIDS.push_back(currentNumOfIDS);
            (*Iter).NumOfm.push_back(currentNumOfm);
            (*Iter).longestContiguousMismatch.push_back(currentlongestContiguousMismatch);
            (*Iter).SWHaplotypePos.push_back(haplotypePosInSWPos);

            if (Totalreadlength > 2*(*Iter).al.Length)
                Bestscore -= (Totalreadlength-2*(*Iter).al.Length)*6.66;

            (*Iter).costF.push_back((-1)*Bestscore);

            // ***********************
            //// Store Best Scores
            // ***********************
            if (prevBestscore < 0)
            {
                prevBestscore = Bestscore;
                (*Iter).bestCigarIndex.push_back(hapindex);
            }
            else if (prevBestscore < Bestscore)
            {
                prevBestscore  = Bestscore;
                (*Iter).bestCigarIndex.clear();
                (*Iter).bestCigarIndex.push_back(hapindex);
            }
            else
                (*Iter).bestCigarIndex.push_back(hapindex);



            //// print
            if ((*Iter).al.Name.find(str) != string::npos)
            {
                if (!haliter->INSREFPOS.empty())
                {
                    for (unsigned int i = 0; i < haliter->INSREFPOS.size(); i++)
                        cerr << "  CHECKHAPfirst: " << hapindex << " " << haliter->INSREFPOS.at(i) << " " << haliter->INSREADPOS.at(i) << " " << haliter->INSLENGTH.at(i) << haliter->TYPE.at(i) << " " << haliter->INSSEQ.at(i) << endl;
                } else
                    cerr << "  CHECKHAPfirst: " << hapindex << " NONE " << endl;

                cerr << "    CHECKHapSize " << m_Haplotypevec.size() << " " << (*Iter).bestCigarIndex.size() << " " << hapindex << " Iter->bestCigarIndex.back() " << (*Iter).bestCigarIndex.back() << endl;
                cerr << "    CHECKINGSW: " << (*Iter).al.Name << " swpos= " << haplotypePosInSWPos << " AlNewCigar= " << newToRealignCigarstr << " size= " << newToRealigncigarinCigarOp.size() << " longestindel= " << (*Iter).longestINDELlength.back();
                cerr << " longest mismatches= " << (*Iter).mismatches.back() << " NumOfIDS= " << (*Iter).NumOfIDS.back() << " NumOfm= " << (*Iter).NumOfm.back() << " longestContiguousMismatch= " << (*Iter).longestContiguousMismatch.back();
                cerr << " cost: " << Bestscore << endl;
                //cerr << "         Hap " << haliter->SEQ << endl;
                //cerr << "         Que " << (*Iter).al.QueryBases << endl;
                cerr << " ----------------------------------------- " << endl;
                cerr << endl;
            }

            hapindex++;
        }

        // ------------------------------------------------------
        // Upcate total read count for Haplotype best aligned
        // ------------------------------------------------------
        float read = (float) 1 / (*Iter).bestCigarIndex.size();
        for (unsigned int i = 0; i < (*Iter).bestCigarIndex.size(); i++) {
            HaplotypeBestMatch.at((*Iter).bestCigarIndex.at(i)) += read;

            //if ( Iter->al.Name.find(str) !=string::npos )
            //    cerr << "   " << Iter->al.Name << " " << Iter->bestCigarIndex.size() << " " << Iter->bestCigarIndex.at(i) << endl;
        }
    }

    // ------------------------------------------------------
    // Find number of "Best" aligned cases
    // ------------------------------------------------------
    for (unsigned int i = 0; i < HaplotypeBestMatch.size(); i++)
    {
        if (i == 0)
            cerr << " RefSeq" << " ; TotalReadBestAligned = " << HaplotypeBestMatch.at(i) << endl;
        else
            cerr << " Hap #" << i << " ; TotalReadBestAligned = " << HaplotypeBestMatch.at(i) << endl;
    }

    // -------------------------------------------------------------------------------------------------------------------------------------
    //Chose haplotype by considering the counts. Possibly need to rethink this
    // -------------------------------------------------------------------------------------------------------------------------------------
    float ratio = 0.3;
    if ((HaplotypeBestMatch.size() > 1) && (HaplotypeBestMatch.at(0) > 0 || (readperfectmatchcount > 0)))
    {
        vector<float> duplicate(HaplotypeBestMatch.size() - 1);
        copy(HaplotypeBestMatch.begin() + 1, HaplotypeBestMatch.end(), duplicate.begin());
        m_ChosenBestHaplotypesIndex = GetIndexForSortedFloatVector(duplicate);

        float HighestReads = duplicate.at(m_ChosenBestHaplotypesIndex.back());

        vector<int>::iterator viter = m_ChosenBestHaplotypesIndex.begin();

        while (viter != m_ChosenBestHaplotypesIndex.end()) {
            if (duplicate.at(*viter) < (HighestReads * ratio))
                m_ChosenBestHaplotypesIndex.erase(viter);
            else
                viter++;
        }

        for (unsigned int i = 0; i < m_ChosenBestHaplotypesIndex.size(); i++)
            m_ChosenBestHaplotypesIndex.at(i)++;

        m_ChosenBestHaplotypesIndex.push_back(0);
    } 
    else
    {
        m_ChosenBestHaplotypesIndex = GetIndexForSortedFloatVector(HaplotypeBestMatch);

        float HighestReads = HaplotypeBestMatch.at(m_ChosenBestHaplotypesIndex.back());

        vector<int>::iterator viter = m_ChosenBestHaplotypesIndex.begin();

        while (viter != m_ChosenBestHaplotypesIndex.end())
        {
            if (HaplotypeBestMatch.at(*viter) < (HighestReads * ratio))
                m_ChosenBestHaplotypesIndex.erase(viter);
            else
                viter++;
        }
    }

    cerr << " Total Selected Haplotypes = " << m_ChosenBestHaplotypesIndex.size() << endl;

    HaplotypeBestMatch.clear();
    cerr << "#Done: Apply SW to realign all Alignments and Search for Best Haplotypes" << endl;
 }


// *********
void RealignFunctionsClass::AdjustCigarsWRTChosenMultipleHaplotypesAndPrepareAlignmentsTobeWrittenOut(vector<SalRealignInfo>& AlGroups, multimap<int, BamAlignment>& SortRealignedAlignmentsMultimap, 
        FastaReference* reference, map<int, string>& RefIDRedName, vector<int>& minmaxRefSeqPos, int nextwinstartpos, int nextwinendpos, float minbaseQ, int& TotalAlR, int& NumAlR, int ploidy)
{
    cerr << "#Start: Adjust Alignments' CIGARs ..." << endl;

    if (!AlGroups.empty())
    {
        vector<SalRealignInfo>::iterator Iter = AlGroups.begin();
        while (Iter != AlGroups.end()) {
            TotalAlR++;

            int swpos, swhappos;
            int startHapPos;
            string Hapseq;
            vector<int> insrefpos;
            vector<int> bestINSlen;
            vector<int> bestINSreadpos;
            vector<string> bestINSseq;
            vector<char> bestType;
            // SW bestCigar
            vector<CigarOp> swbestcigar;

            // Clear all vectors
            insrefpos.clear();
            bestINSlen.clear();
            bestINSreadpos.clear();
            bestINSseq.clear();
            bestType.clear();

            string str = "ZZZZZZZZZZZZZZZZZZ";
            //string str = "tL2";
           

            if ((*Iter).al.Name.find(str) != string::npos)
            {
                cerr << "     " << (*Iter).al.Name;
                cerr << " " << (*Iter).bestCigarIndex.size() << "  newCigars-Best.size()= ";
                cerr << (*Iter).newCigars[m_ChosenBestHaplotypesIndex.back()].size() << " ;  newCigars-2ndBest.size()= " << (*Iter).newCigars[m_ChosenBestHaplotypesIndex.at(m_ChosenBestHaplotypesIndex.size() - 2)].size();
                cerr << " longestINDELlength-Best " << (*Iter).longestINDELlength.at(m_ChosenBestHaplotypesIndex.back()) << " longestINDELlength-2ndBest ";
                cerr << (*Iter).longestINDELlength.at(m_ChosenBestHaplotypesIndex.at(m_ChosenBestHaplotypesIndex.size() - 2));
                cerr << " Mismatch-Best " << (*Iter).mismatches.at(m_ChosenBestHaplotypesIndex.back()) << " Mismatch-2ndBest " << (*Iter).mismatches.at(m_ChosenBestHaplotypesIndex.at(m_ChosenBestHaplotypesIndex.size() - 2));
                cerr << " Cost-Best " << (*Iter).costF.at(m_ChosenBestHaplotypesIndex.back()) << " Cost-2ndBest " << (*Iter).costF.at(m_ChosenBestHaplotypesIndex.at(m_ChosenBestHaplotypesIndex.size() - 2)) << endl;
            }


            if (m_ChosenBestHaplotypesIndex.size() == 1) {
                swbestcigar = (*Iter).newCigars[m_ChosenBestHaplotypesIndex.back()];
                swpos = (*Iter).SWHaplotypePos.at(m_ChosenBestHaplotypesIndex.back());
                bestINSlen = m_Haplotypevec.at(m_ChosenBestHaplotypesIndex.back()).INSLENGTH; // Haplotype's INDEL info
                bestINSreadpos = m_Haplotypevec.at(m_ChosenBestHaplotypesIndex.back()).INSREADPOS;
                bestINSseq = m_Haplotypevec.at(m_ChosenBestHaplotypesIndex.back()).INSSEQ;
                bestType = m_Haplotypevec.at(m_ChosenBestHaplotypesIndex.back()).TYPE;
                Hapseq = m_Haplotypevec.at(m_ChosenBestHaplotypesIndex.back()).SEQ;

                startHapPos = m_Haplotypevec.at(m_ChosenBestHaplotypesIndex.back()).REFSTARTPOS;
                insrefpos = m_Haplotypevec.at(m_ChosenBestHaplotypesIndex.back()).INSREFPOS;

                if ((*Iter).al.Name.find(str) != string::npos)
                    cerr << "   chosen1 " << m_ChosenBestHaplotypesIndex.back() << endl;
            }
            else {
                // Rank the cost
                int SelectedIndex = 0;
                float currentcost;
                unsigned int numHap = 0;

                SelectedIndex = m_ChosenBestHaplotypesIndex.back();
                currentcost = (*Iter).costF.at(SelectedIndex);

                if ((*Iter).al.Name.find(str) != string::npos)
                    cerr << "Test 1" << " : " << SelectedIndex << " ; "<< currentcost << endl;

                if (ploidy == 2)
                    numHap = 3;
                else
                    numHap = m_ChosenBestHaplotypesIndex.size() + 1;

                for (unsigned int i = 2; i < numHap; i++) {
                    if (currentcost > (*Iter).costF.at(m_ChosenBestHaplotypesIndex.at(m_ChosenBestHaplotypesIndex.size() - i))) {
                        SelectedIndex = m_ChosenBestHaplotypesIndex.at(m_ChosenBestHaplotypesIndex.size() - i);
                        currentcost = (*Iter).costF.at(SelectedIndex);

                        if ((*Iter).al.Name.find(str) != string::npos)
                            cerr << "     Test " << i << " : " << SelectedIndex << " ; " << currentcost << endl;
                    }
                }

                swbestcigar = (*Iter).newCigars[SelectedIndex];
                swpos = (*Iter).SWHaplotypePos.at(SelectedIndex);
                bestINSlen = m_Haplotypevec.at(SelectedIndex).INSLENGTH;
                bestINSreadpos = m_Haplotypevec.at(SelectedIndex).INSREADPOS;
                bestINSseq = m_Haplotypevec.at(SelectedIndex).INSSEQ;
                bestType = m_Haplotypevec.at(SelectedIndex).TYPE;
                Hapseq = m_Haplotypevec.at(SelectedIndex).SEQ;

                startHapPos = m_Haplotypevec.at(SelectedIndex).REFSTARTPOS;
                insrefpos = m_Haplotypevec.at(SelectedIndex).INSREFPOS;

                if ((*Iter).al.Name.find(str) != string::npos)
                    cerr << "   chosen2 " << SelectedIndex << endl;
            }


            // ---------------------------------------------------------------------------------------------------------------------
            // ----------------------------------------------------------------------------------------------------
            // IF the read selected the REF Haplotype (Basically no INDELs Haplotype). No need adjustment
            // ----------------------------------------------------------------------------------------------------
            if (bestINSlen.empty()) // NO-indels Haplotype, skip this iteration
            {
                //// Update new cigar in BamAlignment (adjust the cigars
                vector<CigarOp> TobeBAMcigar;
                for (vector<CigarOp>::const_iterator cigarIter = swbestcigar.begin(); cigarIter != swbestcigar.end(); ++cigarIter) {
                    unsigned int L = cigarIter->Length;
                    char T = cigarIter->Type;

                    if (T != 'm')
                    {
                        if (TobeBAMcigar.empty() || (TobeBAMcigar.back().Type != T))
                        {
                            CigarOp op;
                            op.Length = L;
                            op.Type = T;
                            TobeBAMcigar.push_back(op);
                        }
                        else
                            TobeBAMcigar.back().Length += L;
                    }
                    else //"m"
                    {
                        if (TobeBAMcigar.empty() || (TobeBAMcigar.back().Type != 'M'))
                        {
                            CigarOp op;
                            op.Length = L;
                            op.Type = 'M';
                            TobeBAMcigar.push_back(op);
                        }
                        else
                            TobeBAMcigar.back().Length += L;
                    }
                }


                (*Iter).al.CigarData = TobeBAMcigar;

                if ((*Iter).al.Name.find(str) != string::npos)
                {
                    stringstream Newcigar;
                    for (vector<CigarOp>::const_iterator cigarIter = TobeBAMcigar.begin(); cigarIter != TobeBAMcigar.end(); ++cigarIter)
                        Newcigar << (*cigarIter).Length << (*cigarIter).Type;

                    cerr << "   C " << (*Iter).al.Name << " NewCigar " << " C: " << Newcigar.str() << endl;
                }

                //cerr << " TEST: swpos = " << swpos << " startHapPos= " << startHapPos << endl;
                (*Iter).al.Position = swpos + startHapPos;

                if ((*Iter).al.GetEndPosition() < nextwinstartpos)
                {
                    if ((*Iter).al.Name.find(str) != string::npos)
                        cerr << "  CToWrite: " << (*Iter).al.Name << " ; " << (*Iter).al.Position << " < " << nextwinstartpos << " : " << (*Iter).al.GetEndPosition() << endl;

                    SortRealignedAlignmentsMultimap.insert(pair<int, BamAlignment > ((*Iter).al.Position, (*Iter).al));
                    AlGroups.erase(Iter);

                    //cerr << "  ToWrite: DONE " << endl;
                } 
                else
                {
                    if ((*Iter).al.Name.find(str) != string::npos)
                        cerr << "  CToKEEP: " << (*Iter).al.Name << " ; " << (*Iter).al.Position << " < " << nextwinstartpos << " : " << (*Iter).al.GetEndPosition() << endl;

                    string referenceSequence = reference->getSubSequence(RefIDRedName[(*Iter).al.RefID], (*Iter).al.Position, 2*(*Iter).al.Length);
                    SalRealignInfo talr = (*Iter);
                    ParseAlignmentsAndExtractVariantsByBaseQualv(AlGroups, talr, Iter, (*Iter).al, referenceSequence, minmaxRefSeqPos, nextwinstartpos, nextwinendpos, minbaseQ, false);

                    ++Iter; //Increment iterator
                }
                continue;
            }

            //cerr << " N " << (*Iter).al.Name << endl;

            //------------------------------------------------------------------------------------------------------------------------------------------------
            // ------------------------------------------------------------------------------------------------
            // Adjust READ Start Position in BAM file
            // ------------------------------------------------------------------------------------------------
            //cerr << "   Adjust alignment start position " << endl;
            vector<CigarOp> TobeBAMcigar;
            int HapINDELlen;
            char HapINDELtype;
            string HaPINDELseq;
            int startdel = 0;
            int InsCount = -1;
            int qpos = 0;
            int EditDistance = 0;
            //DelOffset       = 0;

            swhappos = swpos + startHapPos;
            unsigned int hapCigarIndex;

            // Find the alignment genome position
            for (unsigned int i = 0; i < bestINSreadpos.size(); i++)
            {
                if (swpos > bestINSreadpos.at(i)) {
                    if (bestType.at(i) == 'I') {
                        int insend = bestINSreadpos.at(i) + bestINSlen.at(i);

                        if (swpos > insend)
                            swhappos -= bestINSlen.at(i);
                        else if (swpos <= insend) // Start with an insertion
                        {
                            // Clear all vectors
                            insrefpos.clear();
                            bestINSlen.clear();
                            bestINSreadpos.clear();
                            bestINSseq.clear();
                            bestType.clear();
                            swbestcigar.clear();

                            swbestcigar = Iter->newCigars[0];
                            swpos = Iter->SWHaplotypePos.at(0);
                            bestINSlen = m_Haplotypevec.at(0).INSLENGTH;
                            bestINSreadpos = m_Haplotypevec.at(0).INSREADPOS;
                            bestINSseq = m_Haplotypevec.at(0).INSSEQ;
                            bestType = m_Haplotypevec.at(0).TYPE;
                            Hapseq = m_Haplotypevec.at(0).SEQ;
                            startHapPos = m_Haplotypevec.at(0).REFSTARTPOS;
                            insrefpos = m_Haplotypevec.at(0).INSREFPOS;

                            swhappos = swpos + startHapPos;
                            break;
                        }
                    } else if (bestType.at(i) == 'D')
                        swhappos += bestINSlen.at(i);

                    hapCigarIndex = i;
                }
                else {
                    hapCigarIndex = i;
                    break;
                }
            }

            // Print to check
            if ((*Iter).al.Name.find(str) != string::npos)
            {
                cerr << "  " << (*Iter).al.Name << " bestINSreadpos.size() " << bestINSreadpos.size() << " beforePos " << (*Iter).al.Position << " swpos " << swpos << " startHapPos " << startHapPos << " AdjustPos " << swhappos << endl;
                cerr << "     swbestcigar.size()= " << swbestcigar.size() << endl;
                for (unsigned int i = 0; i < bestINSreadpos.size(); i++)
                {
                    cerr << "bestINSreadpos.at " << i << " " << bestINSreadpos.at(i) << " " << bestType.at(i) << "  " << bestINSseq.at(i) << " " << bestINSlen.at(i) << endl;
                    cerr << "insrefpos.at " << i << " " << insrefpos.at(i) << endl;
                }
                cerr << Hapseq << endl;
                cerr << hapCigarIndex << endl;
            }


            for (vector<CigarOp>::const_iterator cigarIter = swbestcigar.begin(); cigarIter != swbestcigar.end(); ++cigarIter) {
                unsigned int L = cigarIter->Length;
                char T = cigarIter->Type;
                unsigned int l = 0;

                while (l < L) {
                    // Check Haplotype's base first
                    if ((!bestINSlen.empty()) && (hapCigarIndex < bestINSlen.size())) {
                        if (bestINSreadpos.at(hapCigarIndex) == swpos) {
                            HapINDELlen = bestINSlen.at(hapCigarIndex);
                            if (bestType.at(hapCigarIndex) == 'D') {
                                HapINDELtype = 'D';
                                HaPINDELseq = bestINSseq.at(hapCigarIndex);

                                if (startdel >= HapINDELlen)
                                    startdel = 0;

                                if ((*Iter).al.Name.find(str) != string::npos)
                                    cerr << " D startdel= " << startdel << "  HapINDELlen= " << HapINDELlen << " HaPINDELseq= " << HaPINDELseq << endl;
                            }
                            else if (bestType.at(hapCigarIndex) == 'I') {
                                HapINDELtype = 'I';
                                InsCount = 1;
                                //HaPINDELseq = '';
                            }
                        }
                        else
                        {
                            if (InsCount < 0)
                            {
                                HapINDELlen = 1;
                                HapINDELtype = 'M';
                            }
                        }
                    }
                    else {
                        if (InsCount < 0) {
                            HapINDELlen = 1;
                            HapINDELtype = 'M';
                        }
                    }


                    if ((*Iter).al.Name.find(str) != string::npos)
                    {
                        cerr << "   " << (*Iter).al.Name << " " << qpos << " " << L << T << " " << l << " " << (*Iter).al.QueryBases.substr(qpos, 1) << " : swpos " << swpos;
                        cerr << " HapINDELlen = " << HapINDELlen << " ; HapINDELtype = " << HapINDELtype << " " << Hapseq.substr(swpos, 1);
                    }


                    // Adjust Cigar: Hapbase == D
                    if (HapINDELtype == 'D')
                    {
                        CigarOp op;
                        if ((T == 'M') || (T == 'D') | (T == 'm'))//)
                        {
                            if (TobeBAMcigar.empty() || (TobeBAMcigar.back().Type != 'D')) {
                                op.Length = HapINDELlen - startdel;
                                op.Type = 'D';
                                TobeBAMcigar.push_back(op);
                            }
                            else if (TobeBAMcigar.back().Type == 'D')
                                TobeBAMcigar.back().Length += HapINDELlen - startdel;

                            if ((*Iter).al.Name.find(str) != string::npos)
                                cerr << " NewC " << (HapINDELlen - startdel) << 'D' << endl;

                            EditDistance += HapINDELlen - startdel;
                            startdel = HapINDELlen;
                        }
                        else if (T == 'S') //Check if the seq in S
                        {
                            string qstr = (*Iter).al.QueryBases.substr(qpos, 1);
                            size_t foundSubstrPos = HaPINDELseq.substr(startdel).find(qstr);
                            //
                            if ((foundSubstrPos > 0) && (foundSubstrPos != string::npos))
                            {
                                if (TobeBAMcigar.empty() || (TobeBAMcigar.back().Type != 'D'))
                                {
                                    op.Length = foundSubstrPos + 1;
                                    op.Type = 'D';
                                    TobeBAMcigar.push_back(op);
                                }
                                else if (TobeBAMcigar.back().Type == 'D')
                                    TobeBAMcigar.back().Length += foundSubstrPos + 1;

                                EditDistance += foundSubstrPos + 1;
                            } 
                            else if (foundSubstrPos == string::npos)
                            {
                                if (TobeBAMcigar.empty() || (TobeBAMcigar.back().Type != 'D'))
                                {
                                    op.Length = HapINDELlen - startdel;
                                    op.Type = 'D';
                                    TobeBAMcigar.push_back(op);
                                }
                                else if (TobeBAMcigar.back().Type == 'D')
                                    TobeBAMcigar.back().Length += HapINDELlen - startdel;

                                EditDistance += HapINDELlen - startdel;

                                if ((*Iter).al.Name.find(str) != string::npos)
                                    cerr << " NewC " << 'D' << endl;
                            }

                            //
                            if ((L < (HapINDELlen - foundSubstrPos - startdel)) && (foundSubstrPos != string::npos))
                            {
                                if (TobeBAMcigar.empty() || (TobeBAMcigar.back().Type != 'M')) {
                                    op.Length = L;
                                    op.Type = 'M';
                                    TobeBAMcigar.push_back(op);
                                }
                                else if (TobeBAMcigar.back().Type == 'M')
                                    TobeBAMcigar.back().Length += L;

                                if ((*Iter).al.Name.find(str) != string::npos)
                                    cerr << " NewC " << 'M' << endl;
                            } 
                            else if ((L > (HapINDELlen - foundSubstrPos - startdel)) && (foundSubstrPos != string::npos))
                            {
                                if (TobeBAMcigar.empty() || (TobeBAMcigar.back().Type != 'M'))
                                {
                                    op.Length = (HapINDELlen - foundSubstrPos - startdel);
                                    op.Type = 'M';
                                    TobeBAMcigar.push_back(op);

                                    //update Editdistance !!!
                                }
                                else if (TobeBAMcigar.back().Type == 'M')
                                    TobeBAMcigar.back().Length += (HapINDELlen - foundSubstrPos - startdel);

                                if ((*Iter).al.Name.find(str) != string::npos)
                                    cerr << " NewC " << 'M' << endl;

                                //update 'S'
                                op.Length = L - (HapINDELlen - foundSubstrPos - startdel);
                                op.Type = 'S';
                                TobeBAMcigar.push_back(op);
                            } 
                            else if (foundSubstrPos == string::npos)
                            {
                                if (TobeBAMcigar.empty() || (TobeBAMcigar.back().Type != 'S'))
                                {
                                    op.Length = L;
                                    op.Type = 'S';
                                    TobeBAMcigar.push_back(op);
                                }
                                else if (TobeBAMcigar.back().Type == 'S')
                                    TobeBAMcigar.back().Length += L;
                            }

                            l = L;
                            qpos += L;
                            startdel = HapINDELlen;

                           if ((*Iter).al.Name.find(str) != string::npos)
                                cerr << " NewC " << 'S' << endl;
                        }
                        else if ((T == 'I'))//|| (T == 'm') )
                        {
                            string qstr = (*Iter).al.QueryBases.substr(qpos, L);
                            size_t foundSubstrPos = HaPINDELseq.substr(startdel).find(qstr);
                            //cerr << "  qstr= " << qstr << " foundSubstrPos= " << foundSubstrPos << " L= " << L << " HapINDELlen= " << HapINDELlen << " startdel=" << startdel << endl;

                            
                            if ((foundSubstrPos > 0) && (foundSubstrPos != string::npos)) // Found ???
                            {
                                if (TobeBAMcigar.empty() || (TobeBAMcigar.back().Type != 'D'))
                                {
                                    op.Length = foundSubstrPos + 1;
                                    op.Type = 'D';
                                    TobeBAMcigar.push_back(op);
                                }
                                else if (TobeBAMcigar.back().Type == 'D')
                                    TobeBAMcigar.back().Length += foundSubstrPos + 1;

                                if ((*Iter).al.Name.find(str) != string::npos)
                                    cerr << " NewC " << 'D' << endl;

                               EditDistance += foundSubstrPos + 1;
                            }
                            else if (foundSubstrPos == string::npos)  // Not found
                            {
                                if (TobeBAMcigar.empty() || (TobeBAMcigar.back().Type != 'D'))
                                {
                                    op.Length = HapINDELlen - startdel;
                                    op.Type = 'D';
                                    TobeBAMcigar.push_back(op);
                                }
                                else if (TobeBAMcigar.back().Type == 'D')
                                    TobeBAMcigar.back().Length += HapINDELlen - startdel;

                                if ((*Iter).al.Name.find(str) != string::npos)
                                    cerr << " NewC " << 'D' << endl;

                                EditDistance += HapINDELlen - startdel;
                                startdel = HapINDELlen;
                            }
                            
                            if ((L <= (HapINDELlen - foundSubstrPos - startdel)) && (foundSubstrPos != string::npos)) // FONUD ??
                            {
                                if (TobeBAMcigar.empty() || (TobeBAMcigar.back().Type != 'M'))
                                {
                                    op.Length = L;
                                    op.Type = 'M';
                                    TobeBAMcigar.push_back(op);
                                }
                                else if (TobeBAMcigar.back().Type == 'M')
                                    TobeBAMcigar.back().Length += L;

                                if ((*Iter).al.Name.find(str) != string::npos)
                                    cerr << " NewC " << 'M' << endl;

                                l = L;
                                qpos += L;
                                startdel += L;
                            } 
                            else if ((L > (HapINDELlen - foundSubstrPos - startdel)) && (foundSubstrPos != string::npos))
                            {
                                if (TobeBAMcigar.empty() || (TobeBAMcigar.back().Type != 'M'))
                                {
                                    op.Length = (HapINDELlen - foundSubstrPos - startdel);
                                    op.Type = 'M';
                                    TobeBAMcigar.push_back(op);

                                    //update Editdistance !!!
                                }
                                else if (TobeBAMcigar.back().Type == 'M')
                                    TobeBAMcigar.back().Length += (HapINDELlen - foundSubstrPos - startdel);

                                if ((*Iter).al.Name.find(str) != string::npos)
                                    cerr << " NewC " << 'M' << endl;

                                //update 'I'
                                if (T == 'I')
                                {
                                    op.Length = L - (HapINDELlen - foundSubstrPos - startdel);
                                    op.Type = 'I';
                                    TobeBAMcigar.push_back(op);

                                    l = L;
                                    qpos += L;

                                    if ((*Iter).al.Name.find(str) != string::npos)
                                        cerr << " NewC " << "I" << endl;
                                }
                                else // (T == 'm')
                                {
                                    l += (HapINDELlen - foundSubstrPos - startdel);
                                    qpos += (HapINDELlen - foundSubstrPos - startdel);
                                }

                                startdel = HapINDELlen;


                            } 
                            else if ((foundSubstrPos == string::npos) && (T == 'I'))
                            {
                                if (TobeBAMcigar.empty() || (TobeBAMcigar.back().Type != 'I'))
                                {
                                    op.Length = L;
                                    op.Type = 'I';
                                    TobeBAMcigar.push_back(op);
                                }
                                else if (TobeBAMcigar.back().Type == 'I')
                                    TobeBAMcigar.back().Length += L;

                                l = L;
                                qpos += L;

                                if ((*Iter).al.Name.find(str) != string::npos)
                                    cerr << " NewC " << "I" << endl;
                            }

                            
                        }

                        if (startdel >= HapINDELlen)
                            hapCigarIndex++;
                    }



                    // Adjust Cigar: Hapbase != D
                    if (HapINDELtype == 'I')
                    {
                        CigarOp op;

                        if ((T == 'M') || (T == 'm')) // // Include both Matches and Mismataches, ignore if T == 'D'
                        {
                            if (TobeBAMcigar.empty() || (TobeBAMcigar.back().Type != 'I'))
                            {
                                op.Length = 1;
                                op.Type = 'I';
                                TobeBAMcigar.push_back(op);
                            }
                            else if (TobeBAMcigar.back().Type == 'I')
                                TobeBAMcigar.back().Length++;

                            EditDistance++;
                            InsCount++;
                            swpos++;
                            qpos++;

                            if ((*Iter).al.Name.find(str) != string::npos)
                                cerr << " NewC " << 'I' << endl;

                        } 
                        else if (T == 'S')
                        {
                            if (TobeBAMcigar.empty() || (TobeBAMcigar.back().Type != 'S'))
                            {
                                op.Length = 1;
                                op.Type = 'S';
                                TobeBAMcigar.push_back(op);
                            }
                            else if (TobeBAMcigar.back().Type == 'S')
                                TobeBAMcigar.back().Length++;

                            qpos++;

                            if ((*Iter).al.Name.find(str) != string::npos)
                                cerr << " NewC " << 'S' << endl;

                        } 
                        else if (T == 'I')
                        {
                            if (TobeBAMcigar.empty() || (TobeBAMcigar.back().Type != 'I'))
                            {
                                op.Length = 1;
                                op.Type = 'I';
                                TobeBAMcigar.push_back(op);
                            }
                            else if (TobeBAMcigar.back().Type == 'I')
                                TobeBAMcigar.back().Length++;

                            EditDistance++;
                            qpos++;

                            if ((*Iter).al.Name.find(str) != string::npos)
                                cerr << " NewC " << 'I' << endl;
                        }
                        else if (T == 'D')
                        {
                            InsCount++;
                            swpos++;

                            if ((*Iter).al.Name.find(str) != string::npos)
                                cerr << " NewC " << '-' << endl;
                        }
                        l++;

                        //if (Iter->al.Name.find(str) !=string::npos )
                        //    cerr << "    " << InsCount << endl;

                        if (InsCount > HapINDELlen) {
                            hapCigarIndex++;
                            InsCount = -1;
                        }
                    }
                    else if (HapINDELtype == 'M')
                    {
                        CigarOp op;
                        if ((T == 'M') || (T == 'm')) // Include both Matches and Mismataches
                        {
                            if (TobeBAMcigar.empty() || (TobeBAMcigar.back().Type != 'M')) {
                                op.Length = 1;
                                op.Type = 'M';
                                TobeBAMcigar.push_back(op);
                            } else if (TobeBAMcigar.back().Type == 'M')
                                TobeBAMcigar.back().Length++;

                            if (T == 'm')
                                EditDistance++;

                            swpos++;
                            qpos++;

                            if ((*Iter).al.Name.find(str) != string::npos)
                                cerr << " NewC " << 'M' << endl;
                        }
                        if (T == 'D')
                        {
                            if (TobeBAMcigar.empty() || (TobeBAMcigar.back().Type != T)) {
                                op.Length = 1;
                                op.Type = 'D';
                                TobeBAMcigar.push_back(op);
                            } else if (TobeBAMcigar.back().Type == T)
                                TobeBAMcigar.back().Length++;

                            EditDistance++;
                            swpos++;

                            if ((*Iter).al.Name.find(str) != string::npos)
                                cerr << " NewC " << 'D' << endl;
                        }
                        else if (T == 'S')
                        {
                            if (TobeBAMcigar.empty() || (TobeBAMcigar.back().Type != 'S')) {
                                op.Length = 1;
                                op.Type = 'S';
                                TobeBAMcigar.push_back(op);
                            } else if (TobeBAMcigar.back().Type == 'S')
                                TobeBAMcigar.back().Length++;

                            qpos++;

                            if ((*Iter).al.Name.find(str) != string::npos)
                                cerr << " NewC " << 'S' << endl;
                        }
                        else if (T == 'I')
                        {
                            if (TobeBAMcigar.empty() || (TobeBAMcigar.back().Type != T)) {
                                op.Length = 1;
                                op.Type = 'I';
                                TobeBAMcigar.push_back(op);
                            } else if (TobeBAMcigar.back().Type == T)
                                TobeBAMcigar.back().Length++;

                            EditDistance++;
                            qpos++;

                            if ((*Iter).al.Name.find(str) != string::npos)
                                cerr << " NewC " << 'I' << endl;
                        }
                        l++;
                    }
                }
            }

            //------------------------------------------------------------------------------------------------
            // Check if first type is a deletion
            if (TobeBAMcigar.front().Type == 'D')
            {
                swhappos += TobeBAMcigar.front().Length;
                TobeBAMcigar.erase(TobeBAMcigar.begin());
            }

            ////------------------------------------------------------------------
            //// Check of this alignment is realigned (CIGAR or StartPos changes)
            stringstream Originalcigar, Newcigar;
            for (vector<CigarOp>::const_iterator cigarIter = (*Iter).al.CigarData.begin(); cigarIter != (*Iter).al.CigarData.end(); ++cigarIter)
                Originalcigar << (*cigarIter).Length << (*cigarIter).Type;


            for (vector<CigarOp>::const_iterator cigarIter = TobeBAMcigar.begin(); cigarIter != TobeBAMcigar.end(); ++cigarIter)
                Newcigar << (*cigarIter).Length << (*cigarIter).Type;

            //cerr << "   N " << Iter->al.Name << " NewCigar " << " C: " << Newcigar.str() << endl;

            string OriginalcigarStr = Originalcigar.str();
            string NewcigarStr = Newcigar.str();

            if (((*Iter).al.Position != swhappos) || (OriginalcigarStr.compare(NewcigarStr) != 0))
                NumAlR++;

            //if (((*Iter).al.Position != swhappos)) //|| (cigarstr1.compare(cigarstr2) != 0) )
            //{
            //    NumAlR++;
            //    if ((*Iter).al.Name.find(str) != string::npos)
            //        cerr << "   " << Iter->al.Name << " O: " << (*Iter).al.Position << " C: " << swhappos << endl;
            //}
            //else if (cigarstr1.compare(cigarstr2) != 0)
            //    NumAlR++;

            ////------------------------------------------------------------------------------------------------------------------------------------------------
            ////------------------------------------------------------------
            //// FOR other HAPLOTYPES: Update new cigar in BamAlignment
            ////------------------------------------------------------------
            (*Iter).al.CigarData = TobeBAMcigar;
            (*Iter).al.Position = swhappos;

            //If CIGAR is different
            if (OriginalcigarStr.compare(NewcigarStr) != 0)
            {
                //UPdate BAM TAG:TYPE:VALUE
                (*Iter).al.AddTag("OC", "Z", OriginalcigarStr);
                (*Iter).al.EditTag("NM", "i", EditDistance); //Assume Levenshtein Distance
                (*Iter).HasRealign = true;
            }

            if ((*Iter).al.Position != swhappos) //Original position is different from new
            {
                //UPDATE BAM TAG:TYPE:VALUE
                (*Iter).al.AddTag("OP", "i", swhappos); //The Position is "ZERO-base" in the bAM BUT in the SAM format the position is 1-based leftmost mapping POSition ?
                (*Iter).HasRealign = true;
            }

            //****************
            // Erase alignment s
            if ((*Iter).al.GetEndPosition() < nextwinstartpos)
            {
                if ((*Iter).al.Name.find(str) != string::npos)
                    cerr << "  ToWrite: " << (*Iter).al.Name << " ; " << (*Iter).al.Position << " < " << nextwinstartpos << " : " << (*Iter).al.GetEndPosition() << endl;

                SortRealignedAlignmentsMultimap.insert(pair<int, BamAlignment > ((*Iter).al.Position, (*Iter).al));
                AlGroups.erase(Iter);
                //cerr << "  ToWrite: DONE " << endl;
            }
            else
            {
                if ((*Iter).al.Name.find(str) != string::npos)
                    cerr << "  ToKEEP: " << (*Iter).al.Name << " ; " << (*Iter).al.Position << " < " << nextwinstartpos << " : " << (*Iter).al.GetEndPosition() << endl;
                
                string referenceSequence = reference->getSubSequence(RefIDRedName[(*Iter).al.RefID], (*Iter).al.Position, 2*(*Iter).al.Length);

                if ((*Iter).HasRealign == true)
                {
                    (*Iter).currentReadPosition = 0;
                    (*Iter).currentGenomeSeqPosition = 0;
                    (*Iter).currentAlPosition = (*Iter).al.Position;
                    (*Iter).cigarindex = 0;
                }

                SalRealignInfo talr = (*Iter);
                
                ParseAlignmentsAndExtractVariantsByBaseQualv(AlGroups, talr, Iter, (*Iter).al, referenceSequence, minmaxRefSeqPos, nextwinstartpos, nextwinendpos, (float) minbaseQ, false);

                ++Iter; //Increment iterator
            }
            //****************
        }
    }
    cerr << " DONE: Adjust Alignments CIGARs ..." << endl;
}



