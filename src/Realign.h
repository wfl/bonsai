
#ifndef _REALIGN_H
#define	_REALIGN_H

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "Fasta.h"
#include "api/BamAlignment.h"

struct SalRealignInfo
{
    BamTools::BamAlignment al;
    int currentReadPosition;
    int currentGenomeSeqPosition;
    int currentAlPosition;
    int cigarindex;
    bool HasRealign;
    unsigned int FileIdNo;
    std::vector< std::vector<BamTools::CigarOp> > newCigars;
    std::vector<int> bestCigarIndex;
    std::vector<int> longestINDELlength;
    std::vector<int> NumOfIDS;
    std::vector<int> NumOfm;
    std::vector<int> longestContiguousMismatch;
    std::vector<float> mismatches;
    std::vector<float> costF;
    std::vector<int> SWHaplotypePos;
    int CigarSoftclippingLength;    //Original Cigar
};

bool IsWindowInRepeatRegion(std::string sequence);

class RealignFunctionsClass
{
    public:
    struct SubseqInfo
    {
        int32_t                     REFID;
        std::string                 SEQ;
        int                         REFSTARTPOS;
        std::vector<int>            INSLENGTH;
        std::vector<int>            INSREFPOS;
        std::vector<int>            INSREADPOS;
        std::vector<std::string>    INSSEQ;
        std::vector<char>           TYPE;

        // Construct
        SubseqInfo(void) {}

        SubseqInfo(int32_t refid, std::string seq, int refstartpos, std::vector<int> inslength, std::vector<int> insrefpos, std::vector<int> insreadpos, std::vector<std::string> InsSeq, std::vector<char> type)
            :REFID(refid), SEQ(seq), REFSTARTPOS(refstartpos), INSLENGTH(inslength), INSREFPOS(insrefpos), INSREADPOS(insreadpos), INSSEQ(InsSeq), TYPE(type)
        { }
    };


    struct AlleleNPosition
    {
        double POSITION;
        std::vector<std::string> SEQ;

        AlleleNPosition(void) {}
        AlleleNPosition(double pos, std::vector<std::string> seq)
           :POSITION(pos), SEQ(seq)
        { }

        bool operator<(const AlleleNPosition & AP) const
        {
            return (POSITION < AP.POSITION);
        }

        bool operator==(const AlleleNPosition & AP) const
        {
            //return ((POSITION == AP.POSITION) && (SEQ == AP.SEQ));
            return (POSITION == AP.POSITION);
        }

        //bool CompareGreaterThan(AlleleNPosition& vp1, AlleleNPosition& vp2)
        //{
        //    return (vp1.POSITION > vp2.POSITION);
        //}
    };

    struct MismatchNPosition
    {
        double POSITION;
        std::string SEQ;

        MismatchNPosition(void) {}
        MismatchNPosition(double pos, std::string seq)
           :POSITION(pos), SEQ(seq)
        { }

        bool operator<(const MismatchNPosition & MP) const
        {
            return (POSITION < MP.POSITION);
        }

        bool operator==(const MismatchNPosition & MP) const
        {
            return (POSITION == MP.POSITION);
        }
    };


    struct find_Position : std::unary_function<AlleleNPosition, bool>
    {
        double position;
        find_Position(double p):position(p) { }
        bool operator()(AlleleNPosition const& AP) const
        {
            return (AP.POSITION == position);
        }
    };


    struct stringLength
    {
        bool operator() ( const std::string& str1, const std::string& str2 )
        {
            //cerr << str1 << " ? " << str2 << endl;
            return str1.size() < str2.size();
        }
    };


    RealignFunctionsClass(void);
    ~RealignFunctionsClass(void);

    void Clear();
    float getBaseQualityASCII2float(std::string BaseQualstr);
    float getBaseQualityASCII2float(char BaseQualstr);

    float AlleleCostFunction(float allelecounts, float depth, bool isIndel);
    // Establish minimum start and maximum positions of each alignments
    void UpdateStartandEndPositionForHaplotype(BamTools::BamAlignment& alignment, std::vector<int>& minmaxRefSeqPos);
    
    void ParseAlignmentsAndExtractVariantsByBaseQualv(std::vector<SalRealignInfo>& AlGroups, SalRealignInfo& alr, std::vector<SalRealignInfo>::iterator& Iter,
            BamTools::BamAlignment& alignment, const std::string& referenceSequence, std::vector<int>& minmaxRefSeqPos, int winstartpos, int winendpos, float MinBaseQual, bool IsupdateMap);
    bool PruningByNaiveSelectionProcedureAndConstructHaplotypes2(int winstartpos, int winendpos, int32_t RefID, std::string RefIdstring, std::vector<int> minmaxRefSeqPos, FastaReference* reference);
        
    void SelectHaplotypeCandidates_SmithWatermanBSv(std::vector<SalRealignInfo>& AlGroups, std::vector<int> minmaxRefSeqPos, bool SetLowComplexityRegionSWGapExt);

    void AdjustCigarsWRTChosenMultipleHaplotypesAndPrepareAlignmentsTobeWrittenOut(std::vector<SalRealignInfo>& AlGroups, std::multimap<int, BamTools::BamAlignment>& SortRealignedAlignmentsMultimap, FastaReference* reference, std::map<int, std::string>& RefIDRedName, std::vector<int>& minmaxRefSeqPos, int nextwinstartpos, int nextwinendpos, float minbaseQ, int& TotalAlR, int& NumAlR, int ploidy);



private:
        bool m_isLowComplexityRegion;
        std::vector<SubseqInfo> m_Haplotypevec;
        std::vector<int> m_ChosenBestHaplotypesIndex;

        //std::vector<unsigned int> m_CoverageByPositionInWindow;
        std::map<AlleleNPosition, std::vector<int> > m_Position_IndelCounts;   // Key=Variant Position; Value=Type, count, length
        std::map<MismatchNPosition, std::vector<int> > m_Position_MismatchCounts;

        std::vector<MismatchNPosition> m_MismatchesLociList;
        std::vector<AlleleNPosition> m_DeletionLociList;
        std::vector<AlleleNPosition> m_InsertionLociList;
        std::vector<unsigned int> m_CoverageByPositionInWindow;

};

#endif	/* _REALIGN_H */

