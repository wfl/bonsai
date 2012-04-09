
#include <vector>

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <string>
#include <cstdlib>
#include <map>
#include <utility>

#include <stdlib.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <getopt.h> ///* for getopt_long; standard getopt is in unistd.h */

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"

#include "Fasta.h"
#include "GenFunc.h"
#include "Realign.h"

using namespace std;
using namespace BamTools;

// Help command line
void printSummary(char** argv)
{
 cerr << "usage: " << argv[0] << " [options] bam file [... bam file]" << endl
         << endl
         << "options:" << endl
         << "    -f, --reference             FASTA reference against which alignments have been aligned" << endl
         << "    -r, --region <20:1..500 >   Specify region"<< endl
         << "    -R, --Region <region.txt >  Specify a list of regions in a text file"<< endl
         << "    -s, --stdin                 Set stdin flag" << endl
         << "    -b, --bam <merge.bam>       Provide a bam file" << endl
         << "    -q, --min-base-quality <10> Minimum base quality (in Phred quality scores) required for alleles in a read (Default = 10)" << endl
         << "    -w, --window-size <40>      Window size (Default = 40)" << endl
         << "    -p, --ploidy <2>            Currently works for diploid genome only" << endl
         << "    -E, --Repeat-Extgap         Set this flag means align reads(low complexity) in repeat regions with Smith-Waterman Parameter, Gap_Extention = 0 " << endl
         << "    -l, --LowCompex             Set this flag ONLY to evalute and realign the reads located in repeat(low complexity) regions  " << endl
         << "Outputs a BAM file of realigned reads on stdout." << endl
         << endl;

    exit(0);
}

// ***************************************************************
int IsWithinWindow(BamAlignment& alignment, int winstartpos, int winendpos, int AllowableBasesInWindow);
void GetSampleListAndCountInWindow(BamAlignment& alignment, map<string,unsigned int>& samplesAndAlignmentCountsInWindow);
bool ParseRegionString(const string& regionString, const BamReader& reader, BamRegion& region);
bool MeetIndelDetectorThresholdv(const vector<SalRealignInfo>& AlGroups);
// ********************************************************************************


int main (int argc, char** argv)
{

    // Print Commandline
    string ss(argv[0]);   // convert Char to String
    string commandline = "##Print Command line " + ss;

    int c;

    FastaReference* reference = NULL;
    int minbaseQ        = 10;   //default
    int windowlen       = 40;  //by default
    string regionstr;
    string RegionFile;
    string bamfile;
    bool STdin          = false;
    bool has_region     = false;
    bool has_regionFile = false;
    bool has_bamfile    = false;
    bool has_ref        = false;
    int ploidy         = 2;
    bool SetLowComplexityRegionSWGapExt = false;
    bool SetLowComplexityRegion = false;
   

    if (argc < 2)
    {
        printSummary(argv);
        exit(1);
    }

    while (true)
    {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"ploidy", required_argument, 0, 'p'},
            {"window-size", required_argument, 0, 'w'},
            {"reference", required_argument, 0, 'f'},
            {"min-base-quality", required_argument, 0,'q'},
            {"Region", required_argument, 0, 'R'},
            {"STdin", no_argument, 0, 's'},
            {"bam", required_argument, 0, 'b'},
            {"Repeat-Extgap", no_argument, 0, 'E'},
            {"LowCompex", no_argument, 0, 'l'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hslEf:q:w:s:r:R:p:b:", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
            case 'f':
                reference = new FastaReference(optarg); // will exit on open failure
                commandline = commandline + " -f " + optarg;
                has_ref = true;
                break;

            case 'b':
                has_bamfile = true;
                bamfile = optarg;
                commandline = commandline + " -b " + optarg;
                break;

            case 'r':
                regionstr = optarg;
                has_region = true;
                commandline = commandline + " -r " + optarg;
                break;

             case 'R':
                RegionFile = optarg;
                has_regionFile = true;
                commandline = commandline + " -R " + optarg;
                break;

            case 's':
                STdin = true;
                commandline = commandline + " -s ";
                break;
                
            case 'q':
                minbaseQ = atoi(optarg);
                commandline = commandline + " -q " + optarg;
                break;

            case 'w':
                windowlen = atoi(optarg);
                commandline = commandline + " -w " + optarg;
                break;

            case 'p':
                ploidy = atoi(optarg);
                commandline = commandline + " -p " + optarg;
                break;

            case 'E':
                SetLowComplexityRegionSWGapExt = true;
                commandline = commandline + " -E ";
                break;

            case 'l':
                SetLowComplexityRegion = true;
                commandline = commandline + " -l ";
                break;

            case 'h':
                printSummary(argv);
                commandline = commandline + " -h ";
                exit(0);
                break;

            case '?':
                printSummary(argv);
                exit(1);
                break;

              default:
                abort();
                break;
        }
    }

    //// Open Error log files
    ofstream cerrlog("bonsaiReport.txt");
    streambuf *cerrsave = std::cerr.rdbuf();

    // Redirect stream buffers
    if (cerrlog.is_open())
        cerr.rdbuf(cerrlog.rdbuf());

    cerr << commandline << endl;
    

    //Check for Reference Fasta sequence
    if (!has_ref)
    {
        cerr << "no FASTA reference provided, cannot realign" << endl;
        exit(1);
    }

    ////Check for reader
    BamReader reader;
    if (STdin == true)
    {
        if (!reader.Open("stdin"))
        {
            cerr << "could not open stdin bam for reading" << endl;
            cerr << reader.GetErrorString() << endl;
            reader.Close();
            printSummary(argv);
        }
    }
    else
    {
        if (has_bamfile == true)
        {
            if (!reader.Open(bamfile))
            {
                cerr << "ERROR: could not open bam files from stdin ... Aborting" << endl;
                cerr << reader.GetErrorString() << endl;
                reader.Close();
                printSummary(argv);
            }

            if ( !reader.LocateIndex() )
                reader.CreateIndex();
        }
        else
        {
            cerr << "--bam flag is set but no bamfile is provided... Aborting" << endl;
            reader.Close();
            printSummary(argv);
        }
    }

    //// Check Region Tags
    if ( (has_regionFile == true) && (has_region == true) )
    {
        cerr << "ERROR: You provide both region and has provide a Set Region List... Aborting" << endl;
        exit(1);
    }

    //// store the names of all the reference sequences in the BAM file
    vector<RefData> referencedata = reader.GetReferenceData();
   
    //// Store Region LIST
    vector<BamRegion> regionlist;
    if (has_region == true)
    {
        BamRegion region;
        ParseRegionString(regionstr, reader, region);
        regionlist.push_back(region);
    }
    else if (has_regionFile == true)
    {
        ifstream RG(RegionFile.c_str(), ios_base::in);
        string line;
        while(getline(RG,line))
        {
            BamRegion region;
            ParseRegionString(line, reader, region);
            regionlist.push_back(region);
        }
        RG.close();
    }
    else if ( (has_regionFile == false) && (has_region == false) )
    {
        for (int i= 0; i < (int)referencedata.size(); i++)
        {
            string regionstr = referencedata.at(i).RefName;
            BamRegion region;
            ParseRegionString(regionstr, reader, region);
            if (!reader.SetRegion(region)) // Bam region will get [0,101) = 0 to 100 => [closed, half-opened)
            {
                cerr << "ERROR: set region " << regionstr << " failed. Check that REGION describes a valid range... Aborting" << endl;
                reader.Close();
                exit(1);
            }
            else
                regionlist.push_back(region);
        }
    }

    //// 
    BamWriter writer;
    if (!writer.Open("stdout", reader.GetHeaderText(), reader.GetReferenceData()))
    {
        cerr << "could not open stdout for writing" << endl;
        exit(1);
    }

    //// Smallest start position and Largest end position for Req Seq
    vector<RefData>::iterator refdataIter = referencedata.begin();
    vector<BamRegion>::iterator regionListIter = regionlist.begin();
   

    // CLASS
    RealignFunctionsClass RealignFunction;

    map<int, string> RefIDRedName;
    vector<SalRealignInfo> AlGroups;
    multimap<int, BamAlignment> SortRealignedAlignmentsMultimap;

    int refid               = 0;
    BamAlignment alignment;
    bool IsNextAlignment = reader.GetNextAlignment(alignment);
    //cerr << "   " << alignment.Name << " Chr  " << alignment.RefID << " Startpos: " << alignment.Position << " Endpos: " << alignment.GetEndPosition() << " Length: " << alignment.Length << endl;

    int windowrealigned     = 0;
    int TotalWindowDetected = 0;
    int TotalReadsAligned   = 0;
    int TotalWindow         = 0;
    int TotalReads          = 0;

    while (refdataIter != referencedata.end() )
    {
        string refname = refdataIter->RefName;
        RefIDRedName[refid] = refname;
        int reflength = refdataIter->RefLength;
        int winstartpos, winendpos;
        int AllowableBasesInWindow = 1;
        bool nextChrName = false;

        cerr << "##HeaderINFO: RefID = " << refdataIter->RefName << "\t" << "RefLen = " << reflength << endl;
        
        while (nextChrName == false )
        {
            vector<int> minmaxRefSeqPos;
            bool IsPassDetectorNoRealignment = false;
            minmaxRefSeqPos.push_back(-1);
            minmaxRefSeqPos.push_back(0);
            //cerr << " region: " << (*regionListIter).LeftRefID << " : " << (*regionListIter).LeftPosition << " .. " << (*regionListIter).RightPosition << endl;
            if ((refid == (int)referencedata.size() - 1) && ((*regionListIter).LeftRefID == refid) && ((has_region==true) || (has_regionFile==true)) )
            {
                ////
                if ( (has_region == true) || (has_regionFile == true) )
                {                    
                    winstartpos = (*regionListIter).LeftPosition;
                    winendpos   = winstartpos + windowlen - 1;
                    reflength = (*regionListIter).RightPosition;
                    if (reflength < winendpos)
                        reflength = winendpos;
                                       
                    // Get Next Alignment First
                    if ( (refid == alignment.RefID) && (winstartpos == (*regionListIter).LeftPosition) && (IsNextAlignment == false) )
                        IsNextAlignment = reader.GetNextAlignment(alignment);
                }
                else if (has_region == false)
                {
                    winstartpos = 0;
                    winendpos   = winstartpos + windowlen - 1;

                    // Get Next Alignment First
                    if ( (refid == alignment.RefID) && (winstartpos == 0) && (IsNextAlignment == false) )
                        IsNextAlignment = reader.GetNextAlignment(alignment);
                }
                //cerr << " winstart: " << winstartpos << " ; winend: " << winendpos;
                //cerr << "   " << alignment.Name << " Chr  " << alignment.RefID << " Startpos: " << alignment.Position << " Endpos: " << alignment.GetEndPosition() << " Length: " << alignment.Length << endl;

                ////
                while ((winstartpos < reflength))
                {
                    //// Check window end position
                    if (winendpos > reflength)
                        winendpos = reflength;

                    // Reinitialized
                    unsigned int NewReadMappedcount = 0;
                
                    //// Save and Erase alignments that are outside of window (Deque?)
                    if (!AlGroups.empty())
                    {
                        minmaxRefSeqPos.at(0) = -1;
                        minmaxRefSeqPos.at(1) = 0;

                        //cerr << "#Start: Keep alignments with start position exceed the right end of the window/Region " << endl;
                        vector<SalRealignInfo>::iterator Iter = AlGroups.begin();

                        while (Iter != AlGroups.end())
                        {
                            // Erase alignment s
                            if ((*Iter).al.GetEndPosition() < winstartpos)
                            {
                                //cerr << "  ToWrite: " << (*Iter).second.size() << " ; " << (*Iter).al.Name << " ; " << (*Iter).al.Position << " < " << winstartpos << " : " << (*Iter).al.GetEndPosition() << endl;
                                SortRealignedAlignmentsMultimap.insert(pair<int, BamAlignment > ((*Iter).al.Position, (*Iter).al));
                                AlGroups.erase(Iter);

                                //cerr << "  ToWrite: DONE " << endl;
                            } 
                            else
                            {
                                string referenceSequence = reference->getSubSequence(RefIDRedName[(*Iter).al.RefID], (*Iter).al.Position, 2*(*Iter).al.Length);
                            
                                if ((*Iter).HasRealign == true )
                                {
                                    (*Iter).currentReadPosition = 0;
                                    (*Iter).currentGenomeSeqPosition = 0;
                                    (*Iter).currentAlPosition = (*Iter).al.Position;
                                    (*Iter).cigarindex = 0;
                                }

                                (*Iter).CigarSoftclippingLength = 0;
                                SalRealignInfo talr = (*Iter);
                                //cerr << "  ToKEEP: " << (*Iter).al.Name << " ; " << (*Iter).al.Position << " < " << winstartpos << " : " << (*Iter).al.GetEndPosition() << endl;
                                RealignFunction.ParseAlignmentsAndExtractVariantsByBaseQualv(AlGroups, talr, Iter, (*Iter).al, referenceSequence, minmaxRefSeqPos, winstartpos, winendpos, (float) minbaseQ, false);
                           
                                ++Iter; //Increment iterator
                            }
                        }
                    }
                

                    // Write Sorted Alignments that are outside of window
                    //cerr << "SortRealignedAlignmentsMultimap: " << SortRealignedAlignmentsMultimap.size() << " minmaxRefSeqPos.at(0)= " << minmaxRefSeqPos.at(0) << endl;
                    if (!SortRealignedAlignmentsMultimap.empty()) // && (winWrite < winstartpos ) )
                    {
                        //cerr << "#Start: Write alignments and delete alignments with start position exceed the right end of the window/Region " << endl;
                        multimap<int, BamAlignment>::iterator sraIter = SortRealignedAlignmentsMultimap.begin();

                        while (sraIter != SortRealignedAlignmentsMultimap.end()) 
                        {
                            //cerr << " (*sraIter).first= " <<  (*sraIter).first << " minmaxRefSeqPos.at(0)= " << minmaxRefSeqPos.at(0) << " winstartpos - ((windowlen - 1)*0.9)= " << winstartpos - ((windowlen - 1)*0.9) << endl;
                            if (((float) (*sraIter).first < floor((float) (winstartpos - ((windowlen - 1)*0.9)))) && ((minmaxRefSeqPos.at(0) > 0) && ((*sraIter).first < minmaxRefSeqPos.at(0)))) {
                                //writer.SaveAlignment((*sraIter).second);  // Why sometimes, it doesn't work ?????
                                if (!writer.SaveAlignment((*sraIter).second))
                                    cerr << writer.GetErrorString() << endl;

                                SortRealignedAlignmentsMultimap.erase(sraIter++);
                            } else {
                                ++sraIter;
                            }
                    }
                    //cerr << "#Done: Write alignments and delete alignments with start position exceed the right end of the window/Region " << endl;
                    }

                    //cerr << " winstart: " << winstartpos << " ; winend: " << winendpos;
                    //cerr << "   " << alignment.Name << " Chr  " << alignment.RefID << " Startpos: " << alignment.Position << " Endpos: " << alignment.GetEndPosition() << " Length: " << alignment.Length << endl;
                    //cerr <<  ": " << alignment.RefID << " :" << RefIDRedName[alignment.RefID] << " : " << RefIDRedName[alignment.RefID] << endl;

                    //cerr << "Start: Gather aligmenets that lie (fully or partially) within the window frame and group INDELs if there are ... " << endl;
                    // Gather Reads within a window frame
                  
                    while ((IsNextAlignment) && (refid == alignment.RefID)) // Neeed more conditions
                    {
                        if (SetLowComplexityRegion == true) 
                        {
                            string sequenceInWindow = reference->getSubSequence(RefIDRedName[alignment.RefID], winstartpos, (winendpos-winstartpos+1) );

                            if (IsWindowInRepeatRegion(sequenceInWindow) == true)
                            {
                                if ((IsWithinWindow(alignment, winstartpos, winendpos, AllowableBasesInWindow)) == 0)
                                {
                                    TotalReads++;
                                    if (alignment.IsMapped())
                                    {
                                        string referenceSequence = reference->getSubSequence(RefIDRedName[alignment.RefID], alignment.Position, 2*alignment.Length);
 
                                        vector<SalRealignInfo>::iterator tIter;
                                        SalRealignInfo alr;
                                        alr.al = alignment;
                                        alr.currentReadPosition = 0;
                                        alr.currentGenomeSeqPosition = 0;
                                        alr.currentAlPosition = alignment.Position;
                                        alr.cigarindex = 0;
                                        alr.HasRealign = false;
                                        alr.CigarSoftclippingLength = 0;

                                        string str = "ZZZZZZZZZZZZZZZZZ";
                                        if (alignment.Name.find(str) != string::npos) {
                                            stringstream cigar;
                                            for (vector<CigarOp>::const_iterator cigarIter = alignment.CigarData.begin(); cigarIter != alignment.CigarData.end(); ++cigarIter)
                                                cigar << cigarIter->Length << cigarIter->Type;

                                            string cigarstr = cigar.str();
                                            cerr << "   TRACKING: " << alignment.RefID << " " << alignment.Name << " pos: " << alignment.Position << " cigar: " << cigarstr << endl;
                                        }

                                        RealignFunction.ParseAlignmentsAndExtractVariantsByBaseQualv(AlGroups, alr, tIter, alignment, referenceSequence, minmaxRefSeqPos, winstartpos, winendpos, (float) minbaseQ, true);
                                        NewReadMappedcount++;
                                    } 
                                    else
                                    {
                                        SortRealignedAlignmentsMultimap.insert(pair<int, BamAlignment > (alignment.Position, alignment));
                                        cerr << "UNmapped : " << alignment.Name << endl;
                                    }
                                } 
                                else if ((IsWithinWindow(alignment, winstartpos, winendpos, AllowableBasesInWindow)) == 1)
                                {
                                    SortRealignedAlignmentsMultimap.insert(pair<int, BamAlignment > (alignment.Position, alignment));
                                }
                                else
                                    break;
                            } else {
                                if ((IsWithinWindow(alignment, winstartpos, winendpos, AllowableBasesInWindow)) < 2)
                                    SortRealignedAlignmentsMultimap.insert(pair<int, BamAlignment > (alignment.Position, alignment));
                                else
                                    break;
                            }
                        }
                        else // (SetLowComplexityRegion == false)
                        {
                            if ((IsWithinWindow(alignment, winstartpos, winendpos, AllowableBasesInWindow)) == 0)
                            {
                                TotalReads++;
                                if (alignment.IsMapped())
                                {
                                    string referenceSequence = reference->getSubSequence(RefIDRedName[alignment.RefID], alignment.Position, 2 * alignment.Length);

                                    vector<SalRealignInfo>::iterator tIter;
                                    SalRealignInfo alr;
                                    alr.al = alignment;
                                    alr.currentReadPosition = 0;
                                    alr.currentGenomeSeqPosition = 0;
                                    alr.currentAlPosition = alignment.Position;
                                    alr.cigarindex = 0;
                                    alr.HasRealign = false;
                                    alr.CigarSoftclippingLength = 0;

                                    string str = "ZZZZZZZZZZZZZZZZZ";
                                    if (alignment.Name.find(str) != string::npos)
                                    {
                                        stringstream cigar;
                                        for (vector<CigarOp>::const_iterator cigarIter = alignment.CigarData.begin(); cigarIter != alignment.CigarData.end(); ++cigarIter)
                                            cigar << cigarIter->Length << cigarIter->Type;

                                        string cigarstr = cigar.str();
                                        cerr << "   TRACKING: " << alignment.RefID << " " << alignment.Name << " pos: " << alignment.Position << " cigar: " << cigarstr << endl;
                                    }

                                    RealignFunction.ParseAlignmentsAndExtractVariantsByBaseQualv(AlGroups, alr, tIter, alignment, referenceSequence, minmaxRefSeqPos, winstartpos, winendpos, (float) minbaseQ, true);

                                    //cerr << " winstart: " << winstartpos << " ; winend: " << winendpos;
                                    //cerr << "   INDEL: " << alignment.RefID << " " << alignment.Name << " pos: " << alignment.Position << " Length: " << alignment.Length << " CIGARstr: " << cigarstr << endl;
                                    NewReadMappedcount++;
                                } 
                                else
                                {
                                    SortRealignedAlignmentsMultimap.insert(pair<int, BamAlignment > (alignment.Position, alignment));
                                    cerr << "UNmapped : " << alignment.Name << endl;
                                }
                            }
                            else if ((IsWithinWindow(alignment, winstartpos, winendpos, AllowableBasesInWindow)) == 1) {
                                SortRealignedAlignmentsMultimap.insert(pair<int, BamAlignment > (alignment.Position, alignment));
                            }
                            else
                                break;
                        }

                        ////Get next alignment
                        IsNextAlignment = reader.GetNextAlignment(alignment);
                    }

                   //cerr << "Done: Gather aligmenets that lie (fully or partially) within the window frame and group INDELs if there are ... " << endl;

                    //// Detector Corner
                    bool ToRealign = MeetIndelDetectorThresholdv(AlGroups);
                    cerr << "MeetIndelDetectorThresholdv(AlGroups).size()= " << AlGroups.size() << endl;
                    
                    // **************
                    if (ToRealign)
                    {
                        //cerr << "  ToRealign: " << refdataIter->RefName << "\t" << reflength << "\t" << winstartpos << "\t" << winendpos << "\t" << AlGroups.size() << endl;
                        //cerr << "             minmaxRefSeqPos.at(1)= " << minmaxRefSeqPos.at(1) << " minmaxRefSeqPos.at(0)= " << minmaxRefSeqPos.at(0) << endl;

                        ////// Perform Realign routines
                        int TotalAlR = 0; // Total number of alignments to be realigned
                        int NumAlR = 0; // Now many alignments are aligned
                        TotalWindowDetected++;

                        cerr << "#Start: Meet Threshold, Realigning ... " << endl;

                        if (minmaxRefSeqPos.at(1) < winendpos)
                            minmaxRefSeqPos.at(1) = winendpos;

                        if (minmaxRefSeqPos.at(0) > winstartpos)
                            minmaxRefSeqPos.at(0) = winstartpos;

                        bool IsToRealign = RealignFunction.PruningByNaiveSelectionProcedureAndConstructHaplotypes2(winstartpos, winendpos, refid, refname, minmaxRefSeqPos, reference);

                        if (IsToRealign == true)
                        {
                            RealignFunction.SelectHaplotypeCandidates_SmithWatermanBSv(AlGroups, minmaxRefSeqPos, SetLowComplexityRegionSWGapExt);

                            minmaxRefSeqPos.at(0) = -1;
                            minmaxRefSeqPos.at(1) = 0;

                            int nextwinstartpos = winendpos + 1;
                            int nextwinendpos = winstartpos + windowlen - 1;
                            if (nextwinendpos > reflength)
                                nextwinendpos = reflength;

                            //cerr <<  "   Before Realign : " << SortRealignedAlignmentsMultimap.size() << endl;
                            RealignFunction.AdjustCigarsWRTChosenMultipleHaplotypesAndPrepareAlignmentsTobeWrittenOut(AlGroups, SortRealignedAlignmentsMultimap, reference, RefIDRedName, minmaxRefSeqPos, nextwinstartpos, nextwinendpos, minbaseQ, TotalAlR, NumAlR, ploidy);
                            IsPassDetectorNoRealignment = false; // Set flag to false to deactivate write functions

                            //cerr <<  "   After Realign : " << SortRealignedAlignmentsMultimap.size() << endl;

                            TotalReadsAligned += NumAlR;

                            if (NumAlR > 0) // Realignment done
                                windowrealigned++;
                        } else


                        cerr << "#Done: Meet Threshold, Realigning ... " << endl;
                    }


                    if (NewReadMappedcount > 0)
                        TotalWindow++;

                    RealignFunction.Clear();

                    //// Move the window frame
                    winstartpos = winendpos + 1;
                    winendpos = winstartpos + windowlen - 1;
                }

                //// Save and Erase remaining alignments that are outside of window (Deque?)
                if ((!AlGroups.empty())) {
                    cerr << "#Start: Write Remaining alignments and delete all alignments" << endl;

                    for (vector<SalRealignInfo>::iterator Iter = AlGroups.begin(); Iter != AlGroups.end(); ++Iter) {
                        //cerr << "    Remain alignment start: " << (*Iter).al.Name << " " << Iter->al.Position  << " < " << winstartpos << "  " << winendpos << endl;
                        SortRealignedAlignmentsMultimap.insert(pair<int, BamAlignment > ((*Iter).al.Position, (*Iter).al));
                    }

                    cerr << "#Done: Write Remaining alignments and delete all alignments" << endl;
                }

                AlGroups.clear();


                // Write Sorted remaining Alignments that are outside of window
                if (!SortRealignedAlignmentsMultimap.empty())
                {
                    for (multimap<int, BamAlignment>::iterator sraIter = SortRealignedAlignmentsMultimap.begin(); sraIter != SortRealignedAlignmentsMultimap.end(); ++sraIter)
                    {
                        //writer.SaveAlignment((*sraIter).second);
                        if (!writer.SaveAlignment((*sraIter).second))
                            cerr << writer.GetErrorString() << endl;
                    }
                    SortRealignedAlignmentsMultimap.clear();
                }

            }

            ++regionListIter;
            if ((*regionListIter).LeftRefID > refid)
                nextChrName = true;
        }

        //// If End of the chromosome position
        //// increament iterator
        ++refdataIter;
        ++refid;
    }


    reader.Close();
    writer.Close();

    cerr << "##-Completed- " << endl;
    cerr << " Total Reads processed =  " << TotalReads << endl;
    cerr << " Total Reads Aligned =    " << TotalReadsAligned << endl;
    cerr << " Total Window processed = " << TotalWindow << endl;
    cerr << " Total Window Detected =  " << TotalWindowDetected << endl;
    cerr << " Total Windows Aligned =  " << windowrealigned << endl;


    // Restore cerr's stream buffer before terminating
    if (cerrlog.is_open())
        cerr.rdbuf(cerrsave);

    commandline.clear();
    return 0;
}


// Reads are in window

int IsWithinWindow(BamAlignment& alignment, int winstartpos, int winendpos, int AllowableBasesInWindow) {
    //int allowableFracLen = ceil((float) PercentReadLengthInWindow * alignment.Length);

    if ((alignment.Position >= winstartpos) && ((alignment.GetEndPosition() - 1) <= winendpos)) // reads that are exactly in window
        return 0;
    else if (((alignment.Position >= winstartpos) && (alignment.Position <= winendpos))
            || (((alignment.GetEndPosition() - 1) <= winendpos) && ((alignment.GetEndPosition() - AllowableBasesInWindow - 1) >= winstartpos))) // certain reads that are partially in window
        return 0;
    else if ((alignment.Position < winstartpos) && ((alignment.GetEndPosition() - 1) >= winendpos)) // window is in read length
        return 0;
    else if ((alignment.Position < winstartpos) && ((alignment.GetEndPosition() - 1)< winstartpos) )
        return 1;
    else if (alignment.Position > winendpos)
        return 2;
}

void GetSampleListAndCountInWindow(BamAlignment& alignment, map<string, unsigned int>& samplesAndAlignmentCountsInWindow) {

    if (samplesAndAlignmentCountsInWindow.empty())
        samplesAndAlignmentCountsInWindow[alignment.Name.substr(0, 7)] = 1;
    else {
        map<string, unsigned int>::iterator it = samplesAndAlignmentCountsInWindow.find(alignment.Name.substr(0, 7));

        if (it != samplesAndAlignmentCountsInWindow.end())
            it->second++;
        else
            samplesAndAlignmentCountsInWindow[alignment.Name.substr(0, 7)] = 1;
    }
}

bool MeetIndelDetectorThresholdv(const vector<SalRealignInfo>& AlGroups) {
    if ((!AlGroups.empty()) && (AlGroups.size() > 1)) {
        return true;
    } else
        return false;
}

// Same as ParseRegionString() above, but accepts a BamMultiReader
bool ParseRegionString(const string& regionString,
                                  const BamReader& reader,
                                  BamRegion& region)
{
    // -------------------------------
    // parse region string

    // check first for empty string
    if ( regionString.empty() )
        return false;

    //cerr << "ParseRegionString Input: " << regionString << endl;

    // non-empty string, look for a colom
    size_t foundFirstColon = regionString.find(':');

    // store chrom strings, and numeric positions
    string startChrom;
    string stopChrom;
    int startPos;
    int stopPos;

    // no colon found
    // going to use entire contents of requested chromosome
    // just store entire region string as startChrom name
    // use BamReader methods to check if its valid for current BAM file
    if ( foundFirstColon == string::npos ) {
        startChrom = regionString;
        startPos   = 0;
        stopChrom  = regionString;
        stopPos    = -1;
    }

    // colon found, so we at least have some sort of startPos requested
    else {
        // store start chrom from beginning to first colon
        startChrom = regionString.substr(0,foundFirstColon);

        // look for ".." after the colon
        size_t foundRangeDots = regionString.find("..", foundFirstColon+1);

        // no dots found
        // so we have a startPos but no range
        // store contents before colon as startChrom, after as startPos
        if ( foundRangeDots == string::npos )
        {
            startPos   = atoi( regionString.substr(foundFirstColon+1).c_str() );
            stopChrom  = startChrom;
            stopPos    = -1;
        }

        // ".." found, so we have some sort of range selected
        else {

            // store startPos between first colon and range dots ".."
            startPos = atoi( regionString.substr(foundFirstColon+1, foundRangeDots-foundFirstColon-1).c_str() );

            // look for second colon
            size_t foundSecondColon = regionString.find(':', foundRangeDots+1);

            // no second colon found
            // so we have a "standard" chrom:start..stop input format (on single chrom)
            if ( foundSecondColon == string::npos ) {
                stopChrom  = startChrom;
                stopPos    = atoi( regionString.substr(foundRangeDots+2).c_str() );
            }

            // second colon found
            // so we have a range requested across 2 chrom's
            else {
                stopChrom  = regionString.substr(foundRangeDots+2, foundSecondColon-(foundRangeDots+2));
                stopPos    = atoi( regionString.substr(foundSecondColon+1).c_str() );
            }
        }
    }
    
    // -------------------------------
    // validate reference IDs & genomic positions
    const RefVector references = reader.GetReferenceData();

    // if startRefID not found, return false
    int startRefID = reader.GetReferenceID(startChrom);
    if ( startRefID == -1 ) return false;
    // startPos cannot be greater than or equal to reference length
    const RefData& startReference = references.at(startRefID);
    if ( startPos >= startReference.RefLength ) return false;

    // if stopRefID not found, return false
    int stopRefID = reader.GetReferenceID(stopChrom);
    if ( stopRefID == -1 ) return false;

    // stopPosition cannot be larger than reference length
    const RefData& stopReference = references.at(stopRefID);
    if ( stopPos > stopReference.RefLength ) return false;

    // if no stopPosition specified, set to reference end
    if ( stopPos == -1 ) stopPos = stopReference.RefLength;

    // -------------------------------
    // set up Region struct & return

    region.LeftRefID     = startRefID;
    region.LeftPosition  = startPos;
    region.RightRefID    = stopRefID;;
    region.RightPosition = stopPos;

    //cerr << "ParseRegionString " << region.LeftRefID <<  " " << region.LeftPosition << " " << region.RightPosition << endl;
    return true;
}


