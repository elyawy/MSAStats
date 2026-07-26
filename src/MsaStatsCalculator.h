#ifndef _MSASTATS_CALCULATOR_H
#define _MSASTATS_CALCULATOR_H
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <map>
#include "read_seqs.h"
#include "summary_stats_defs.h"
using namespace std;

class MsaStatsCalculator
{
public:

	// MsaStatsCalculator(const MSA &msa); 

	// This construction also computes all the summary statistics.
	MsaStatsCalculator(const vector<string> & seqArray); 
	
	// This construction reads a FASTA file into the MSA and computes summary statistics
	MsaStatsCalculator(string filename);

	// This construction accepts a run-length encoded (sparse) MSA representation directly,
	// e.g. from a simulator's internal MSA object, avoiding building full aligned strings.
	// sparseRuns[row] is a list of run lengths for that row: positive = non-gap run,
	// negative = gap run (abs value = run length). Same convention as Sailfish's SparseMSA.
	MsaStatsCalculator(const vector<vector<int>> & sparseRuns, int numberOfSequences, int msaLength);

	// MSA(){};

	int getMSAlength() const {return _isSparseInput ? _msaLengthSparse : (int)_originalAlignedSeqs[0].size();}
	int getNumberOfSequences() const {return _numberOfSequences;} 
	int getTotalNumberOfIndels() const {return _totalNumberOfIndels;}
	int getTotalNumberOfUniqueIndels() const {return _totalNumberOfUniqueIndels;}
	int getNumberOfIndelsOfLengthOne() const {return _numberOfIndelsOfLengthOne;}
	int getNumberOfIndelsOfLengthTwo() const {return _numberOfIndelsOfLengthTwo;}
	int getNumberOfIndelsOfLengthThree() const {return _numberOfIndelsOfLengthThree;}
	int getNumberOfIndelsOfLengthAtLeastFour() const {return _numberOfIndelsOfLengthAtLeastFour;}
	double getAverageIndelSize() const {return _aveIndelLength;}
	double getAverageUniqueIndelSize() const {return _aveUniqueIndelLength;}
	int getMSALongestSeqLength()const {return _longestSeqLength;}
	int getMSAshortestSeqLength()const {return _shortestSeqLength;}
	int getnumberOfIndelsOfLengthOneInOnePosition()const {return _numberOfIndelsOfLengthOneInOnePosition;}
	int getnumberOfIndelsOfLengthOneInTwoPositions()const {return _numberOfIndelsOfLengthOneInTwoPositions;}
	int getnumberOfIndelsOfLengthOneInNMinus1Positions()const {return _numberOfIndelsOfLengthOneInNMinus1Positions;}
	int getnumberOfIndelsOfLengthTwoInOnePosition()const {return _numberOfIndelsOfLengthTwoInOnePosition;}
	int getnumberOfIndelsOfLengthTwoInTwoPositions()const {return _numberOfIndelsOfLengthTwoInTwoPositions;}
	int getnumberOfIndelsOfLengthTwoInNMinus1Positions()const {return _numberOfIndelsOfLengthTwoInNMinus1Positions;}
	int getnumberOfIndelsOfLengthThreeInOnePosition()const {return _numberOfIndelsOfLengthThreeInOnePosition;}
	int getnumberOfIndelsOfLengthThreeInTwoPositions()const {return _numberOfIndelsOfLengthThreeInTwoPositions;}
	int getnumberOfIndelsOfLengthThreeInNMinus1Positions()const {return _numberOfIndelsOfLengthThreeInNMinus1Positions;}
	int getnumberOfIndelsOfLengthAtLeastFourInOnePosition()const {return _numberOfIndelsOfLengthAtLeastFourInOnePosition;}
	int getnumberOfIndelsOfLengthAtLeastFourInTwoPositions()const {return _numberOfIndelsOfLengthAtLeastFourInTwoPositions;}
	int getnumberOfIndelsOfLengthAtLeastFourInNMinus1Positions()const {return _numberOfIndelsOfLengthAtLeastFourInNMinus1Positions;}
	size_t getNumberOfMSA_position_with_0_gaps() const {return _numberOfMSA_position_with_0_gaps;}
	size_t getNumberOfMSA_position_with_1_gaps() const { return _numberOfMSA_position_with_1_gaps;}
	size_t getNumberOfMSA_position_with_2_gaps() const { return _numberOfMSA_position_with_2_gaps;}
	size_t getNumberOfMSA_position_with_n_minus_1_gaps() const { return _numberOfMSA_position_with_n_minus_1_gaps;}

	double getStatValByType(stat_type statTypeToGet);
	vector<string> getUnalignedSeqs() const;
	vector<string> getAlignedSeqs () const {return _originalAlignedSeqs;}
	void printMSA();
    void recomputeStats();

	vector<double> getStatVec();

	
	~MsaStatsCalculator();

private:
	const vector<string> _originalAlignedSeqs; //The aligned sequences

	vector<string> _alignedSeqs; //The aligned sequences
	int _numberOfSequences; // NUMBER OF SEQUENCES IN THE MSA

	// --- sparse (run-length encoded) input path ---
	bool _isSparseInput = false;
	int _msaLengthSparse = 0;
	vector<vector<int>> _originalSparseSeqs; //run-length rows as given (no trimming needed: no all-gap columns can occur)
	double _aveIndelLength;
	int _totalNumberOfIndels;
	int _longestSeqLength;
	int _shortestSeqLength;
	//vector<string> _unalignedSeqs; //The unaligned sequences

	int _numberOfIndelsOfLengthOne; //counts the number of indels of size 1 over all sequences
	int _numberOfIndelsOfLengthTwo; //counts the number of indels of size 2 over all sequences
	int _numberOfIndelsOfLengthThree; //counts the number of indels of size 3 over all sequences
	int _numberOfIndelsOfLengthAtLeastFour; //counts the number of indels of size 4+ over all sequences
	
	int _numberOfIndelsOfLengthOneInOnePosition;
	int _numberOfIndelsOfLengthOneInTwoPositions;
	int _numberOfIndelsOfLengthOneInNMinus1Positions;
	int _numberOfIndelsOfLengthTwoInOnePosition;
	int _numberOfIndelsOfLengthTwoInTwoPositions;
	int	_numberOfIndelsOfLengthTwoInNMinus1Positions;
	int _numberOfIndelsOfLengthThreeInOnePosition;
	int _numberOfIndelsOfLengthThreeInTwoPositions;
	int	_numberOfIndelsOfLengthThreeInNMinus1Positions;
	int _numberOfIndelsOfLengthAtLeastFourInOnePosition;
	int _numberOfIndelsOfLengthAtLeastFourInTwoPositions;
	int	_numberOfIndelsOfLengthAtLeastFourInNMinus1Positions;

	size_t _numberOfMSA_position_with_0_gaps;
	size_t _numberOfMSA_position_with_1_gaps;
	size_t _numberOfMSA_position_with_2_gaps;
	size_t _numberOfMSA_position_with_n_minus_1_gaps;
	

	map<pair<int,int>,vector<int>> _uniqueIndelMap;
	vector<int> _indelCounter;
	// unique indels summary statistics
	double _aveUniqueIndelLength;
	int _totalNumberOfUniqueIndels;

	//methods
	void setValuesOfIndelSummStats();
	void fillUniqueGapsMap();
	//void unalignSeqs();
	void setLongestAndShortestSequenceLengths();
	void trimMSAFromAllIndelPositionAndgetSummaryStatisticsFromIndelCounter();

	void initializeAllVariables();

	// sparse-path equivalents
	void setValuesOfIndelSummStatsSparse();
	void buildIndelMapAndPositionStatsSparse();
	void setLongestAndShortestSequenceLengthsSparse();

};
#endif