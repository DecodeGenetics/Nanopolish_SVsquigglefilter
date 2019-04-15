//---------------------------------------------------------
// Copyright 2018 deCODE genetics
// Written by Doruk Beyter (dorukb@decode.is)
//---------------------------------------------------------
//
// nanopolish_call_variants -- find variants wrt a reference
//
#ifndef __LIKELIHOODS_CORE_H
#define __LIKELIHOODS_CORE_H



#include "nanopolish_alignment_db.h"
#include "nanopolish_haplotype.h"
#include "nanopolish_variant.h"


#include <string.h>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>



std::vector<SequenceAlignmentRecordInfo> get_min_leftflankrefpos_and_seqrecords(const AlignmentDB& alignments, const std::string& contig, int variant_start, 
	int SV_LEN_CUTOFF_TO_SEPARATE, const std::string& bam_filename, int& min_required_ref_coord);

std::vector<SequenceAlignmentRecordInfo> get_max_rightflankrefpos_and_seqrecords(const AlignmentDB& alignments, const std::string& contig, int variant_start, int variant_ref_len, 
	int SV_LEN_CUTOFF_TO_SEPARATE, const std::string& bam_filename, int& max_required_ref_coord);




HMMInputData get_events_and_read_coords_for_leftflank(const AlignmentDB& alignments, const SequenceAlignmentRecordInfo& rec, int flanking_sequence_amount,
									const int variant_start, int& variant_readPos, int& readpos_start_for_event, int& readpos_end_for_event, bool& successful);

HMMInputData get_events_and_read_coords_for_rightflank(const AlignmentDB& alignments, const SequenceAlignmentRecordInfo& rec, int flanking_sequence_amount,
									const int variant_start, const int variant_ref_len, int& variant_readPos, int& readpos_start_for_event, int& readpos_end_for_event, bool& successful);



Haplotype base_haplotype_for_leftflank(const AlignmentDB& alignments, const SequenceAlignmentRecordInfo& rec, const std::string& contig, const int contig_length,
		const int readpos_start_for_event, const int readpos_end_for_event, const int variant_start, const int variant_readPos, int& ref_beginPos, int& ref_endPos, bool& is_left_flanking);

Haplotype base_haplotype_for_rightflank(const AlignmentDB& alignments, const SequenceAlignmentRecordInfo& rec, const std::string& contig, const int contig_length,
		const int readpos_start_for_event, const int readpos_end_for_event, const int variant_start, const int variant_ref_len, const int variant_readPos, 
		int& ref_beginPos, int& ref_endPos, bool& is_right_flanking);



Haplotype alt_haplotype_for_leftflank_DEL(const AlignmentDB& alignments, const std::string& contig, const int contig_length, Variant test_variant,
		const int readpos_start_for_event, const int readpos_end_for_event, const int variant_start, const int variant_readPos, const int variant_len, const int base_ref_beginPos, bool& is_left_flanking);

Haplotype alt_haplotype_for_rightflank_DEL(const AlignmentDB& alignments, const std::string& contig, const int contig_length, Variant test_variant,
		const int readpos_start_for_event, const int readpos_end_for_event, const int variant_start, const int variant_readPos, const int base_ref_endPos, bool& is_right_flanking);



Haplotype alt_haplotype_for_leftflank_INS(const AlignmentDB& alignments, const std::string& contig, Variant test_variant,
		const int base_ref_beginPos, const int base_ref_endPos);

Haplotype alt_haplotype_for_rightflank_INS(const AlignmentDB& alignments, const std::string& contig, Variant test_variant,
		const int base_ref_beginPos, const int base_ref_endPos);



void left_flank_analysis(const AlignmentDB& alignments, const uint32_t& alignment_flags, const std::string& contig, const int contig_length,
	const std::vector<SequenceAlignmentRecordInfo>& seq_records, 
	const Variant& test_variant, const int flanking_sequence_amount, const std::vector<std::string>& methylation_types, FILE* out_fp);


void right_flank_analysis(const AlignmentDB& alignments, const uint32_t& alignment_flags, const std::string& contig, const int contig_length,
	const std::vector<SequenceAlignmentRecordInfo>& seq_records, 
	const Variant& test_variant, const int flanking_sequence_amount, const std::vector<std::string>& methylation_types, FILE* out_fp);


void print_result_to_csv(FILE* out_fp, const std::string& contig, const SequenceAlignmentRecordInfo& rec, 
	const Variant& test_variant, const int variant_start, const int variant_len, const int event_length, const int flanking_sequence_amount, 
	const double base_score, const double alt_score, const double likelihood_diff, const std::string& flank);


void print_error_to_csv(FILE* out_fp, const std::string& contig, const SequenceAlignmentRecordInfo& rec, 
	const Variant& test_variant, const int variant_start, const int variant_len, const int flanking_sequence_amount, const std::string& flank, const std::string& error_message);

#endif
