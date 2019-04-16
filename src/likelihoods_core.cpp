//---------------------------------------------------------
// Copyright 2018 deCODE genetics
// Written by Doruk Beyter (dorukb@decode.is)
//---------------------------------------------------------
//
// nanopolish_call_variants -- find variants wrt a reference
//

#include "likelihoods_core.h"
#include <climits>
#include <limits>       // std::numeric_limits



std::vector<SequenceAlignmentRecordInfo> get_min_leftflankrefpos_and_seqrecords(const AlignmentDB& alignments, const std::string& contig, int variant_start, 
	int SV_LEN_CUTOFF_TO_SEPARATE, const std::string& bam_filename, int& min_required_ref_coord)
{
	std::vector<SequenceAlignmentRecordInfo> seq_records = alignments.load_region_sequences_info(contig, variant_start-SV_LEN_CUTOFF_TO_SEPARATE/2, variant_start+SV_LEN_CUTOFF_TO_SEPARATE/2, bam_filename);
	min_required_ref_coord = std::numeric_limits<int>::max();
	for (const auto& rec : seq_records)
	{
		if (min_required_ref_coord > rec.beginPos)
		{
			min_required_ref_coord = rec.beginPos;
		}

	}
	// min required ref coord is found. No need for max on left flank because that is hand given.
	return seq_records;
}


std::vector<SequenceAlignmentRecordInfo> get_max_rightflankrefpos_and_seqrecords(const AlignmentDB& alignments, const std::string& contig, int variant_start, int variant_ref_len, 
	int SV_LEN_CUTOFF_TO_SEPARATE, const std::string& bam_filename, int& max_required_ref_coord)
{
	int variant_right_flank = variant_start + variant_ref_len;
	std::vector<SequenceAlignmentRecordInfo> seq_records = alignments.load_region_sequences_info(contig, variant_right_flank-SV_LEN_CUTOFF_TO_SEPARATE/2, variant_right_flank+SV_LEN_CUTOFF_TO_SEPARATE/2, bam_filename);
	max_required_ref_coord = 0;
	for (const auto& rec : seq_records)
	{
		if (max_required_ref_coord < rec.endPos)
		{
			max_required_ref_coord = rec.endPos;
		}

	}
	// max required ref coord is found. No need for max on left flank because that is hand given.
	return seq_records;
}


 

HMMInputData get_events_and_read_coords_for_leftflank(const AlignmentDB& alignments, const SequenceAlignmentRecordInfo& rec, int flanking_sequence_amount,
									const int variant_start, int& variant_readPos, int& readpos_start_for_event, int& readpos_end_for_event, bool& success)
{
	int read_beginPos = 0;
	alignments._find_read_pos_from_ref_pos(rec.aligned_bases, rec.beginPos, read_beginPos);
	//std::cout << "read_beginPos: " << read_beginPos << std::endl;

	alignments._find_read_pos_from_ref_pos(rec.aligned_bases, variant_start, variant_readPos);
	//std::cout << variant_start << " variant_readPos: " << variant_readPos << std::endl;

	int variant_leftflank_inread = 0;
	alignments._find_read_pos_from_ref_pos(rec.aligned_bases, (variant_start - flanking_sequence_amount), variant_leftflank_inread);
	//std::cout  << " variant_leftflank_inread: " << variant_leftflank_inread << std::endl;
	//std::cout  << " just the read's variant pos - flanking: " <<  (variant_readPos-flanking_sequence_amount) << std::endl;

	//TODO This max should actually be removed.
	readpos_start_for_event = std::max(read_beginPos, variant_leftflank_inread);
	//std::cout << "readpos_start_for_event: " << readpos_start_for_event << std::endl;

	// now the end.
	int read_endPos = rec.sequence_len -1;
	//std::cout << "read_endPos: " << read_endPos << std::endl;

	// int readpos_end_for_event = std::min(read_endPos, (variant_readPos + flanking_sequence_amount));
	readpos_end_for_event = std::min(read_endPos, (readpos_start_for_event + 2*flanking_sequence_amount));
	//std::cout << "readpos_end_for_event: " << readpos_end_for_event << std::endl;


	// Make sure somehow that begin end read positions are meaningful and covers the deletion or insertion.
	//TODO:

	// now get the event begin and end.
	int event_begin_idx, event_end_idx;
	success = alignments.find_event_coords_for_given_read_coords(rec, readpos_start_for_event, readpos_end_for_event, event_begin_idx, event_end_idx);

	if (success)
	{
		//std::cout << "event begin and end: " << event_begin_idx << " - " << event_end_idx << std::endl;
		HMMInputData event_sequence_scrappie = alignments.get_given_event_subsequences_for_record_for_events(event_begin_idx, event_end_idx, rec);
						
		return event_sequence_scrappie;
	}
	else
	{
		HMMInputData event_sequence_scrappie;
		return event_sequence_scrappie;
	}
	
}


HMMInputData get_events_and_read_coords_for_rightflank(const AlignmentDB& alignments, const SequenceAlignmentRecordInfo& rec, int flanking_sequence_amount,
									const int variant_start, const int variant_ref_len, int& variant_readPos, int& readpos_start_for_event, int& readpos_end_for_event, bool& success)
{

	// first fix the end site.
	int read_endPos;
	alignments._find_read_pos_from_ref_pos(rec.aligned_bases, rec.endPos, read_endPos);
	//std::cout << "read_endPos: " << read_endPos << std::endl;

	// get the variant pos
	alignments._find_read_pos_from_ref_pos(rec.aligned_bases, variant_start + variant_ref_len, variant_readPos);
	//std::cout << " variant_readPos: " << variant_readPos << std::endl;

	// if it was spanning the variant, this would be the right flank place.
	int variant_rightflank_inread = 0;
	alignments._find_read_pos_from_ref_pos(rec.aligned_bases, (variant_start  + variant_ref_len + flanking_sequence_amount), variant_rightflank_inread);
	//std::cout  << " variant_rightflank_inread: " << variant_rightflank_inread << std::endl;

	// pick the min.
	readpos_end_for_event = std::min(read_endPos, variant_rightflank_inread);
	//std::cout << "readpos_end_for_event: " << readpos_end_for_event << std::endl;


	// now the begin site.
	int read_beginPos = 0;
	//std::cout << "read_beginPos: " << read_beginPos << std::endl;

	readpos_start_for_event = std::max(read_beginPos, (readpos_end_for_event - 2*flanking_sequence_amount));
	//std::cout << "readpos_start_for_event: " << readpos_start_for_event << std::endl;

	// now get the event begin and end.
	int event_begin_idx, event_end_idx;
	success = alignments.find_event_coords_for_given_read_coords(rec, readpos_start_for_event, readpos_end_for_event, event_begin_idx, event_end_idx);

	if (success)
	{
		//std::cout << "event begin and end: " << event_begin_idx << " - " << event_end_idx << std::endl;
		HMMInputData event_sequence_scrappie = alignments.get_given_event_subsequences_for_record_for_events(event_begin_idx, event_end_idx, rec);
		return event_sequence_scrappie;
	}
	else
	{
		HMMInputData event_sequence_scrappie;
		return event_sequence_scrappie;
	}
}






Haplotype base_haplotype_for_leftflank(const AlignmentDB& alignments, const SequenceAlignmentRecordInfo& rec, const std::string& contig, const int contig_length,
		const int readpos_start_for_event, const int readpos_end_for_event, const int variant_start, const int variant_readPos, int& ref_beginPos, int& ref_endPos, bool& is_left_flanking)
{

	//int ref_beginPos;
	alignments._find_ref_pos_from_read_pos(rec.aligned_bases, readpos_start_for_event, ref_beginPos);

	is_left_flanking = true;
	if (variant_start <= ref_beginPos)
	{
		//std::cout << " Assumption violated, the breakpoint is more off than the flank site in LEFT FLANK." << std::endl;
		is_left_flanking = false;
		return Haplotype();
	}

	ref_endPos = ref_beginPos + (readpos_end_for_event - readpos_start_for_event);
	//std::cout << "ref_endPos: " << ref_endPos << std::endl;

	if ((variant_start >= ref_endPos) || (ref_endPos >= contig_length))
	{
		//std::cout << " Assumption violated, the breakpoint is more off than the flank site in LEFT FLANK." << std::endl;
		is_left_flanking = false;
		return Haplotype();
	}

	Haplotype base_haplotype(contig, ref_beginPos, alignments.get_reference_substring(contig, ref_beginPos, ref_endPos));
	
	return base_haplotype;
}


Haplotype base_haplotype_for_rightflank(const AlignmentDB& alignments, const SequenceAlignmentRecordInfo& rec, const std::string& contig, const int contig_length,
		const int readpos_start_for_event, const int readpos_end_for_event, const int variant_start, const int variant_ref_len, const int variant_readPos, int& ref_beginPos, int& ref_endPos, bool& is_right_flanking)
{

	alignments._find_ref_pos_from_read_pos(rec.aligned_bases, readpos_end_for_event, ref_endPos);
	//std::cout << "ref_endPos: " << ref_endPos << std::endl;


	is_right_flanking = true;
	if ( (variant_start + variant_ref_len)  >= ref_endPos)
	{
		//std::cout << " Assumption violated, the breakpoint is more off than the flank site in RIGHT FLANK." << std::endl;
		is_right_flanking = false;
		return Haplotype();
	}


	//ref_beginPos = variant_start + variant_ref_len - (variant_readPos - readpos_start_for_event);
	ref_beginPos = ref_endPos - (readpos_end_for_event - readpos_start_for_event);
	//std::cout << "ref_beginPos: " << ref_beginPos << std::endl;

	if (( (variant_start + variant_ref_len)  <= ref_beginPos) || (ref_beginPos < 0))
	{
		//std::cout << " Assumption violated, the breakpoint is more off than the flank site in RIGHT FLANK." << std::endl;
		is_right_flanking = false;
		return Haplotype();
	}


	Haplotype base_haplotype(contig, ref_beginPos, alignments.get_reference_substring(contig, ref_beginPos, ref_endPos));
	//std::cout << "base_haplotype: " << base_haplotype.get_sequence() << std::endl;

	return base_haplotype;
}





Haplotype alt_haplotype_for_leftflank_DEL(const AlignmentDB& alignments, const std::string& contig, const int contig_length, Variant test_variant,
		const int readpos_start_for_event, const int readpos_end_for_event, const int variant_start, const int variant_readPos, const int variant_len, const int base_ref_beginPos, bool& is_left_flanking)
{
	//std::cout << "variant_len: " << variant_len << std::endl;
	//int alt_ref_endPos = variant_start + variant_len - 1 + (readpos_end_for_event - variant_readPos);
	int alt_ref_endPos = base_ref_beginPos + test_variant.ref_length - 1 + (readpos_end_for_event - readpos_start_for_event);

	if (alt_ref_endPos >= contig_length)
	{
		//std::cout << "the alternate reference end position exceeds the contig length." << std::endl;
		is_left_flanking = false;
		return Haplotype();
	}

	//std::cout << "base_ref_beginPos: " << base_ref_beginPos << std::endl;
	//std::cout << "alt_ref_endPos: " << alt_ref_endPos << std::endl;

	Haplotype alt_haplotype(contig, base_ref_beginPos, alignments.get_reference_substring(contig, base_ref_beginPos, alt_ref_endPos));
	alt_haplotype.apply_variant(test_variant);
	//std::cout << "alt_haplotype: " << alt_haplotype.get_sequence() << std::endl;
	is_left_flanking = true;
	return alt_haplotype;
}

Haplotype alt_haplotype_for_rightflank_DEL(const AlignmentDB& alignments, const std::string& contig, const int contig_length, Variant test_variant,
		const int readpos_start_for_event, const int readpos_end_for_event, const int variant_start, const int variant_readPos, const int base_ref_endPos, bool& is_right_flanking)
{

	//std::cout << "variant_len: " << test_variant.ref_length << std::endl;
	// int alt_ref_beginPos = variant_start - (variant_readPos - readpos_start_for_event);
	int alt_ref_beginPos = base_ref_endPos - test_variant.ref_length + 1 - (readpos_end_for_event - readpos_start_for_event);

	if (alt_ref_beginPos < 0)
	{
		is_right_flanking = false;
		return Haplotype();
	}

	//std::cout << "base_ref_endPos: " << base_ref_endPos << std::endl;
	//std::cout << "alt_ref_beginPos: " << alt_ref_beginPos << std::endl;

	Haplotype alt_haplotype(contig, alt_ref_beginPos, alignments.get_reference_substring(contig, alt_ref_beginPos, base_ref_endPos));
	alt_haplotype.apply_variant(test_variant);
	//std::cout << "alt_haplotype: " << alt_haplotype.get_sequence() << std::endl;
	is_right_flanking = true;
	return alt_haplotype;

}


Haplotype alt_haplotype_for_leftflank_INS(const AlignmentDB& alignments, const std::string& contig, Variant test_variant,
		const int base_ref_beginPos, const int base_ref_endPos)
{
	
	Haplotype alt_haplotype(contig, base_ref_beginPos, alignments.get_reference_substring(contig, base_ref_beginPos, base_ref_endPos));
	alt_haplotype.apply_variant(test_variant);
	alt_haplotype.truncate_seq_from_right_end(test_variant.alt_length-1);
	//std::cout << "alt_haplotype: " << alt_haplotype.get_sequence() << std::endl;

	return alt_haplotype;

}


Haplotype alt_haplotype_for_rightflank_INS(const AlignmentDB& alignments, const std::string& contig, Variant test_variant,
		const int base_ref_beginPos, const int base_ref_endPos)
{

	//std::cout << "variant_len: " << variant_len << std::endl;
	Haplotype alt_haplotype(contig, base_ref_beginPos, alignments.get_reference_substring(contig, base_ref_beginPos, base_ref_endPos));
	alt_haplotype.apply_variant(test_variant);
	alt_haplotype.apply_variant(test_variant);
	alt_haplotype.truncate_seq_from_left_end(test_variant.alt_length-1);
	//std::cout << "alt_haplotype: " << alt_haplotype.get_sequence() << std::endl;

	return alt_haplotype;

}



void left_flank_analysis(const AlignmentDB& alignments, const uint32_t& alignment_flags, const std::string& contig, const int contig_length,
	const std::vector<SequenceAlignmentRecordInfo>& seq_records, 
	const Variant& test_variant, const int flanking_sequence_amount, const std::vector<std::string>& methylation_types, FILE* out_fp)
{

	int variant_start = test_variant.ref_position;
	int variant_len = std::max(test_variant.ref_length, test_variant.alt_length);
	for (const auto& rec : seq_records)
	{
		try
		{
			// Make sure you don't get alignments that softclipped there.
			int read_beginPos = 0;
			bool found = alignments._find_read_pos_from_ref_pos(rec.aligned_bases, rec.beginPos, read_beginPos);
			if (!found)
			{
				print_error_to_csv(out_fp, contig, rec, test_variant, variant_start, variant_len, flanking_sequence_amount, "LEFT", "BEGINPOS_ERROR");
				continue;

			}
			
			int variant_leftflank_inread = 0;
			found = alignments._find_read_pos_from_ref_pos(rec.aligned_bases, (variant_start - flanking_sequence_amount), variant_leftflank_inread);
			if (!found)
			{
				print_error_to_csv(out_fp, contig, rec, test_variant, variant_start, variant_len, flanking_sequence_amount, "LEFT", "LEFTFLANKPOS_NOTFOUND");
				continue;
			}

			//std::cout  << "If those are equal then the alignment is left soft clipped and we don't want it." << std::endl;
			//std::cout << "read_beginPos: " << read_beginPos << std::endl;
			//std::cout << "variant_leftflank_inread: " << variant_leftflank_inread << std::endl;
			if (std::binary_search(test_variant.supporting_reads.begin(), test_variant.supporting_reads.end(), rec.read_name) && read_beginPos != variant_leftflank_inread)
			{
				std::cout << rec.read_name << std::endl;
				//std::cout << "Read supports the SV!" << std::endl;
			
				// get the event sequence.
				int variant_readPos, readpos_start_for_event, readpos_end_for_event;
				bool successful;
				HMMInputData event_sequence_scrappie = get_events_and_read_coords_for_leftflank(alignments, rec, flanking_sequence_amount,
																				variant_start, variant_readPos, readpos_start_for_event, readpos_end_for_event, successful);

				if (!successful)
				{	
					print_error_to_csv(out_fp, contig, rec, test_variant, variant_start, variant_len, flanking_sequence_amount, "LEFT", "BAD_QUALITY");
					continue;
				}

				int event_length =  event_sequence_scrappie.event_stop_idx - event_sequence_scrappie.event_start_idx + 1;
				// get the base haplotype for left flank
				int base_ref_beginPos, base_ref_endPos;
				bool is_left_flanking;
				Haplotype base_haplotype = base_haplotype_for_leftflank(alignments, rec, contig, contig_length, readpos_start_for_event, readpos_end_for_event, 
					variant_start, variant_readPos, base_ref_beginPos, base_ref_endPos, is_left_flanking);

				if (!is_left_flanking)
				{
					print_error_to_csv(out_fp, contig, rec, test_variant, variant_start, variant_len, flanking_sequence_amount, "LEFT", "POS_TOO_AWAY");
					continue;
				}


				// get the ref likelihood.
				int base_seqlen, base_hmm_input_len;
				double base_score = score_haplotype_for_record(base_haplotype, event_sequence_scrappie, alignment_flags, methylation_types, base_seqlen, base_hmm_input_len);
				
				std::cout << "base_score: " << base_score << std::endl;
				std::cout << "base_seqlen: " << base_seqlen << std::endl;
				

				// get the alternative haplotype for left flank
				Haplotype alt_haplotype;
				if (test_variant.alt_seq == "<DEL>")
				{
					alt_haplotype = alt_haplotype_for_leftflank_DEL(alignments, contig, contig_length, test_variant, readpos_start_for_event, readpos_end_for_event, variant_start, variant_readPos, variant_len, base_ref_beginPos, is_left_flanking);
					if (!is_left_flanking)
					{
						print_error_to_csv(out_fp, contig, rec, test_variant, variant_start, variant_len, flanking_sequence_amount, "LEFT", "POS_TOO_AWAY");
						continue;
					}

				}
				else if (test_variant.alt_seq == "<INS>")
				{
					alt_haplotype = alt_haplotype_for_leftflank_INS(alignments, contig, test_variant, base_ref_beginPos, base_ref_endPos);				
				}
				

				// get the alt likelihood.
				int alt_seqlen, alt_hmm_input_len;
				double alt_score = score_haplotype_for_record(alt_haplotype, event_sequence_scrappie, alignment_flags, methylation_types, alt_seqlen, alt_hmm_input_len);
				std::cout << "alt_score: " << alt_score << std::endl;
				std::cout << "alt_seqlen: " << alt_seqlen << std::endl;	

				double likelihood_diff = alt_score - base_score;

				print_result_to_csv(out_fp, contig, rec, test_variant, variant_start, variant_len, event_length, flanking_sequence_amount, base_score, alt_score, likelihood_diff, "LEFT");
				

			}
			else
			{
				std::cout << "Do not compute likelihood for now, because read does not support the SV or is softclipping on the left flank." << std::endl;
			}

		}
		catch (std::runtime_error const & e)
		{
			print_error_to_csv(out_fp, contig, rec, test_variant, variant_start, variant_len, flanking_sequence_amount, "LEFT", "NAN_SCORE");
			continue;

		}

	}
	return;
}



void right_flank_analysis(const AlignmentDB& alignments, const uint32_t& alignment_flags, const std::string& contig, const int contig_length,
	const std::vector<SequenceAlignmentRecordInfo>& seq_records, 
	const Variant& test_variant, const int flanking_sequence_amount, const std::vector<std::string>& methylation_types, FILE* out_fp)
{

	int variant_start = test_variant.ref_position;
	int variant_len = std::max(test_variant.ref_length, test_variant.alt_length);
	for (const auto& rec : seq_records)
	{
		try
		{
			int read_beginPos = 0;
			bool found = alignments._find_read_pos_from_ref_pos(rec.aligned_bases, rec.beginPos, read_beginPos);
			
			if (!found)
			{	print_error_to_csv(out_fp, contig, rec, test_variant, variant_start, variant_len, flanking_sequence_amount, "RIGHT", "BEGINPOS_ERROR");
				continue;
			}
		
			int variant_rightflank_inread = 0;
			bool right_flank_is_valid = alignments._find_read_pos_from_ref_pos(rec.aligned_bases, (variant_start + test_variant.ref_length + flanking_sequence_amount), variant_rightflank_inread);

			if (std::binary_search(test_variant.supporting_reads.begin(), test_variant.supporting_reads.end(), rec.read_name) && right_flank_is_valid)
			{
			
				// get the event sequence.
				int variant_readPos, readpos_start_for_event, readpos_end_for_event;
				bool successful;
				HMMInputData event_sequence_scrappie = get_events_and_read_coords_for_rightflank(alignments, rec, flanking_sequence_amount,
																				variant_start, test_variant.ref_length, variant_readPos, 
																				readpos_start_for_event, readpos_end_for_event, successful);

				if (!successful)
				{
					print_error_to_csv(out_fp, contig, rec, test_variant, variant_start, variant_len, flanking_sequence_amount, "RIGHT", "BAD_QUALITY");
					continue;
				}

				int event_length =  event_sequence_scrappie.event_stop_idx - event_sequence_scrappie.event_start_idx + 1;
				// get the base haplotype for left flank
				int base_ref_beginPos, base_ref_endPos;
				bool is_right_flanking;
				Haplotype base_haplotype = base_haplotype_for_rightflank(alignments, rec, contig, contig_length, readpos_start_for_event, readpos_end_for_event, 
					variant_start, test_variant.ref_length, variant_readPos, base_ref_beginPos, base_ref_endPos, is_right_flanking);

				if (!is_right_flanking)
				{
					print_error_to_csv(out_fp, contig, rec, test_variant, variant_start, variant_len, flanking_sequence_amount, "RIGHT", "POS_TOO_AWAY");
					continue;
				}

				// get the ref likelihood.
				int base_seqlen, base_hmm_input_len;
				double base_score = score_haplotype_for_record(base_haplotype, event_sequence_scrappie, alignment_flags, methylation_types, base_seqlen, base_hmm_input_len);
				std::cout << "base_score: " << base_score << std::endl;
				std::cout << "base_seqlen: " << base_seqlen << std::endl;
				
				Haplotype alt_haplotype;
				if (test_variant.alt_seq == "<DEL>")
				{
					alt_haplotype = alt_haplotype_for_rightflank_DEL(alignments, contig, contig_length, test_variant, readpos_start_for_event, 
						readpos_end_for_event, variant_start, variant_readPos, base_ref_endPos, is_right_flanking);
					if (!is_right_flanking)
					{
						print_error_to_csv(out_fp, contig, rec, test_variant, variant_start, variant_len, flanking_sequence_amount, "RIGHT", "POS_TOO_AWAY");
						continue;
					}
				}
				else if (test_variant.alt_seq == "<INS>")
				{
					alt_haplotype = alt_haplotype_for_rightflank_INS(alignments, contig, test_variant, base_ref_beginPos, base_ref_endPos);		
				}


				// get the alt likelihood.
				int alt_seqlen, alt_hmm_input_len;
				double alt_score = score_haplotype_for_record(alt_haplotype, event_sequence_scrappie, alignment_flags, methylation_types, alt_seqlen, alt_hmm_input_len);
				std::cout << "alt_score: " << alt_score << std::endl;
				std::cout << "alt_seqlen: " << alt_seqlen << std::endl;	

				double likelihood_diff = alt_score - base_score;
				
				print_result_to_csv(out_fp, contig, rec, test_variant, variant_start, variant_len, event_length, flanking_sequence_amount, base_score, alt_score, likelihood_diff, "RIGHT");
				

			}
			else
			{
				std::cout << rec.read_name << std::endl;
				std::cout << "It's read_beginPos was: " << read_beginPos << std::endl;
				std::cout << "Do not compute likelihood for now, because read does not support the SV or is softclipping on the left flank." << std::endl;
			}

		}
		catch (std::runtime_error const & e)
		{
			print_error_to_csv(out_fp, contig, rec, test_variant, variant_start, variant_len, flanking_sequence_amount, "RIGHT", "NAN_SCORE");
			continue;
		}

	}
	
}


void print_result_to_csv(FILE* out_fp, const std::string& contig, const SequenceAlignmentRecordInfo& rec, 
	const Variant& test_variant, const int variant_start, const int variant_len, const int event_length, const int flanking_sequence_amount, 
	const double base_score, const double alt_score, const double likelihood_diff, const std::string& flank)
{
	int breakpoint = variant_start;
	if (flank == "RIGHT")
	{
		breakpoint = variant_start + test_variant.ref_length;
	}

	fprintf(out_fp, "%s\t%s\t%s\t%d\t%d\t%s\t%d\t", contig.c_str(), test_variant.SV_ID.c_str(), flank.c_str(), variant_start, breakpoint, test_variant.alt_seq.c_str(), variant_len);
	fprintf(out_fp, "%s\t%d\t%d\t%d\t%d\t%f\t%f\t%f\n", rec.read_name.c_str(), rec.sequence_len, event_length, rec.map_quality, flanking_sequence_amount, base_score, alt_score, likelihood_diff);


}

void print_error_to_csv(FILE* out_fp, const std::string& contig, const SequenceAlignmentRecordInfo& rec, 
	const Variant& test_variant, const int variant_start, const int variant_len, const int flanking_sequence_amount, const std::string& flank, const std::string& error_message)
{

	int breakpoint = variant_start;
	if (flank == "RIGHT")
	{
		breakpoint = variant_start + test_variant.ref_length;
	}

	fprintf(out_fp, "%s\t%s\t%s\t%d\t%d\t%s\t%d\t", contig.c_str(), test_variant.SV_ID.c_str(), flank.c_str(), variant_start, breakpoint, test_variant.alt_seq.c_str(), variant_len);
	fprintf(out_fp, "%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n", rec.read_name.c_str(), rec.sequence_len, -1, rec.map_quality, flanking_sequence_amount, "XXX", "XXX", error_message.c_str());
}