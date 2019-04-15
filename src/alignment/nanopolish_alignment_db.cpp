//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_alignment_db -- abstraction for working
// with sets of reads/events aligned to a reference genome
//
#include <assert.h>
#include <algorithm>
#include "nanopolish_alignment_db.h"
#include "htslib/faidx.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "nanopolish_methyltrain.h"


#include <iostream>
#include <string>


// Various file handle and structures
// needed to traverse a bam file
struct BamHandles
{
	htsFile* bam_fh;
	bam1_t* bam_record;
	hts_itr_t* itr;
	bam_hdr_t* bam_hdr;
};


char complement(char n)
{   
	switch(n)
	{   
	case 'A':
		return 'T';
	case 'T':
		return 'A';
	case 'G':
		return 'C';
	case 'C':
		return 'G';
	}   
	assert(false);
	return ' ';
}  


bool is_homopolymer(const std::string& s)
{
	return ( (s.find_first_not_of(s[0]) == std::string::npos) && (s[0] != 'N') ) ;
}

void fill_indices_to_missed_bases(std::vector<int>& base_event_indices, int kmer_size)
{

	// for(int k = 0; k < base_event_indices.size(); k++)
	// {
	// 	std::cout << base_event_indices[k]<< std::endl;
	// }


	// first, fill the last k-1 base indices with the latest that is filled.
	if (base_event_indices[base_event_indices.size() - kmer_size] != -1)
	{
		for (int i = 1; i < kmer_size; i++)
		{
			base_event_indices[base_event_indices.size() - kmer_size + i] = base_event_indices[base_event_indices.size() - kmer_size];
		}

	}
	else
	{
		std::cout << " Last substr must be filled. Aborting." <<std::endl;
		std::abort();
	}

	//For any other -1 in the middle, take the average and assign the prev one (floor of avg)
	int last_seen_base = 0;
	int first_tobeseen_base = 0;
	for (int i = 0; i < base_event_indices.size(); i++)
	{
		if (base_event_indices[i] >= 0)
		{
			last_seen_base = base_event_indices[i];
			first_tobeseen_base = base_event_indices[i];
		}
		else
		{
			int d = i;
			while (d < base_event_indices.size() && base_event_indices[d] < 0)
			{
				d++;
				first_tobeseen_base = base_event_indices[d];
			}

			// now take the avg of the two.
			base_event_indices[i] = (last_seen_base + first_tobeseen_base) / 2;
			last_seen_base = base_event_indices[i];
		}
	}


	// CHECK IF EVENTS ARE MONOTONICALLY INCREASING.
	int cur_event = -1;
	for (int k = 0; k < base_event_indices.size(); k++)
	{
		int next_event = base_event_indices[k];
		if (next_event < cur_event)
		{
			std::cout << " The event indices here must be monotonically increasing or remaining the same." << std::endl;
			std::abort();
		}
		else
		{
			cur_event = next_event;
		}

	}

	return;

}


//functions
void get_begin_end_event_indices_for_read_region(const SequenceAlignmentRecordInfo& seq_record, const std::vector<int>& event_inds_for_bases, 
                                    const int k_mer_length, const int read_start, const int read_end, int& event_begin_idx, int& event_end_idx)
{
    if (seq_record.rc)
    {
    	
        //std::cout << "rc: " << (bool)seq_record.rc << std::endl;

        event_begin_idx = event_inds_for_bases[seq_record.sequence.length() -1 - read_start - (k_mer_length-1)];
        event_end_idx = event_inds_for_bases[seq_record.sequence.length() -1 - read_end];
    }
    else
    {
        //std::cout << "nay" << std::endl;
        event_begin_idx = event_inds_for_bases[read_start];
        event_end_idx = event_inds_for_bases[read_end - (k_mer_length-1)];
    }

    return;
}



BamHandles _initialize_bam_itr(const std::string& bam_filename,
							   const std::string& contig,
							   int start_position,
							   int stop_position)
{
	BamHandles handles;

	// load bam file
	handles.bam_fh = sam_open(bam_filename.c_str(), "r");
	assert(handles.bam_fh != NULL);

	// load bam index file
	hts_idx_t* bam_idx = sam_index_load(handles.bam_fh, bam_filename.c_str());
	if(bam_idx == NULL) {
		bam_index_error_exit(bam_filename);
	}

	// read the bam header to get the contig ID
	bam_hdr_t* hdr = sam_hdr_read(handles.bam_fh);
	handles.bam_hdr = hdr;

	int contig_id = bam_name2id(hdr, contig.c_str());


	
	// Initialize iteration
	handles.bam_record = bam_init1();
	handles.itr = sam_itr_queryi(bam_idx, contig_id, start_position, stop_position);

	hts_idx_destroy(bam_idx);
	//bam_hdr_destroy(hdr);
	return handles;
}



SequenceAlignmentRecordInfo::SequenceAlignmentRecordInfo(const bam1_t* record, const bam_hdr_t *bamHdr)
{
	this->read_name = bam_get_qname(record);
	this->beginPos = record->core.pos;
	this->endPos = bam_endpos(record) - 1;
	//char * chrom = m_bamHdr->target_name[record->core.tid];
	this->chromosome = bamHdr->target_name[record->core.tid];
	//this->rc = bam_is_rev(record);
	//uint32_t len = aln->core.l_qseq; //length of the read.

	this->rc = bam_is_rev(record);
	

	// copy sequence out of the record
	uint8_t* pseq = bam_get_seq(record);
	this->sequence.resize(record->core.l_qseq);
	for(int i = 0; i < record->core.l_qseq; ++i) {
		this->sequence[i] = seq_nt16_str[bam_seqi(pseq, i)];
	}
	this->sequence_len = this->sequence.length();
	// copy read base-to-reference alignment
	std::vector<AlignedSegment> alignments = get_aligned_segments(record);
	if(alignments.size() > 1) {
		fprintf(stderr, "Errorr: spliced alignments detected when loading read %s\n", this->read_name.c_str());
		fprintf(stderr, "Please align the reads to the genome using a non-spliced aligner\n");
		exit(EXIT_FAILURE);
	}
	assert(!alignments.empty());
	this->aligned_bases = alignments[0];
	this->map_quality = record->core.qual;

	//this->base_to_event_map = ;

}


//added by dorukb
bool AlignmentDB::find_event_coords_for_region_for_record(const SequenceAlignmentRecordInfo& sequence_record, 
										const int start_pos, const int stop_pos, int& event_begin_idx, int& event_end_idx)
{

	int r1 = 0, r2 = 0;
	bool bounded = _find_by_ref_bounds(sequence_record.aligned_bases, start_pos, stop_pos, r1, r2); //r1 and r2 (read start/end updated.)

	if (!bounded)
	{
		std::cout << "Must be bounded. Aborting." << std::endl;
		std::abort();
	}

	// std::cout << "read start pos: " << r1 << std::endl;
	// std::cout << "read stop pos: " << r2 << std::endl;
	// std::cout << "sequence: " << sequence_record.sequence.substr(r1,r2-r1+1) << std::endl;
	std::vector<int> event_inds_for_bases;
	bool success = find_scrappie_events_for_basecall(sequence_record, event_inds_for_bases);

	if (success)
	{
		event_begin_idx = 0;
		event_end_idx = 0;
		int k_mer_length = 5;
		get_begin_end_event_indices_for_read_region(sequence_record, event_inds_for_bases, k_mer_length, r1, r2, event_begin_idx, event_end_idx);
		return success;
	}
	else
	{
		return success;
	}
}



//added by dorukb
bool AlignmentDB::find_event_coords_for_given_read_coords(const SequenceAlignmentRecordInfo& sequence_record, 
										const int r1, const int r2, int& event_begin_idx, int& event_end_idx) const
{
	
	//std::cout << "sequence: " << sequence_record.sequence.substr(r1,r2-r1+1) << std::endl;
	//std::vector<int> event_inds_for_bases = find_scrappie_events_for_basecall(sequence_record);
	std::vector<int> event_inds_for_bases;
	bool success = find_scrappie_events_for_basecall(sequence_record, event_inds_for_bases);

	if (success)
	{
		event_begin_idx = 0;
		event_end_idx = 0;
		int k_mer_length = 5;
		get_begin_end_event_indices_for_read_region(sequence_record, event_inds_for_bases, k_mer_length, r1, r2, event_begin_idx, event_end_idx);
		return success;
	}
	else
	{
		return success;

	}

}


// added by dorukb
bool AlignmentDB::find_event_coords_for_given_read_coords_nanopolish(const SequenceAlignmentRecordInfo& seq_align_info, 
										const int r1, const int r2, int& event_begin_idx, int& event_end_idx) const
{

	SquiggleRead* sr;
	if (m_squiggle_read_map.find(seq_align_info.read_name) != m_squiggle_read_map.end())
	{
		//std::cout << "read is found when preparing HMMInput data." << std::endl;
		sr = m_squiggle_read_map.at(seq_align_info.read_name);
	}
	else
	{
		std::cout << "read is NOT found when calling for preparing HMMinput data." << std::endl;
		std::abort();
	}



	const EventAlignmentRecord& record = EventAlignmentRecord(sr, 0, seq_align_info);
	if(record.aligned_events.empty()) {
		std::cout << "record's aligned events is empty." << std::endl;
		//continue;
		std::abort();
	}

	if(!record.sr->has_events_for_strand(record.strand)) {
		std::cout << "record doesn't have events for the strand." << std::endl;
		//continue;
		std::abort();
	}


	//int e1,e2;
	//std::cout << "record.aligned_events.size(): " << record.aligned_events.size() << std::endl;

	AlignedPairRefLBComp lb_comp;
	auto start_iter = std::lower_bound(record.aligned_readpos_events.begin(), record.aligned_readpos_events.end(),
								  r1, lb_comp);

	auto stop_iter = std::lower_bound(record.aligned_readpos_events.begin(), record.aligned_readpos_events.end(),
								 r2, lb_comp);


	if(start_iter == record.aligned_readpos_events.end() || stop_iter == record.aligned_readpos_events.end())
	{
		std::cout << "iter ends are not bounded." << std::endl;
		return false;
	}




	event_begin_idx = start_iter->read_pos;
	event_end_idx = stop_iter->read_pos;
	// std::cout << "NP event_begin_idx: " <<  event_begin_idx << std::endl;
	// std::cout << "NP event_end_idx: " <<  event_end_idx << std::endl;
	// bool bounded = _find_by_ref_bounds(record.aligned_readpos_events, 
	// 								   r1, 
	// 								   r2,
	// 								   event_begin_idx,
	// 								   event_end_idx);

	return true;


}




bool AlignmentDB::find_scrappie_events_for_basecall(const SequenceAlignmentRecordInfo& sequence_record, std::vector<int>& event_indices_for_bases) const
{


	SquiggleRead* sr;
	if (m_squiggle_read_map.find(sequence_record.read_name) != m_squiggle_read_map.end())
	{
		//std::cout << "read is found in the squiggle." << std::endl;
		sr = m_squiggle_read_map.at(sequence_record.read_name);
	}
	else
	{
		std::cout << "read is NOT found in the squiggle." << std::endl;
		std::abort();
	}



	std::string basecall = sr->read_sequence;
	// std::cout <<  basecall.substr(0, 20)  <<  " IS THE FIRST PREFIX." << std::endl;
	// if (sequence_record.rc)
	// {
	// 	transform(begin(basecall), end(basecall), begin(basecall), complement);
	// 	reverse(basecall.begin(), basecall.end());
	// }
	// std::cout << basecall.substr(0, 20) <<  " IS THE AFTER PREFIX." << std::endl;
	// std::cout << basecall.length() <<  " IS THE len of the basecall." << std::endl;





	//std::vector<int> event_indices_for_bases(basecall.length(), -1);
	event_indices_for_bases.resize(basecall.length(), -1);

	int kmer_size = 5;
	//SquiggleRead* sr = new SquiggleRead(sequence_record.read_name, m_read_db);

	
	
	//std::cout << "NUM_STRANDS: " << (int) NUM_STRANDS << std::endl;

	bool success = false;
	for(size_t si = 0; si < NUM_STRANDS; ++si) 
	{
		//std::cout << "si: " << si << std::endl;
 		// skip complement
 		//std::cout << "C_IDX: " << (int) C_IDX << std::endl;
		//std::cout <<  "sr->has_events_for_strand(si)" << (bool) sr->has_events_for_strand(si) << std::endl;
		if( (si != C_IDX) && (sr->has_events_for_strand(si)) )
		{
			//std::vector<SquigglePosAndState> event_pos_and_states = sr->event_pos_and_states;

			//std::cout << "basecalled event size: " << sr->basecalled_scrappie_events[si].size() << std::endl;
			auto et = sr->basecalled_scrappie_events[si];
			int kmer_size = et[0].kmer.length();

			unsigned int pos_missed_offset = 0;
			unsigned int string_idx = 0;
			int i = 0;
			success = true;
			//std::cout << "et.size: " << et.size() << std::endl;
			while (i < et.size())
			{
				if (et[i].pos == 0)
				{
					if (et[i].kmer == basecall.substr(0, kmer_size))
					{
						event_indices_for_bases[0] = et[i].idx;
					}
					else
					{
						std::cout << "first event must match " << std::endl;
						std::cout << et[i].kmer << " is the et[i].kmer" <<  std::endl;
						std::cout << basecall.substr(0, kmer_size) << " is the basecall.substr(0, kmer_size)" <<  std::endl;
						std::cout << sequence_record.read_name << " is the read name."  <<  std::endl;
						std::cout << sequence_record.chromosome << " is the chromosome."  <<  std::endl;
						std::cout << sequence_record.beginPos << " is the begin pos"  <<  std::endl;
						std::cout <<  basecall.substr(0, 100)  <<  " IS THE seq." << std::endl;
						std::cout << sequence_record.sequence_len << " is the seq len"  <<  std::endl;
						std::abort();
					}
				}
				//for the non-initial events non-stay events,
				else
				{
					// if matches
					// assign event to here
					assert((et[i].pos + pos_missed_offset + kmer_size) <= basecall.length());

					std::string cur_substr = basecall.substr(et[i].pos + pos_missed_offset, kmer_size);
					if ( et[i].kmer == cur_substr )
					{
						event_indices_for_bases[et[i].pos + pos_missed_offset] = et[i].idx;
						
						// if this is the last homopolymer in the stretch (next one is not)

						// find the next non NNNNN event.
						// find the next kmer.
						std::string next_event = "";
						int next_event_idx = i + 1;
						int next_substr_idx =  et[i].pos + pos_missed_offset + 1;
						
						//if both done at the same time, event alignment is done and the world is a pretty and peaceful place.

						if ( (next_substr_idx + kmer_size) > basecall.length() &&  next_event_idx >= et.size() )
						{
							//std::cout << "event alignment peacefully completed." << std::endl;

							break;
						}
						else if ( !( (next_substr_idx + kmer_size) <= basecall.length() &&  next_event_idx < et.size()) )
						{

							//If the reason for ending is that the sequence ended
							if (!( (next_substr_idx + kmer_size) <= basecall.length()))
							{
								// if current substr and last event kmer are homopolymers and matching.
								if (is_homopolymer(cur_substr) && 
									is_homopolymer(et[et.size()-1].kmer) &&
									cur_substr == et[et.size()-1].kmer )
								{
									event_indices_for_bases[et[i].pos + pos_missed_offset] = et[et.size()-1].idx;
									//std::cout << "event alignment successfully completed. as the sequence ended." << std::endl;
									break;
								}
							}
							else if (!(next_event_idx < et.size()))
							{
								//this is the last event, means the sequence has still homopolymers to go, check.
								// and assign all remaining substr to last event idx.
								int advance = 1;
								//std::cout << "event indices are done. but the sequence is not over." << std::endl;
								//std::cout << "et[i].pos + advance + pos_missed_offset + kmer_size: " << et[i].pos + advance + pos_missed_offset + kmer_size << std::endl;
								//std::cout << "basecall.length(): " << basecall.length() << std::endl;
								while ( (et[i].pos + advance + pos_missed_offset + kmer_size) <= basecall.length() )
								{
									//std::cout << "et[i].pos + advance + pos_missed_offset: " << et[i].pos + advance + pos_missed_offset << std::endl;
									//std::cout << "event_indices_for_bases.size(): " << event_indices_for_bases.size() << std::endl;

									//std::cout << "event_indices_for_bases.size(): " << event_indices_for_bases << std::endl;


									event_indices_for_bases[et[i].pos + advance + pos_missed_offset] = et[i].idx;
									advance++;
								}
								break;


							}
							else
							{
								std::cout << " i: " <<  i << std::endl;
								std::cout << " next_substr_idx: " <<  next_substr_idx << std::endl;
								std::cout << "(next_substr_idx + kmer_size): " << (next_substr_idx + kmer_size) << std::endl;
								std::cout << " basecall.length(): " <<  basecall.length() << std::endl;
								std::cout << "next_event_idx: " <<  next_event_idx << std::endl;
								std::cout << "et.size(): " <<  et.size() << std::endl;
								std::cout << "partially done. Error" << std::endl;
								std::abort();

							}

						}
						// else, keep on.
						// std::cout << "next_event_idx: " << next_event_idx << std::endl;
						// std::cout << "et.size(): " << et.size() << std::endl;

						// std::cout << "next_substr_idx: " << next_substr_idx << std::endl;
						// std::cout << "basecall.length(): " << basecall.length() << std::endl;


						next_event = et[next_event_idx].kmer;
						std::string next_substr = basecall.substr(next_substr_idx, kmer_size);
						
						// std::cout << "Base: " << et[i].pos + pos_missed_offset << " matches event: " << i << std::endl;
						// std::cout << "POS: " << et[i].pos  << " and offset: " << pos_missed_offset  << std::endl;
						// std::cout << "et[i].kmer: " << et[i].kmer << std::endl;
						// std::cout << "next_event: " << next_event << std::endl;
						// std::cout << "next_substr: " << next_substr << std::endl;
						// std::cout << "curr_substr: " << basecall.substr(et[i].pos + pos_missed_offset, kmer_size) << std::endl;


						//sequence has the last homopolymer in the strech, but next event is still homopolymer.
						if ( is_homopolymer(cur_substr) && !is_homopolymer(next_substr))
						{
							//advance in the events until you match the last homopolymer in the events.
							while ( (i+1) < et.size() &&  is_homopolymer(et[i+1].kmer) )
							{
								i++;
								pos_missed_offset--;
								//std::cout << "in while; pos, offset: " << et[i].pos << ",  " << pos_missed_offset <<std::endl;
							}
							if (i >= et.size() )
							{
								//std::cout << " not all substrings are aligned. Abort." << std::endl;
								std::abort();
							}
							//std::cout << "out while, pos, offset: " << et[i].pos << ",  " << pos_missed_offset <<std::endl;
							event_indices_for_bases[et[i].pos + pos_missed_offset] = et[i].idx;
							
						}




						// event is last homopolymer in the strecth
						// make sure to advance in the sequnce and assign the its index as much as you can.
						else if ( is_homopolymer(et[i].kmer) && ( !is_homopolymer(next_event)) )
						{
							//std::cout << "IN LOOP mer: " << et[i].kmer << " and next kmer: " << next_event << std::endl;
							// assign this event index to the homopolymer stretch in the basecalls until it ends.
							while ( ((et[i].pos + pos_missed_offset + 1 + kmer_size) <= basecall.length() )
								&& is_homopolymer(basecall.substr(et[i].pos + pos_missed_offset + 1, kmer_size)) )
							{
								pos_missed_offset++;
								event_indices_for_bases[et[i].pos + pos_missed_offset] = et[i].idx;
								//std::cout << "et[i].pos + pos_missed_offset: " << et[i].pos + pos_missed_offset << std::endl;
								//std::cout << "event_indices_for_bases.size(): " << event_indices_for_bases.size() << std::endl;
							}
							if ((et[i].pos + pos_missed_offset + 1 + kmer_size) > basecall.length() )
							{
								// std::cout << "index substr: " << (et[i].pos + pos_missed_offset + 1 + kmer_size) << std::endl;
								// std::cout << " not all events are aligned. Abort." << std::endl;
								std::abort();
							}
							
						}
					}
					else
					{
						
						std::cout << "i: " << i  << " kmer: " << et[i].kmer << " pos: "<< et[i].pos 
							<< " pos_missed_offset: " << pos_missed_offset 
							<< " substr: " << basecall.substr(et[i].pos + pos_missed_offset, kmer_size) << std::endl;
						std::cout << "No mismatch of events and substrings should occur. Aborting." << std::endl;
						std::abort();
					}
				}
				
			i++;
			}
			
			break;
		}

	 }

	 if (success)
	 {
	 	fill_indices_to_missed_bases(event_indices_for_bases, kmer_size);
	 	//return event_indices_for_bases;
	 	return success;
	 }

	 return success;

	 


}


//added by dorukb
std::vector<SequenceAlignmentRecordInfo> AlignmentDB::load_region_sequences_info(const std::string& contig, 
												int start_position, int stop_position, 
												const std::string& sequence_bam) const
{
	assert(!contig.empty());
	assert(start_position >= 0);
	assert(stop_position >= 0);

	BamHandles handles = _initialize_bam_itr(sequence_bam, contig, start_position, stop_position);
	//bam_hdr_t* hdr = sam_hdr_read(handles.bam_fh);



	std::vector<SequenceAlignmentRecordInfo> records;

	int result;
	while((result = sam_itr_next(handles.bam_fh, handles.itr, handles.bam_record)) >= 0) {

		// skip ambiguously mapped reads
		//std::cout << "load region sequences INFO now" << std::endl;
		if(handles.bam_record->core.qual < 1) {
		    continue;
		}

		records.emplace_back(handles.bam_record, handles.bam_hdr);

	}

	// cleanup
	sam_itr_destroy(handles.itr);
	bam_destroy1(handles.bam_record);
	sam_close(handles.bam_fh);
	return records;
}

//added by dorukb
std::vector<SequenceAlignmentRecordInfo> AlignmentDB::load_all_sequences_info(const std::string& sequence_bam)
{
	samFile *fp_in = hts_open(sequence_bam.c_str(),"r"); //open bam file
	bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header
	bam1_t *aln = bam_init1(); //initialize an alignment

	std::vector<SequenceAlignmentRecordInfo> records;
	while(sam_read1(fp_in, bamHdr, aln) > 0)
	{
		if(aln->core.qual < 1) {
		    continue;
		}
		//int32_t pos = aln->core.pos; //left most position of alignment in zero based coordianate (+1)
		//char *chr = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
		//uint32_t len = aln->core.l_qseq; //length of the read.
		records.emplace_back(aln, bamHdr);
	}

	// cleanup
	//sam_itr_destroy(handles.itr);
	bam_destroy1(aln);
	sam_close(fp_in);
	
	return records;
}











//
// SequenceAlignmentRecord
//
SequenceAlignmentRecord::SequenceAlignmentRecord(const bam1_t* record)
{
	this->read_name = bam_get_qname(record);
	this->rc = bam_is_rev(record);
	

	// copy sequence out of the record
	uint8_t* pseq = bam_get_seq(record);
	this->sequence.resize(record->core.l_qseq);
	for(int i = 0; i < record->core.l_qseq; ++i) {
		this->sequence[i] = seq_nt16_str[bam_seqi(pseq, i)];
	}
	
	// copy read base-to-reference alignment
	std::vector<AlignedSegment> alignments = get_aligned_segments(record);
	if(alignments.size() > 1) {
		fprintf(stderr, "Error: spliced alignments detected when loading read %s\n", this->read_name.c_str());
		fprintf(stderr, "Please align the reads to the genome using a non-spliced aligner\n");
		exit(EXIT_FAILURE);
	}
	assert(!alignments.empty());
	this->aligned_bases = alignments[0];
}

//
// EventAlignmentRecord
//
EventAlignmentRecord::EventAlignmentRecord(SquiggleRead* sr,
										   const int strand_idx,
										   const SequenceAlignmentRecord& seq_record)
{ 
	this->sr = sr;
	size_t k = this->sr->get_model_k(strand_idx);
	size_t read_length = this->sr->read_sequence.length();
	
	for(size_t i = 0; i < seq_record.aligned_bases.size(); ++i) {
		// skip positions at the boundary
		if(seq_record.aligned_bases[i].read_pos < k) {
			continue;
		}

		if(seq_record.aligned_bases[i].read_pos + k >= read_length) {
			continue;
		}

		size_t kmer_pos_ref_strand = seq_record.aligned_bases[i].read_pos;
		size_t kmer_pos_read_strand = seq_record.rc ? this->sr->flip_k_strand(kmer_pos_ref_strand, k) : kmer_pos_ref_strand;
		size_t event_idx = this->sr->get_closest_event_to(kmer_pos_read_strand, strand_idx);
		this->aligned_events.push_back( { seq_record.aligned_bases[i].ref_pos, (int)event_idx });
		this->aligned_readpos_events.push_back( { seq_record.aligned_bases[i].read_pos, (int)event_idx });
		// #if DEBUG_PRINT_STATS
		// std::cout << "kmer_pos_ref_strand: " << kmer_pos_ref_strand << std::endl;
		// std::cout << "kmer_pos_read_strand: " << kmer_pos_read_strand << std::endl; 
		// std::cout << "event_idx: " << event_idx << std::endl; 
		// #endif
	}
	this->rc = strand_idx == 0 ? seq_record.rc : !seq_record.rc;
	this->strand = strand_idx;
	this->stride = this->aligned_events.front().read_pos < this->aligned_events.back().read_pos ? 1 : -1;
	
	this->read_name = seq_record.read_name;
}


EventAlignmentRecord::EventAlignmentRecord(SquiggleRead* sr,
										   const int strand_idx,
										   const SequenceAlignmentRecordInfo& seq_record)
{
	this->sr = sr;
	size_t k = this->sr->get_model_k(strand_idx);
	size_t read_length = this->sr->read_sequence.length();
	
	for(size_t i = 0; i < seq_record.aligned_bases.size(); ++i) {
		// skip positions at the boundary
		if(seq_record.aligned_bases[i].read_pos < k) {
			continue;
		}

		if(seq_record.aligned_bases[i].read_pos + k >= read_length) {
			continue;
		}

		size_t kmer_pos_ref_strand = seq_record.aligned_bases[i].read_pos;
		size_t kmer_pos_read_strand = seq_record.rc ? this->sr->flip_k_strand(kmer_pos_ref_strand, k) : kmer_pos_ref_strand;
		size_t event_idx = this->sr->get_closest_event_to(kmer_pos_read_strand, strand_idx);
		this->aligned_events.push_back( { seq_record.aligned_bases[i].ref_pos, (int)event_idx });
		this->aligned_readpos_events.push_back( { seq_record.aligned_bases[i].read_pos, (int)event_idx });
		// #if DEBUG_PRINT_STATS
		// std::cout << "kmer_pos_ref_strand: " << kmer_pos_ref_strand << std::endl;
		// std::cout << "kmer_pos_read_strand: " << kmer_pos_read_strand << std::endl; 
		// std::cout << "event_idx: " << event_idx << std::endl; 
		// #endif


	}
	this->rc = strand_idx == 0 ? seq_record.rc : !seq_record.rc;
	this->strand = strand_idx;
	this->stride = this->aligned_events.front().read_pos < this->aligned_events.back().read_pos ? 1 : -1;
	//std::cout << "STRIDE IS: " << this->stride << std::endl;
	this->read_name = seq_record.read_name;
}



//
// AlignmentDB
//

AlignmentDB::AlignmentDB(const std::string& reads_file,
						 const std::string& reference_file,
						 const std::string& sequence_bam,
						 const std::string& event_bam) :
							m_reference_file(reference_file),
							m_sequence_bam(sequence_bam),
							m_event_bam(event_bam)
{   
	// added by dorukb
	// samFile *fp_in = hts_open(sequence_bam.c_str(),"r"); //open bam file
	// m_bamHdr = sam_hdr_read(fp_in); //read header
	// added by dorukb

	m_read_db.load(reads_file);
	_clear_region();
}

AlignmentDB::~AlignmentDB()
{
	_clear_region();
}

void AlignmentDB::set_alternative_basecalls_bam(const std::string& alternative_basecalls_bam)
{
	m_alternative_basecalls_bam = alternative_basecalls_bam;
}

bool AlignmentDB::are_coordinates_valid(const std::string& contig,
										int start_position,
										int stop_position) const
{
	return m_region_contig == contig &&
		   start_position >= m_region_start &&
		   stop_position <= m_region_end;
}

std::string AlignmentDB::get_reference_substring(const std::string& contig,
												 int start_position,
												 int stop_position) const
{
	if(!are_coordinates_valid(contig, start_position, stop_position)) {
		fprintf(stderr, "[alignmentdb] error: requested coordinates "
				"[%s %d %d] is outside of region boundary [%s %d %d]\n",
				contig.c_str(), start_position, stop_position,
				m_region_contig.c_str(), m_region_start, m_region_end);
		exit(EXIT_FAILURE);
	}

	return m_region_ref_sequence.substr(start_position - m_region_start, stop_position - start_position + 1);
}


std::vector<std::string> AlignmentDB::get_read_substrings(const std::string& contig,
														  int start_position,
														  int stop_position) const
{
	assert(m_region_contig == contig);
	assert(m_region_start <= start_position);
	assert(m_region_end >= stop_position);

	std::vector<std::string> out;
	for(size_t i = 0; i < m_sequence_records.size(); ++i) {
		const SequenceAlignmentRecord& record = m_sequence_records[i];
		if(record.aligned_bases.empty())
			continue;

		int r1, r2;
		bool bounded = _find_by_ref_bounds(record.aligned_bases,
										   start_position,
										   stop_position,
										   r1,
										   r2);

		if(bounded) {
			out.push_back(record.sequence.substr(r1, r2 - r1 + 1));
		}
	}
	return out;
}


//added by dorukb
HMMInputData AlignmentDB::get_event_subsequence_for_record(const std::string& contig,
															  int& start_position,
															  int& stop_position, 
															  const SequenceAlignmentRecordInfo& seq_align_info) const
{
	
	bool hmm_data_ready = false;
	HMMInputData data;


	assert(m_region_contig == contig); 
	assert(m_region_start <= start_position);
	assert(m_region_end >= stop_position);
	// std::cout << "m_region_contig: " << m_region_contig << std::endl;
	// std::cout << "contig: " << contig << std::endl;

	// std::cout << "m_region_start: " << m_region_start << std::endl;
	// std::cout << "start_position: " << start_position << std::endl;

	// std::cout << "m_region_end: " << m_region_end << std::endl;
	// std::cout << "stop_position: " << stop_position << std::endl;

	// std::cout << "m_event_records.size() = " << m_event_records.size() << std::endl;
	for(size_t i = 0; i < m_event_records.size(); ++i) {
		const EventAlignmentRecord& record = m_event_records[i];
		if (record.read_name != seq_align_info.read_name)
		{
			continue;
		}

		if(record.aligned_events.empty()) {
			std::cout << "record's aligned events is empty." << std::endl;
			continue;
		}

		if(!record.sr->has_events_for_strand(record.strand)) {
			std::cout << "record doesn't have events for the strand." << std::endl;
			continue;
		}

		if (start_position < record.aligned_events[0].ref_pos)
		{
			start_position = record.aligned_events[0].ref_pos;
		}
		if (stop_position > record.aligned_events[record.aligned_events.size()-1].ref_pos)
		{
			stop_position = record.aligned_events[record.aligned_events.size()-1].ref_pos;
		}
		
	

		data.read = record.sr;
		data.pore_model = record.sr->get_base_model(record.strand);
		data.strand = record.strand;
		data.rc = record.rc;
		data.event_stride = record.stride;
		
		int e1,e2;
		//std::cout << "record.aligned_events.size(): " << record.aligned_events.size() << std::endl;
		bool bounded = _find_by_ref_bounds(record.aligned_events, 
										   start_position, 
										   stop_position,
										   e1,
										   e2);


		//std::cout << "bounded is: " << bounded << std::endl;

		if(bounded) {
			double ratio = fabs(e1 - e2) / fabs(stop_position - start_position);
			 //std::cout << "ratio is: " << ratio << std::endl;
			 //std::cout << "MAX_EVENT_TO_BP_RATIO is: " << MAX_EVENT_TO_BP_RATIO << std::endl;
			// Some low quality reads appear to have "stuck" states where you get 
			// a long run of consecutive stays. They can cause an assertion in the HMM
			// so we added this heuristic to catch these.
			if(ratio < MAX_EVENT_TO_BP_RATIO) {
				assert(e1 >= 0);
				assert(e2 >= 0);
				//std::cout << "e1 is: " << e1 << std::endl;
				//std::cout << "e2 is: " << e2 << std::endl;

				data.event_start_idx = e1;
				data.event_stop_idx = e2;
				//out.push_back(data);
				hmm_data_ready = true;
				break;

			}  
		}
	}

	if (!hmm_data_ready)
	{
		std::cout << "The hmm data could not have been prepared correctly. Aborting." << std::endl;
		std::abort();
	}

	return(data);

}







//added by dorukb
HMMInputData AlignmentDB::get_given_event_subsequences_for_record(int e1, int e2, 
						const SequenceAlignmentRecordInfo& seq_align_info) const
{
	bool hmm_data_ready = false;
	HMMInputData data;
	//std::cout << "m_event_records.size() = " << m_event_records.size() << std::endl;
	for(size_t i = 0; i < m_event_records.size(); ++i) {
		const EventAlignmentRecord& record = m_event_records[i];
		if(record.aligned_events.empty()) {
			std::cout << "record's aligned events is empty." << std::endl;
			continue;
		}

		if(!record.sr->has_events_for_strand(record.strand)) {
			std::cout << "record doesn't have events for the strand." << std::endl;
			continue;
		}

		data.read = record.sr;
		data.pore_model = record.sr->get_base_model(record.strand);
		data.strand = record.strand;
		data.rc = record.rc;
		data.event_stride = record.stride;
		
		//int e1,e2;
		//std::cout << "record.aligned_events.size(): " << record.aligned_events.size() << std::endl;
	  
		assert(e1 >= 0);
		assert(e2 >= 0);
		//std::cout << "e1 is: " << e1 << std::endl;
		//std::cout << "e2 is: " << e2 << std::endl;

		data.event_start_idx = e1;
		data.event_stop_idx = e2;
		hmm_data_ready = true;
		break;		
		
	}


	if (!hmm_data_ready)
	{
		std::cout << "The hmm data could not have been prepared correctly. Aborting." << std::endl;
		std::abort();
	}

	return(data);

}

//added by dorukb
HMMInputData AlignmentDB::get_given_event_subsequences_for_record_for_events(int e1, int e2, 
						const SequenceAlignmentRecordInfo& seq_align_info) const
{
	bool hmm_data_ready = false;
	HMMInputData data;
	//std::cout << "m_event_records.size() = " << m_event_records.size() << std::endl;
	//for(size_t i = 0; i < m_event_records.size(); ++i) {
		//const EventAlignmentRecord& record = m_event_records[i];

		SquiggleRead* sr;
		if (m_squiggle_read_map.find(seq_align_info.read_name) != m_squiggle_read_map.end())
		{
			//std::cout << "read is found when preparing HMMInput data." << std::endl;
			sr = m_squiggle_read_map.at(seq_align_info.read_name);
		}
		else
		{
			std::cout << "read is NOT found when calling for preparing HMMinput data." << std::endl;
			std::abort();
		}



		const EventAlignmentRecord& record = EventAlignmentRecord(sr, 0, seq_align_info);
		if(record.aligned_events.empty()) {
			std::cout << "record's aligned events is empty." << std::endl;
			//continue;
			std::abort();
		}

		if(!record.sr->has_events_for_strand(record.strand)) {
			std::cout << "record doesn't have events for the strand." << std::endl;
			//continue;
			std::abort();
		}

		data.read = record.sr;
		data.pore_model = record.sr->get_base_model(record.strand);
		data.strand = record.strand;
		data.rc = record.rc;
		data.event_stride = record.stride;
		
		//int e1,e2;
		//std::cout << "record.aligned_events.size(): " << record.aligned_events.size() << std::endl;
	  
		assert(e1 >= 0);
		assert(e2 >= 0);
		//std::cout << "e1 is: " << e1 << std::endl;
		//std::cout << "e2 is: " << e2 << std::endl;

		data.event_start_idx = e1;
		data.event_stop_idx = e2;
		hmm_data_ready = true;
		//break;		
		
	//}


	if (!hmm_data_ready)
	{
		std::cout << "The hmm data could not have been prepared correctly. Aborting." << std::endl;
		std::abort();
	}

	return(data);

}




// added by dorukb for testing purposes.

std::vector<HMMInputData> AlignmentDB::get_given_event_subsequences(int e1, int e2) const
{
   
	std::vector<HMMInputData> out;
	//std::cout << "m_event_records.size() = " << m_event_records.size() << std::endl;
	for(size_t i = 0; i < m_event_records.size(); ++i) {
		const EventAlignmentRecord& record = m_event_records[i];
		if(record.aligned_events.empty()) {
			std::cout << "record's aligned events is empty." << std::endl;
			continue;
		}

		if(!record.sr->has_events_for_strand(record.strand)) {
			std::cout << "record doesn't have events for the strand." << std::endl;
			continue;
		}

		HMMInputData data;
		data.read = record.sr;
		data.pore_model = record.sr->get_base_model(record.strand);
		data.strand = record.strand;
		data.rc = record.rc;
		data.event_stride = record.stride;
		
		//int e1,e2;
		//std::cout << "record.aligned_events.size(): " << record.aligned_events.size() << std::endl;
	  
		assert(e1 >= 0);
		assert(e2 >= 0);
		//std::cout << "e1 is: " << e1 << std::endl;
		//std::cout << "e2 is: " << e2 << std::endl;

		data.event_start_idx = e1;
		data.event_stop_idx = e2;
		out.push_back(data);
			
		
	}

	return out;
}














std::vector<HMMInputData> AlignmentDB::get_event_subsequences(const std::string& contig,
															  int start_position,
															  int stop_position) const
{
	assert(m_region_contig == contig);
	assert(m_region_start <= start_position);
	assert(m_region_end >= stop_position);
	// std::cout << "m_region_contig: " << m_region_contig << std::endl;
	// std::cout << "contig: " << contig << std::endl;

	// std::cout << "m_region_start: " << m_region_start << std::endl;
	// std::cout << "start_position: " << start_position << std::endl;

	// std::cout << "m_region_end: " << m_region_end << std::endl;
	// std::cout << "stop_position: " << stop_position << std::endl;

	std::vector<HMMInputData> out;
	//std::cout << "m_event_records.size() = " << m_event_records.size() << std::endl;
	for(size_t i = 0; i < m_event_records.size(); ++i) {
		const EventAlignmentRecord& record = m_event_records[i];
		if(record.aligned_events.empty()) {
			std::cout << "record's aligned events is empty." << std::endl;
			continue;
		}

		if(!record.sr->has_events_for_strand(record.strand)) {
			std::cout << "record doesn't have events for the strand." << std::endl;
			continue;
		}

		HMMInputData data;
		data.read = record.sr;
		data.pore_model = record.sr->get_base_model(record.strand);
		data.strand = record.strand;
		data.rc = record.rc;
		data.event_stride = record.stride;
		
		int e1,e2;
		//std::cout << "record.aligned_events.size(): " << record.aligned_events.size() << std::endl;
		bool bounded = _find_by_ref_bounds(record.aligned_events, 
										   start_position, 
										   stop_position,
										   e1,
										   e2);
		//std::cout << "bounded is: " << bounded << std::endl;

		// // HARD=CODING
		// e1 = 65360;
		// e2 = 65319;

		//  e1 = 22185;
		// e2 = 22134;

		if(bounded) {
			double ratio = fabs(e1 - e2) / fabs(stop_position - start_position);
			// std::cout << "ratio is: " << ratio << std::endl;
			// std::cout << "MAX_EVENT_TO_BP_RATIO is: " << MAX_EVENT_TO_BP_RATIO << std::endl;
			// Some low quality reads appear to have "stuck" states where you get 
			// a long run of consecutive stays. They can cause an assertion in the HMM
			// so we added this heuristic to catch these.
			if(ratio < MAX_EVENT_TO_BP_RATIO) {
				assert(e1 >= 0);
				assert(e2 >= 0);
				// std::cout << "e1 is: " << e1 << std::endl;
				//  std::cout << "e2 is: " << e2 << std::endl;

				data.event_start_idx = e1;
				data.event_stop_idx = e2;
				out.push_back(data);
			}  
		}
	}



	return out;
}

std::vector<HMMInputData> AlignmentDB::get_events_aligned_to(const std::string& contig,
															 int position) const
{
	assert(m_region_contig == contig);
	assert(m_region_start <= position);
	assert(m_region_end >= position);

	std::vector<HMMInputData> out;
	for(size_t i = 0; i < m_event_records.size(); ++i) {
		const EventAlignmentRecord& record = m_event_records[i];
		if(record.aligned_events.empty()) {
			continue;
		}

		if(!record.sr->has_events_for_strand(record.strand)) {
			continue;
		}

		HMMInputData data;
		data.read = record.sr;
		data.strand = record.strand;
		data.rc = record.rc;
		data.event_stride = record.stride;
	
		AlignedPairConstIter start_iter;
		AlignedPairConstIter stop_iter;
		bool bounded = _find_iter_by_ref_bounds(record.aligned_events, position, position, start_iter, stop_iter);
		if(bounded && start_iter->ref_pos == position) {
			data.event_start_idx = start_iter->read_pos;
			data.event_stop_idx = start_iter->read_pos;
			out.push_back(data);
		}
	}
	return out;
}

std::vector<Variant> AlignmentDB::get_variants_in_region(const std::string& contig,
														 int start_position,
														 int stop_position,
														 double min_frequency,
														 int min_depth) const
{
	std::vector<Variant> variants;
	std::map<std::string, std::pair<Variant, int>> map;
	assert(stop_position >= start_position);
	const int MIN_DISTANCE_TO_REGION_END = 1;

	std::vector<int> depth(stop_position - start_position + 1, 0);

	for(size_t i = 0; i < m_sequence_records.size(); ++i) {
		const SequenceAlignmentRecord& record = m_sequence_records[i];
		if(record.aligned_bases.empty())
			continue;

		AlignedPairConstIter start_iter;
		AlignedPairConstIter stop_iter;
		_find_iter_by_ref_bounds(record.aligned_bases, start_position, stop_position, start_iter, stop_iter);

		// Increment the depth over this region
		int depth_start = start_iter->ref_pos;
		int depth_end = stop_iter == record.aligned_bases.end() ?
			record.aligned_bases.back().ref_pos : stop_iter->ref_pos;

		// clamp
		depth_start = std::max(depth_start, start_position);
		depth_end = std::min(depth_end, stop_position);

		for(; depth_start < depth_end; ++depth_start) {
			assert(depth_start >= start_position);
			assert(depth_start - start_position < depth.size());
			depth[depth_start - start_position]++;
		}

		//printf("[%zu] iter: [%d %d] [%d %d] first: %d last: %d\n", i, start_iter->ref_pos, start_iter->read_pos, stop_iter->ref_pos, stop_iter->read_pos,
		//            record.aligned_bases.front().ref_pos, record.aligned_bases.back().ref_pos);

		// Find the boundaries of a matching region
		while(start_iter != stop_iter) {
			// skip out-of-range
			int rp = start_iter->ref_pos;
			if(rp < start_position || rp > stop_position) {
				continue;
			}

			char rb = m_region_ref_sequence[start_iter->ref_pos - m_region_start];
			char ab = record.sequence[start_iter->read_pos];

			bool is_mismatch = rb != ab;
			auto next_iter = start_iter + 1;

			bool is_gap = next_iter != stop_iter &&
							(next_iter->ref_pos != start_iter->ref_pos + 1 ||
								next_iter->read_pos != start_iter->read_pos + 1);

			if(is_gap) {
				// advance the next iterator until a match is found
				while(next_iter != stop_iter) {
					char n_rb = m_region_ref_sequence[next_iter->ref_pos - m_region_start];
					char n_ab = record.sequence[next_iter->read_pos];
					if(n_rb == n_ab) {
						break;
					}
					++next_iter;
				}
			}

			// Make sure this is a variant, that it did not go off the end of the reference and that
			// it is not too close to the end of the region
			if(next_iter != stop_iter && (is_mismatch || is_gap) && next_iter->ref_pos < stop_position - MIN_DISTANCE_TO_REGION_END) {
				Variant v;
				v.ref_name = contig;
				v.ref_position = start_iter->ref_pos;

				size_t ref_sub_start = start_iter->ref_pos - m_region_start;
				size_t ref_sub_end = next_iter->ref_pos - m_region_start;
				v.ref_seq = m_region_ref_sequence.substr(ref_sub_start, ref_sub_end - ref_sub_start);
				v.alt_seq = record.sequence.substr(start_iter->read_pos, next_iter->read_pos - start_iter->read_pos);

				std::string key = v.key();
				auto iter = map.find(key);
				if(iter == map.end()) {
					map.insert(std::make_pair(key, std::make_pair(v, 1)));
				} else {
					iter->second.second += 1;
				}
			}
			start_iter = next_iter;
		}
	}

	for(auto iter = map.begin(); iter != map.end(); ++iter) {
		Variant& v = iter->second.first;
		size_t count = iter->second.second;
		size_t d = depth[v.ref_position - start_position];
		double f = (double)count / d;
		if(f >= min_frequency && (int)d >= min_depth) {
			v.add_info("BaseCalledReadsWithVariant", count);
			v.add_info("BaseCalledFraction", f);
			variants.push_back(v);
		}
	}

	std::sort(variants.begin(), variants.end(), sortByPosition);
	return variants;
}

//modified by dorukb
void AlignmentDB::load_region(const std::string& contig,
							  int start_position,
							  int stop_position)
{
	// load reference fai file
	faidx_t *fai = fai_load(m_reference_file.c_str());

	// Adjust end position to make sure we don't go out-of-range
	m_region_contig = contig;
	m_region_start = start_position;
	int contig_length = faidx_seq_len(fai, contig.c_str());
	if(contig_length == -1) {
		fprintf(stderr, "Error: could not retrieve length of contig %s from the faidx.\n", contig.c_str());
		exit(EXIT_FAILURE);
	}
 
	m_region_end = std::min(stop_position, contig_length);
	
	assert(!m_region_contig.empty());
	assert(m_region_start >= 0);
	assert(m_region_end >= 0);

	// load the reference sequence for this region
	// its ok to use the unthreadsafe fetch_seq here since we have our own fai
	int fetched_len = 0;
	char* ref_segment = faidx_fetch_seq(fai, m_region_contig.c_str(), m_region_start, m_region_end, &fetched_len);
	m_region_ref_sequence = ref_segment;
	
	// load base-space alignments
	m_sequence_records = _load_sequence_by_region(m_sequence_bam);
	//m_sequence_records_info = _load_sequence_info_by_region(m_sequence_bam);

	// load event-space alignments, possibly inferred from the base-space alignments
	if(m_event_bam.empty()) {
		m_event_records = _load_events_by_region_from_read(m_sequence_records);
	} else {
		m_event_records = _load_events_by_region_from_bam(m_event_bam);
	}

	// If an alternative basecall set was provided, load it
	// intentially overwriting the current records
	if(!m_alternative_basecalls_bam.empty()) {
		m_sequence_records = _load_sequence_by_region(m_alternative_basecalls_bam);
	}

	//_debug_print_alignments();

	free(ref_segment);
	fai_destroy(fai);
}



//added by dorukb
void AlignmentDB::load_region_select_reads(const std::string& contig,
							  int start_position,
							  int stop_position, 
							  std::vector<std::string> candidate_readnames)
{
	// load reference fai file
	faidx_t *fai = fai_load(m_reference_file.c_str());

	// Adjust end position to make sure we don't go out-of-range
	m_region_contig = contig;
	m_region_start = start_position;
	int contig_length = faidx_seq_len(fai, contig.c_str());
	if(contig_length == -1) {
		fprintf(stderr, "Error: could not retrieve length of contig %s from the faidx.\n", contig.c_str());
		exit(EXIT_FAILURE);
	}
 
	m_region_end = std::min(stop_position, contig_length);
	
	assert(!m_region_contig.empty());
	assert(m_region_start >= 0);
	assert(m_region_end >= 0);

	// load the reference sequence for this region
	// its ok to use the unthreadsafe fetch_seq here since we have our own fai
	int fetched_len = 0;
	char* ref_segment = faidx_fetch_seq(fai, m_region_contig.c_str(), m_region_start, m_region_end, &fetched_len);
	m_region_ref_sequence = ref_segment;
	
	// load base-space alignments
	m_sequence_records = _load_sequence_by_region(m_sequence_bam);
	//m_sequence_records_info = _load_sequence_info_by_region(m_sequence_bam);

	// load event-space alignments, possibly inferred from the base-space alignments
	if(m_event_bam.empty()) {
		//std::cout << "Seriously m+event+bam is empty and we load from read." << std::endl;
		//m_event_records = _load_events_by_region_from_read(m_sequence_records);
		m_event_records = _load_events_by_region_from_select_reads(m_sequence_records, candidate_readnames);
	} else {
		m_event_records = _load_events_by_region_from_bam(m_event_bam);
	}

	// If an alternative basecall set was provided, load it
	// intentially overwriting the current records
	if(!m_alternative_basecalls_bam.empty()) {
		m_sequence_records = _load_sequence_by_region(m_alternative_basecalls_bam);
	}

	//_debug_print_alignments();

	free(ref_segment);
	fai_destroy(fai);
}

void AlignmentDB::_clear_region()
{
	// Delete the SquiggleReads
	for(SquiggleReadMap::iterator iter = m_squiggle_read_map.begin();
		iter != m_squiggle_read_map.end(); ++iter) 
	{
		delete iter->second;
		iter->second = NULL;
	}

	m_squiggle_read_map.clear();
	m_sequence_records.clear();
	m_event_records.clear();

	m_region_contig = "";
	m_region_start = -1;
	m_region_end = -1;
}



std::vector<SequenceAlignmentRecord> AlignmentDB::_load_sequence_by_region(const std::string& sequence_bam)
{
	assert(!m_region_contig.empty());
	assert(m_region_start >= 0);
	assert(m_region_end >= 0);

	BamHandles handles = _initialize_bam_itr(sequence_bam, m_region_contig, m_region_start, m_region_end);
	std::vector<SequenceAlignmentRecord> records;

	int result;
	while((result = sam_itr_next(handles.bam_fh, handles.itr, handles.bam_record)) >= 0) {


		// commented out by dorukb
		// skip ambiguously mapped reads
		//std::cout << "loading a sequence... " << std::endl;
		if(handles.bam_record->core.qual < 1) {
			continue;
		}

		records.emplace_back(handles.bam_record);

		//added by dorukb
		//m_sequence_records_info.emplace_back(handles.bam_record, m_bamHdr);
	}

	// cleanup
	sam_itr_destroy(handles.itr);
	bam_destroy1(handles.bam_record);
	sam_close(handles.bam_fh);
	
	return records;
}

std::vector<EventAlignmentRecord> AlignmentDB::_load_events_by_region_from_bam(const std::string& event_bam)
{
	//std::cout << "loading events by region from bam" << std::endl;
	BamHandles handles = _initialize_bam_itr(event_bam, m_region_contig, m_region_start, m_region_end);

	std::vector<EventAlignmentRecord> records;

	int result;
	while((result = sam_itr_next(handles.bam_fh, handles.itr, handles.bam_record)) >= 0) {
		EventAlignmentRecord event_record;

		std::string full_name = bam_get_qname(handles.bam_record);
		
		// Check for the template/complement suffix
		bool is_template = true;
		size_t suffix_pos = 0;
		suffix_pos = full_name.find(".template");
		if(suffix_pos == std::string::npos) {
			suffix_pos = full_name.find(".complement");
			assert(suffix_pos != std::string::npos);
			is_template = false;
		}

		std::string read_name = full_name.substr(0, suffix_pos);
		_load_squiggle_read(read_name);
		event_record.sr = m_squiggle_read_map[read_name];

		// extract the event stride tag which tells us whether the
		// event indices are increasing or decreasing
		assert(bam_aux_get(handles.bam_record, "ES") != NULL);
		int event_stride = bam_aux2i(bam_aux_get(handles.bam_record, "ES"));

		// copy event alignments
		std::vector<AlignedSegment> alignments = get_aligned_segments(handles.bam_record, event_stride);
		assert(alignments.size() > 0);
		event_record.aligned_events = alignments[0];

		event_record.rc = bam_is_rev(handles.bam_record);
		event_record.stride = event_stride;
		event_record.strand = is_template ? T_IDX : C_IDX;
		records.push_back(event_record);

		/*
		printf("event_record[%zu] name: %s stride: %d align bounds [%d %d] [%d %d]\n", 
			m_event_records.size() - 1,
			bam_get_qname(handles.bam_record),
			event_stride,
			m_event_records.back().aligned_events.front().ref_pos,
			m_event_records.back().aligned_events.front().read_pos,
			m_event_records.back().aligned_events.back().ref_pos,
			m_event_records.back().aligned_events.back().read_pos);
		*/
	}

	// cleanup
	sam_itr_destroy(handles.itr);
	bam_destroy1(handles.bam_record);
	sam_close(handles.bam_fh);
	
	return records;
}


// EventAlignmentRecord AlignmentDB::load_events_from_single_read(const SequenceAlignmentRecord& seq_record)
// {
//     EventAlignmentRecord record;

//     // conditionally load the squiggle read if it hasn't been loaded already
//     _load_squiggle_read(seq_record.read_name);
//     SquiggleRead* sr = m_squiggle_read_map[seq_record.read_name];

//     bool events_loaded = false;
//     for(size_t si = 0; si < NUM_STRANDS; ++si) {

//         if ( (si != C_IDX) && (sr->has_events_for_strand(si)) )
//         {
//             record = EventAlignmentRecord(sr, si, seq_record);
//             events_loaded = true;
//             break;
//         }
//     }

//     if (!events_loaded)
//     {
//         std::cout << "event have not been loaded. Abort." << std::endl;
//         std::abort();
//     }

//     return(record);


// }

std::vector<EventAlignmentRecord> AlignmentDB::_load_events_by_region_from_read(const std::vector<SequenceAlignmentRecord>& sequence_records)
{
	//std::cout << "loading events by region from READ" << std::endl;
	std::vector<EventAlignmentRecord> records;
	for(size_t i = 0; i < sequence_records.size(); ++i) {
		const SequenceAlignmentRecord& seq_record = sequence_records[i];

		// conditionally load the squiggle read if it hasn't been loaded already
		_load_squiggle_read(seq_record.read_name);
		for(size_t si = 0; si < NUM_STRANDS; ++si) {
			
			// skip complement
			if(si == C_IDX) {
				continue;
			}
	
			// skip reads that do not have events here
			SquiggleRead* sr = m_squiggle_read_map[seq_record.read_name];
			if(!sr->has_events_for_strand(si)) {
				continue;
			}

			records.emplace_back(sr, si, seq_record);
		}
	}

	return records;
}


std::vector<EventAlignmentRecord> AlignmentDB::_load_events_by_region_from_select_reads(const std::vector<SequenceAlignmentRecord>& sequence_records, const std::vector<std::string>& candidate_readnames )
{
	//std::cout << "loading events by region from SELECT READS" << std::endl;
	std::vector<EventAlignmentRecord> records;
	for(size_t i = 0; i < sequence_records.size(); ++i) {

		const SequenceAlignmentRecord& seq_record = sequence_records[i];
		// conditionally load the squiggle read if it hasn't been loaded already
		if (std::binary_search(candidate_readnames.begin(), candidate_readnames.end(), seq_record.read_name))
		{
			_load_squiggle_read(seq_record.read_name);
			for(size_t si = 0; si < NUM_STRANDS; ++si) 
			{
				// skip complement
				if(si == C_IDX) {
					continue;
				}
		
				// skip reads that do not have events here
				SquiggleRead* sr = m_squiggle_read_map[seq_record.read_name];
				if(!sr->has_events_for_strand(si)) {
					continue;
				}
				records.emplace_back(sr, si, seq_record);
			}
			
		}

		





	}

	return records;
}




void AlignmentDB::_debug_print_alignments()
{
	// Build a map from a squiggle read to the middle base of the reference region it aligned to
	std::map<std::string, int> read_to_ref_middle;
	for(size_t i = 0; i < m_sequence_records.size(); ++i) {
		const SequenceAlignmentRecord& record = m_sequence_records[i];
		int middle = record.aligned_bases[record.aligned_bases.size()/2].ref_pos;
		read_to_ref_middle.insert(std::make_pair(record.read_name, middle));
	}

	for(size_t i = 0; i < m_event_records.size(); ++i) {
		const EventAlignmentRecord& record = m_event_records[i];
		AlignedPair first = record.aligned_events.front();

		int ref_middle = read_to_ref_middle[record.sr->read_name];
		int event_middle_start, event_middle_end;
		_find_by_ref_bounds(record.aligned_events, ref_middle, ref_middle, event_middle_start, event_middle_end);
		AlignedPair last = record.aligned_events.back();
		printf("event_record[%zu] name: %s strand: %d stride: %d rc: %d align bounds [%d %d] [%d %d] [%d %d]\n", 
				i,
				record.sr->read_name.c_str(),
				record.strand,
				record.stride,
				record.rc,
				first.ref_pos,
				first.read_pos,
				ref_middle,
				event_middle_start,
				last.ref_pos,
				last.read_pos);
	}
}

void AlignmentDB::_load_squiggle_read(const std::string& read_name)
{
	// Do we need to load this fast5 file?
	//std::cout << "loading the squiggle?" << std::endl;
	if(m_squiggle_read_map.find(read_name) == m_squiggle_read_map.end()) {

		SquiggleRead* sr = new SquiggleRead(read_name, m_read_db);
		//std::cout << "Squiggle loaded." << std::endl;
		m_squiggle_read_map[read_name] = sr;
	}
	//std::cout << "done." << std::endl;
} 

std::vector<EventAlignment> AlignmentDB::_build_event_alignment(const EventAlignmentRecord& event_record) const
{
	std::vector<EventAlignment> alignment;
	const SquiggleRead* sr = event_record.sr;
	size_t k = sr->get_model_k(event_record.strand);

	for(const auto& ap : event_record.aligned_events) { 

		EventAlignment ea;
		ea.ref_position = ap.ref_pos;
		if(ea.ref_position < m_region_start || ea.ref_position >= m_region_end - k) {
			continue;
		}

		ea.event_idx = ap.read_pos;

		std::string kmer = get_reference_substring(m_region_contig, ea.ref_position, ea.ref_position + k - 1);
		assert(kmer.size() == k);

		// ref data
		ea.ref_name = "read"; // not needed
		ea.read_idx = -1; // not needed
		ea.ref_kmer = kmer;
		ea.strand_idx = event_record.strand;
		ea.rc = event_record.rc;
		ea.model_kmer = kmer;
		ea.hmm_state = 'M';
		alignment.push_back(ea);
	}

	return alignment;
}

bool AlignmentDB::_find_read_pos_from_ref_pos(const std::vector<AlignedPair>& pairs,
									  int ref_pos,
									  int& read_pos)
{
	AlignedPairRefLBComp lb_comp;
	//AlignedPairRefUBComp ub_comp;

	auto start_iter = std::lower_bound(pairs.begin(), pairs.end(),
								  ref_pos, lb_comp);

	// auto stop_iter = std::upper_bound(pairs.begin(), pairs.end(),
	// 							 ref_pos, ub_comp);


	

	if (start_iter == pairs.end())
	{
		return false;
	}
	read_pos = start_iter->read_pos;
	return true;



	//std::cout << "start_iter->read_pos: " << start_iter->read_pos << std::endl;
	//std::cout << "stop_iter->ref_pos: " << stop_iter->ref_pos << std::endl;
	// else
	// {
	// 	read_pos = -1;
	// 	return false;
	// }
	
}

bool AlignmentDB::_find_ref_pos_from_read_pos(const std::vector<AlignedPair>& pairs,
									  int read_pos,
									  int& ref_pos)
{
	AlignedPairReadLBComp lb_comp;
	//AlignedPairReadUBComp ub_comp;

	auto start_iter = std::lower_bound(pairs.begin(), pairs.end(),
								  read_pos, lb_comp);

	// auto stop_iter = std::upper_bound(pairs.begin(), pairs.end(),
	// 							 read_pos, ub_comp);


	if (start_iter == pairs.end())
	{
		std::cout << "the read position must exist. Aborting." << std::endl;
		std::abort();
	}

	ref_pos = start_iter->ref_pos;
	return true;
	


	// std::cout << "just passed read_pos: " << read_pos << std::endl;
	// std::cout << "start_iter->read_pos: " << start_iter->read_pos << std::endl;
	// std::cout << "start_iter->ref_pos: " << start_iter->ref_pos << std::endl;
	// std::cout << "stop_iter->read_pos: " << stop_iter->read_pos << std::endl;
	// else
	// {
	// 	read_pos = -1;
	// 	return false;
	// }
	
}




bool AlignmentDB::_find_iter_by_ref_bounds(const std::vector<AlignedPair>& pairs,
									  int ref_start,
									  int ref_stop,
									  AlignedPairConstIter& start_iter,
									  AlignedPairConstIter& stop_iter)
{
	AlignedPairRefLBComp lb_comp;
	start_iter = std::lower_bound(pairs.begin(), pairs.end(),
								  ref_start, lb_comp);

	stop_iter = std::lower_bound(pairs.begin(), pairs.end(),
								 ref_stop, lb_comp);

	// std::cout << "start_iter->ref_pos: " << start_iter->ref_pos << std::endl;
	// std::cout << "(start_iter - 1)->ref_pos: " << (start_iter - 1)->ref_pos << std::endl;
	// std::cout << "ref_start: " <<ref_start << std::endl;
	// std::cout << "pairs[0].ref_pos: " <<pairs[0].ref_pos << std::endl;


	// std::cout << "stop_iter->ref_pos: " << stop_iter->ref_pos << std::endl;
	// std::cout << "(stop_iter + 1)->ref_pos: " << (stop_iter + 1)->ref_pos << std::endl;
	// std::cout << "ref_stop: " << ref_stop << std::endl;
	// std::cout << "pairs[pairs.size()-1].ref_pos: " << pairs[pairs.size()-1].ref_pos << std::endl;
	
	if(start_iter == pairs.end() || stop_iter == pairs.end())
	{
		std::cout << "iter ends are not bounded." << std::endl;
		return false;
	}



	// require at least one aligned reference base at or outside the boundary
	bool left_bounded = start_iter->ref_pos <= ref_start ||
						(start_iter != pairs.begin() && (start_iter - 1)->ref_pos <= ref_start);
	
	bool right_bounded = stop_iter->ref_pos >= ref_stop ||
						(stop_iter != pairs.end() && (stop_iter + 1)->ref_pos >= ref_start);

	return left_bounded && right_bounded;
}


bool AlignmentDB::_find_by_ref_bounds(const std::vector<AlignedPair>& pairs,
									  int ref_start,
									  int ref_stop,
									  int& read_start,
									  int& read_stop)
{
	AlignedPairConstIter start_iter;
	AlignedPairConstIter stop_iter;
	bool bounded = _find_iter_by_ref_bounds(pairs, ref_start, ref_stop, start_iter, stop_iter);
	if(bounded) {
		read_start = start_iter->read_pos;
		read_stop = stop_iter->read_pos;
		return true;
	} else {
		return false;
	}
}
