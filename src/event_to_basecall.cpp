

#include <assert.h>
#include <algorithm>
#include <iostream>
#include <string>

#include "nanopolish_alignment_db.h"
#include "event_to_basecall.h"




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


  

bool map_events_to_basecall(const SquiggleRead * sr, std::vector<int>& event_indices_for_bases)
{


	std::string basecall = sr->read_sequence;  
	
	event_indices_for_bases.resize(basecall.length(), -1);

	int kmer_size = 5;
	bool success = false;
	for(size_t si = 0; si < NUM_STRANDS; ++si) 
	{
		
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
						//std::cout << sequence_record.read_name << " is the read name."  <<  std::endl;
						//std::cout << sequence_record.chromosome << " is the chromosome."  <<  std::endl;
						//std::cout << sequence_record.beginPos << " is the begin pos"  <<  std::endl;
						std::cout <<  basecall.substr(0, 100)  <<  " IS THE seq." << std::endl;
						//std::cout << sequence_record.sequence_len << " is the seq len"  <<  std::endl;
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