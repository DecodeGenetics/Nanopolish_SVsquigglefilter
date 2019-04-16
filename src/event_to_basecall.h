//---------------------------------------------------------
// Copyright 2018 deCODE genetics
// Written by Doruk Beyter (dorukb@decode.is)
//---------------------------------------------------------
//
// nanopolish_call_variants -- find variants wrt a reference
//
#ifndef __EVENT_TO_BASECALL_H
#define __EVENT_TO_BASECALL_H


#include <string.h>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

char complement(char n);
bool is_homopolymer(const std::string& s);
void fill_indices_to_missed_bases(std::vector<int>& base_event_indices, int kmer_size);
void get_begin_end_event_indices_for_read_region(const SequenceAlignmentRecordInfo& seq_record, const std::vector<int>& event_inds_for_bases, 
                                    const int k_mer_length, const int read_start, const int read_end, int& event_begin_idx, int& event_end_idx);
bool map_events_to_basecall(const SquiggleRead * sr, std::vector<int>& event_indices_for_bases);







#endif
