//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_alignment_db -- abstraction for working
// with sets of reads/events aligned to a reference genome
//
#ifndef ALIGNMENT_DB
#define ALIGNMENT_DB

#include <string>
#include <vector>
#include <map>
#include "nanopolish_anchor.h"
#include "nanopolish_variant.h"

#define MAX_EVENT_TO_BP_RATIO 20

// structs


struct SequenceAlignmentRecordInfo
{
    SequenceAlignmentRecordInfo(const bam1_t* record, const bam_hdr_t *bamHdr);

    long beginPos;
    long endPos;
    std::string chromosome;
    std::string read_name;

    std::string sequence;
    int sequence_len;
    std::vector<AlignedPair> aligned_bases;
    
    uint8_t rc; // with respect to reference genome
    std::vector<unsigned int> base_to_event_map;

    int map_quality;

};


struct SequenceAlignmentRecord
{
    SequenceAlignmentRecord(const bam1_t* record);

    std::string read_name;
    std::string sequence;
    std::vector<AlignedPair> aligned_bases;
    
    uint8_t rc; // with respect to reference genome
};

struct EventAlignmentRecord
{
    EventAlignmentRecord() {}
    EventAlignmentRecord(SquiggleRead* sr,
                         const int strand_idx,
                         const SequenceAlignmentRecord& seq_record);
    EventAlignmentRecord(SquiggleRead* sr,
                         const int strand_idx,
                         const SequenceAlignmentRecordInfo& seq_record);

    SquiggleRead* sr;
    uint8_t rc; // with respect to reference genome
    uint8_t strand; // 0 = template, 1 = complement
    int stride; // whether event indices increase or decrease along the reference
    std::vector<AlignedPair> aligned_events;
    std::vector<AlignedPair> aligned_readpos_events;
    std::string read_name;
};

// typedefs
typedef std::map<std::string, SquiggleRead*> SquiggleReadMap;






class AlignmentDB
{
    public:
        AlignmentDB(const std::string& reads_file,
                    const std::string& reference_file,
                    const std::string& sequence_bam,
                    const std::string& event_bam);

        ~AlignmentDB();

        void load_region(const std::string& contig,
                         int start_position,
                         int stop_position);


        void load_region_select_reads(const std::string& contig,
                              int start_position,
                              int stop_position, 
                              std::vector<std::string> candidate_readnames);
    
        // Some high quality basecallers, like scrappie, may not output event
        // annotations. This call is to support using scrappie basecalls
        // with metrichor events.
        void set_alternative_basecalls_bam(const std::string& alternative_sequence_bam);
        
        const std::string& get_reference() const { return m_region_ref_sequence; }

        bool are_coordinates_valid(const std::string& contig,
                                   int start_position,
                                   int stop_position) const;

        std::string get_reference_substring(const std::string& contig,
                                            int start_position,
                                            int stop_position) const;

        std::vector<std::string> get_read_substrings(const std::string& contig,
                                                     int start_position,
                                                     int stop_position) const;

        std::vector<HMMInputData> get_event_subsequences(const std::string& contig,
                                                         int start_position,
                                                         int stop_position) const;

        //added by dorukb
        HMMInputData get_given_event_subsequences_for_record_for_events(int e1, int e2, 
                        const SequenceAlignmentRecordInfo& seq_align_info) const;

        //added by dorukb
        bool find_event_coords_for_region_for_record(const SequenceAlignmentRecordInfo& sequence_record, 
                                        const int start_pos, const int stop_pos, int& event_begin_idx, int& event_end_idx);



        //added by dorukb
        bool find_event_coords_for_given_read_coords(const SequenceAlignmentRecordInfo& sequence_record, 
                                        const int r1, const int r2, int& event_begin_idx, int& event_end_idx) const;

     
        std::vector<HMMInputData> get_events_aligned_to(const std::string& contig, int position) const;

        std::vector<Variant> get_variants_in_region(const std::string& contig,
                                                    int start_position,
                                                    int stop_position,
                                                    double min_frequency,
                                                    int min_depth) const;

        const std::vector<EventAlignmentRecord>& get_eventalignment_records() const { return m_event_records; }

        // reference metadata
        std::string get_region_contig() const { return m_region_contig; }
        int get_region_start() const { return m_region_start; }
        int get_region_end() const { return m_region_end; }
        
        void set_alternative_model_type(const std::string model_type_string) { m_model_type_string = model_type_string; }
        
        // Search the vector of AlignedPairs using lower_bound/upper_bound
        // and the input reference coordinates. If the search succeeds,
        // set read_start/read_stop to be the read_pos of the bounding elements
        // and return true. 
        static bool _find_by_ref_bounds(const std::vector<AlignedPair>& pairs,
                                 int ref_start,
                                 int ref_stop,
                                 int& read_start,
                                 int& read_stop);

        static bool _find_iter_by_ref_bounds(const std::vector<AlignedPair>& pairs,
                                      int ref_start,
                                      int ref_stop,
                                      AlignedPairConstIter& start_iter,
                                      AlignedPairConstIter& stop_iter);

        // added by dorukb
        static bool _find_read_pos_from_ref_pos(const std::vector<AlignedPair>& pairs,
                                      int ref_pos,
                                      int& read_pos);
        // added by dorukb
        static bool _find_ref_pos_from_read_pos(const std::vector<AlignedPair>& pairs,
                                      int read_pos,
                                      int& ref_pos);
        // added by dorukb
        std::vector<SequenceAlignmentRecordInfo> load_all_sequences_info(const std::string& sequence_bam);
        
        // added by dorukb
        std::vector<SequenceAlignmentRecordInfo> load_region_sequences_info(const std::string& contig, 
                                                int start_position, int stop_position,
                                                const std::string& sequence_bam) const;


        bool find_scrappie_events_for_basecall(const SequenceAlignmentRecordInfo& sequence_record, std::vector<int>& event_indices_for_bases) const;



    private:
       
        std::vector<SequenceAlignmentRecord> _load_sequence_by_region(const std::string& sequence_bam);
        std::vector<EventAlignmentRecord> _load_events_by_region_from_bam(const std::string& event_bam);
        std::vector<EventAlignmentRecord> _load_events_by_region_from_read(const std::vector<SequenceAlignmentRecord>& sequence_records);
        std::vector<EventAlignmentRecord> _load_events_by_region_from_select_reads(const std::vector<SequenceAlignmentRecord>& sequence_records, 
                                                                                    const std::vector<std::string>& candidate_readnames );

        void _load_squiggle_read(const std::string& read_name);

        void _clear_region();

        void _debug_print_alignments();

        std::vector<EventAlignment> _build_event_alignment(const EventAlignmentRecord& event_record) const;


        //
        // data
        //
        std::string m_reference_file;
        std::string m_sequence_bam;
        std::string m_event_bam;
        std::string m_alternative_basecalls_bam;

        // parameters

        // loaded region
        std::string m_region_ref_sequence;
        std::string m_region_contig;
        int m_region_start;
        int m_region_end;

        // cached alignments for a region
        ReadDB m_read_db;
        std::vector<SequenceAlignmentRecord> m_sequence_records;
        std::vector<EventAlignmentRecord> m_event_records;
        SquiggleReadMap m_squiggle_read_map;
        std::string m_model_type_string;

        // added by dorukb 
        // cached for the header.
        //std::vector<SequenceAlignmentRecordInfo> m_sequence_records_info;
        //bam_hdr_t *m_bamHdr;
};

#endif
