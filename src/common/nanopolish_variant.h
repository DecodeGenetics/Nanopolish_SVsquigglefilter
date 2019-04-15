//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_variant -- tools for calling variants
//
#ifndef NANOPOLISH_VARIANT_H
#define NANOPOLISH_VARIANT_H

#include <sstream>
#include <string>
#include <vector>
#include <iterator>

#include "stdaln.h"
#include <regex>
#include "nanopolish_common.h"
#include "nanopolish_hmm_input_sequence.h"
#include <limits>


// forward declare
class Haplotype;
class VariantGroup;
class AlignmentDB;







// void regex_split(const std::string &s, std::string delim, std::vector <std::string> &elems)
// {
//  regex re(delim);
//  //string s = "The White Rabbit,  is very,late.";
//  std::sregex_token_iterator it(s.begin(), s.end(), re, -1);
//  std::sregex_token_iterator reg_end;
//  for (; it != reg_end; ++it) 
//  {
//       std::cout << "Regex: " << it->str() << std::endl;
//  }

// }

struct Variant
{


  std::string getFastqSuitableReadname(std::string const & record_name)
  {
    std::string refined_record_name;
    std::size_t found = record_name.find("_Basecall");
    if(found == std::string::npos)
    {
      refined_record_name = record_name;
    }
    else
    {
      refined_record_name = record_name.substr(0, found);
    }
    
    return(refined_record_name);
  }


  void split(const std::string &s, char delim, std::vector <std::string> &elems) 
  {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
      if (item.length() > 0) {
        elems.push_back(item);  
      }
    }
    return;
  }







  static void write_vcf_header(FILE* fp,  const std::vector<std::string>& tag_lines = std::vector<std::string>());

  static std::string make_vcf_tag_string(const std::string& tag,
                      const std::string& id,
                      int count,
                      const std::string& type,
                      const std::string& description);

  Variant() { }   

  //modified by dorukb
  Variant(const std::string& line, bool isSVvcf=false) 
  {  
    if (!isSVvcf)
    {
      read_vcf(line);
    }
    else
    {
      read_SVvcf(line);
    }
    
  }

  // added by dorukb
  Variant(const std::string& ref_name, const size_t& ref_position, const std::string& alt_seq, const std::string& var_seq);

  // generate a unique identifier for this variant
  std::string key() const
  {
      std::stringstream out;
      out << ref_name << ':' << ref_position << ':' << ref_seq << ':' << alt_seq;
      return out.str();
  }

  void write_vcf(FILE* fp) const
  {
      assert(fp != NULL);
      const char* gt_def = genotype.empty() ? NULL : "GT";
      const char* gt_str = genotype.empty() ? NULL : genotype.c_str();

      fprintf(fp, "%s\t%zu\t%s\t", ref_name.c_str(), ref_position + 1, ".");
      fprintf(fp, "%s\t%s\t%.1lf\t", ref_seq.c_str(), alt_seq.c_str(), quality);
      fprintf(fp, "%s\t%s\t%s\t%s\n", "PASS", info.c_str(), gt_def, gt_str);
  }

  void read_vcf(const std::string& line)
  {
    std::stringstream ss(line);
    std::string dummy;
    ss >> ref_name;
    ss >> ref_position;
    ss >> dummy; // ID, not used
    ss >> ref_seq;
    ss >> alt_seq;
    ss >> quality;
    ss >> dummy; // FILTER, not used
    ss >> info;
    ss >> dummy; // GT tag
    ss >> genotype;

    // VCF is 1-based but we internally represent a variant as 0-based
    ref_position -= 1;

    assert(!ref_name.empty());
    assert(!ref_seq.empty());
    assert(!alt_seq.empty());
    //assert(ref_position >= 0);
    assert(quality >= 0.0f);
  }

        // added by dorukb
  void read_SVvcf(const std::string& line)
  {
    std::stringstream ss(line);
    std::string dummy;
    std::string quality_str;

    ss >> ref_name;
    ss >> ref_position;
    ss >> SV_ID; 

    ss >> ref_seq;
    ss >> alt_seq;
    ss >> quality_str;

    ss >> dummy; // FILTER, not used
    ss >> info;
    ss >> dummy; // GT tag
    ss >> genotype;

    // VCF is 1-based but we internally represent a variant as 0-based
    ref_position -= 1;

    // Learn the readnames supporting the SV event. 
    std::regex rnames_regex(";RNAMES=(.*);SUPTYPE");

    std::string seq_identifier = ";SEQ=";
    std::string re_identifier = ";RE=";
    // I had to strip the INS sequence and only accept some prefix and suffix because the regex fails.
    int info_prefix_end_index =  info.find(seq_identifier) + seq_identifier.length();
    std::string info_prefix = info.substr(0, info.find(seq_identifier));

    auto regex_it = std::sregex_iterator(info_prefix.begin(), info_prefix.end(), rnames_regex);
    if (regex_it != std::sregex_iterator()) // there is a match (it should be)
    { 

      std::string rnames = ((*regex_it)[1]).str();
      std::vector<std::string> read_ids;
      split(rnames, ',', read_ids);
      for (int i = 0; i < read_ids.size(); i++)
      {
        std::string parsed_read_id = getFastqSuitableReadname(read_ids[i]);
        supporting_reads.push_back(parsed_read_id);
        //std::cout << "rname er: " << parsed_read_id << std::endl;
      } 
      std::sort(supporting_reads.begin(), supporting_reads.end());
      
    }
    else
    {
      std::cerr<< "[SV " << ref_position << "] could not find RNAMES in the info field" << std::endl;
      std::abort();
    }

    if (alt_seq == "<DEL>")
    { 
      std::regex regex_seq_info(";END=(.*);");
      auto regex_it = std::sregex_iterator(info.begin(), info.end(), regex_seq_info);
      
      
      if (regex_it != std::sregex_iterator()) // there is a match (it should be)
      {    
        //std::cout << " regex worked." << std::endl;
        std::cout << std::abs( atoi( (((*regex_it)[1]).str()).c_str()) - ref_position) << std::endl;
        std::cout << ref_position << std::endl;
        size_t end_site = atoi( (((*regex_it)[1]).str()).c_str());
        if (end_site >= ref_position)
        {
          ref_length = end_site - ref_position;
        }
        else
        {
          ref_length = ref_position - end_site;
        }
        ref_length += 1;
        
        alt_length = 1; // vcf format specification for a deletion ALT has length 1
        std::cout << " ref_length" << ref_length<< std::endl;
        std::cout << " alt_length" << alt_length << std::endl;
      } 
      else
      {
        //std::cout<< "[SV " << ref_position << "] could not find SVLEN in the info field" << std::endl;
        std::cerr<< "[SV " << ref_position << "] could not find SVLEN in the info field" << std::endl;
      }
    }
    else if (alt_seq == "<INS>")
    {

      if ( (info.find(seq_identifier) != std::string::npos) && (info.find(re_identifier) != std::string::npos))
      {
        ref_length = 1;
        int seq_begin_site = info.find(seq_identifier) + seq_identifier.length();
        int seq_length = info.find(re_identifier) - seq_begin_site;
        var_seq = info.substr(seq_begin_site, seq_length);
        alt_length = var_seq.length() + 1;
      }
      else
      {
        std::cerr << "[SV " << ref_position << "] I could not find SEQ in the info field" << std::endl;
      }
    }

    assert(!ref_name.empty());
    assert(!ref_seq.empty());
    assert(!alt_seq.empty());
    //assert(ref_position >= 0);
    //assert(quality >= 0.0f);
  }

  template<typename T>
  void add_info(const std::string& key, T value)
  {
      std::stringstream ss;
      ss << key << "=" << value;
      if(info.empty()) {
          info = ss.str();
      } else {
          info.append(1, ';');
          info.append(ss.str());
      }
  }

  bool is_snp() const { return ref_seq.length() == 1 && alt_seq.length() == 1; }

  std::string ref_name;
  size_t ref_position;
  //size_t end_position;
  std::string SV_ID;
  std::string ref_seq;
  std::string alt_seq;
  double quality = 0;
  std::string info;
  std::string genotype;
  size_t ref_length = 0;
  size_t alt_length = 0;
  std::string var_seq = "";
  std::vector <std::string> supporting_reads;
};

inline bool sortByPosition(const Variant& a, const Variant& b) 
{ 
    return a.ref_name == b.ref_name ? 
        a.ref_position < b.ref_position : 
        a.ref_name < b.ref_name; 
}

class VariantKeyComp
{
    public: 
        inline bool operator()(const Variant& a, const Variant& b)
        {
            return a.key() < b.key();
        }
};

// Read a collection of variants from a VCF file
std::vector<Variant> read_variants_from_file(const std::string& filename, bool isSVvcf=false);
std::vector<Variant> read_variants_for_region(const std::string& filename,
                                              const std::string& contig,
                                              int region_start,
                                              int region_end,
                                              bool isSVvcf=false);

// Remove variants that are in the vector fewer than min_count times
void filter_variants_by_count(std::vector<Variant>& variants, int min_count);

// Remove snps or indels 
void filter_out_non_snp_variants(std::vector<Variant>& variants);

// Expand the HMMInputSequence to contain alternatives representing
// methylated versions with the same basic sequence. The output includes
// the unmethylated version of the sequence over the nucleotide alphabet
std::vector<HMMInputSequence> generate_methylated_alternatives(const HMMInputSequence& sequence, 
                                                                const std::vector<std::string>& methylation_types);

// Score the variants contained within the input group using the nanopolish HMM
void score_variant_group(VariantGroup& variant_group,
                          Haplotype base_haplotype, 
                          const std::vector<HMMInputData>& input,
                          const int max_haplotypes,
                          const int ploidy,
                          const bool genotype_all_input_variants,
                          const uint32_t alignment_flags,
                          const std::vector<std::string>& methylation_types);

// Call genotypes for the variants in this group using a simple model
std::vector<Variant> simple_call(VariantGroup& variant_group,
                                  const int ploidy,
                                  const bool genotype_all_input_variants);

// Call genotypes for the variants in this group using support from nearby variants
std::vector<Variant> multi_call(VariantGroup& variant_group,
                                std::vector<const VariantGroup*> neighbor_groups,
                                const int ploidy,
                                const bool genotype_all_input_variants);

// Score a single variant, stopping when the absolute value of the score relative
// to the reference meets a threshold
Variant score_variant_thresholded(const Variant& input_variant,
                                  Haplotype base_haplotype, 
                                  const std::vector<HMMInputData>& input,
                                  const uint32_t alignment_flags,
                                  const uint32_t score_threshold,
                                  const std::vector<std::string>& methylation_types);

// added by dorukb
double score_variant_for_record(const Variant& input_variant,
                                  Haplotype base_haplotype, 
                                  const HMMInputData& input,
                                  const uint32_t alignment_flags,
                                  const std::vector<std::string>& methylation_types,
                                  double& base_score,
                                  double& variant_score,
                                  int& ref_seq_len,
                                  int& alt_seq_len,
                                  int& hmm_input_len);
// added by dorukb
double score_haplotype_for_record(Haplotype a_haplotype, 
                    const HMMInputData& input,
                    const uint32_t alignment_flags,
                    const std::vector<std::string>& methylation_types,
                    int& seq_len,
                    int& hmm_input_len);

// Annotate each SNP variant in the input set with the fraction of reads supporting every possible base at the position
void annotate_variants_with_all_support(std::vector<Variant>& input, const AlignmentDB& alignments, int min_flanking_sequence, const uint32_t alignment_flags);


#endif
