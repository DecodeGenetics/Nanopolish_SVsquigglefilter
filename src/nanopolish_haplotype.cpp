//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_haplotype - a haplotype derived from 
// a reference sequence and a set of variants
//
#include "nanopolish_haplotype.h"
 
// Definitions
const size_t Haplotype::INSERTED_POSITION = std::string::npos;


Haplotype::Haplotype()
{
}

Haplotype::Haplotype(const std::string& ref_name,
					 const size_t ref_position,
					 const std::string& ref_sequence) : 
						m_ref_name(ref_name),
						m_ref_position(ref_position),
						m_reference(ref_sequence)
{
	// m_sequence = m_reference;


	// line modified by dorukb
	m_sequence = replace_M_with_N_in_reference(m_reference);
	m_coordinate_map.resize(m_reference.size());
	for(size_t i = 0; i < m_coordinate_map.size(); ++i) {
		m_coordinate_map[i] = m_ref_position + i;
	}
}

Haplotype::~Haplotype()
{
}
 
// added by dorukb due to a failing assertion in nanopolish_alphabet.h
std::string Haplotype::replace_M_with_N_in_reference(const std::string& ref_sequence)
{
	std::string modified_seq = ref_sequence;
	std::replace(modified_seq.begin(), modified_seq.end(), 'M', 'N');
	return(modified_seq);

}

bool Haplotype::truncate_seq_from_right_end(int bp)
{
	if (m_sequence.length() > bp){
		size_t derived_idx = m_sequence.length() - bp;
		m_sequence.erase(derived_idx, bp);

		// fix the coordinate map accordingly now.
		// make a pair of iterators that bound the changed sequence
		std::vector<size_t>::iterator fi = m_coordinate_map.begin() + derived_idx;
		std::vector<size_t>::iterator li = fi + bp;

		// erase the positions of the changed bases
		std::vector<size_t>::iterator ii = m_coordinate_map.erase(fi, li);

		// insert new positions for the alt bases with invalid indices
		//m_coordinate_map.insert(ii, al, INSERTED_POSITION);

		// sanity check
		assert(m_coordinate_map.size() == m_sequence.size());
		return true;

	}
	else
	{
		std::cout << " You are deleting more than the haplotype's length." << std::endl;
		return false;
	}
}


bool Haplotype::truncate_seq_from_left_end(int bp)
{
	if (m_sequence.length() > bp){
		size_t derived_idx = 0;
		m_sequence.erase(derived_idx, bp);

		// fix the coordinate map accordingly now.
		// make a pair of iterators that bound the changed sequence
		std::vector<size_t>::iterator fi = m_coordinate_map.begin() + derived_idx;
		std::vector<size_t>::iterator li = fi + bp;

		// erase the positions of the changed bases
		std::vector<size_t>::iterator ii = m_coordinate_map.erase(fi, li);

		// insert new positions for the alt bases with invalid indices
		//m_coordinate_map.insert(ii, al, INSERTED_POSITION);

		// sanity check
		assert(m_coordinate_map.size() == m_sequence.size());
		return true;

	}
	else
	{
		std::cout << " You are deleting more than the haplotype's length." << std::endl;
		return false;
	}
}




//       
bool Haplotype::apply_variant(const Variant& v)
{
	// Search the coordinate map for the reference position
	size_t derived_idx = _find_derived_index_by_ref_lower_bound(v.ref_position);

	// if we could not find the reference position in the map
	// this variant is incompatable with the haplotype, do nothing
	if(derived_idx == m_coordinate_map.size() ||
	   m_coordinate_map[derived_idx] != v.ref_position)
	{
		return false;
	}
   
	// Check that the string matches
	size_t rl = v.ref_seq.length();
	size_t al = v.alt_seq.length();

	// std::cout << "v.var_seq (haplotype): " << v.var_seq << std::endl;
	// std::cerr << "SV " << v.ref_position << " ALT=" << v.alt_seq << " REF="<< v.ref_seq << std::endl;

	if (v.alt_seq == "<DEL>")
	{
		// Replace reference sequence [derived_idx ... derived_idx+rl] by [derived_idx]
		rl = v.ref_length; // ref_length is set when reading vcf line (see Variant constructor)
		al = 1;

		m_sequence.replace(derived_idx, rl, m_sequence.substr(derived_idx, 1));
		std::cerr << "---------- [DEL] length = " << rl << std::endl;
	}
	else if (v.alt_seq == "<INS>")
	{
		//std::regex regex_seq_info(";SEQ=([A|C|G|T|a|c|g|t]*);");
		//auto regex_it = std::sregex_iterator(v.info.begin(), v.info.end(), regex_seq_info);
		std::string insertion_seq{};

		//if (regex_it != std::sregex_iterator())
		//{ // there is a match (it should be)
			insertion_seq = m_sequence.substr(derived_idx, 1) + v.var_seq;
			//insertion_seq = m_sequence.substr(derived_idx, 1) + ((*regex_it)[1]).str();
		//}
		//else
		//   return false;
		//insertion_seq = v.var_seq;
		rl = 1;
		al = insertion_seq.length();
		std::cerr << "---------- [INS] length and v.alt_length = " << al << "and" << v.alt_length << std::endl;
		assert(al == v.alt_length); // regex should work the same when reading the variant

		m_sequence.replace(derived_idx, rl, insertion_seq);
	}
	else
	{
		// no match, variant conflicts with haplotype sequence
		if (m_sequence.substr(derived_idx, rl) != v.ref_seq) {
			return false;
		}

		// update sequence
		m_sequence.replace(derived_idx, rl, v.alt_seq);
	}

	// update coordinate map

	// make a pair of iterators that bound the changed sequence
	// std::cout << "m_coordinate_map.size(): " << m_coordinate_map.size() << std::endl;
	// std::cout << "derived_idx: " << derived_idx << std::endl;
	// std::cout << "rl: " << rl << std::endl;

	std::vector<size_t>::iterator fi = m_coordinate_map.begin() + derived_idx;
	std::vector<size_t>::iterator li = fi + rl;

	// erase the positions of the changed bases
	//assert((derived_idx+rl) < m_coordinate_map.size());
	std::vector<size_t>::iterator ii = m_coordinate_map.erase(fi, li);

	// insert new positions for the alt bases with invalid indices
	m_coordinate_map.insert(ii, al, INSERTED_POSITION);

	// sanity check
	assert(m_coordinate_map.size() == m_sequence.size());

	m_variants.push_back(v);
	return true;
}

bool Haplotype::apply_variants(const std::vector<Variant>& variants)
{
	bool good = true;
	for(const auto& v : variants) {
		good = good && apply_variant(v);
	}
	return good;
}


// return a new haplotype subsetted by reference coordinates
Haplotype Haplotype::substr_by_reference(size_t start, size_t end) const
{
	assert(start >= m_ref_position);
	assert(start <= m_ref_position + m_reference.length());
	
	assert(end >= m_ref_position);
	assert(end <= m_ref_position + m_reference.length());

	size_t derived_base_start = _find_derived_index_by_ref_lower_bound(start);
	size_t derived_base_end = _find_derived_index_by_ref_lower_bound(end);
	
	// Bump out the reference coordinate to encompass the complete range (start, end)
	while(m_coordinate_map[derived_base_start] > start ||
		  m_coordinate_map[derived_base_start] == INSERTED_POSITION)
	{ 
		derived_base_start -= 1;
	}

	assert(derived_base_start != m_coordinate_map.size());

	// JTS: temporary debug dump for #199
	if(derived_base_end == m_coordinate_map.size()) {
		print_debug_info();
	}

	assert(derived_base_end != m_coordinate_map.size());
	assert(m_coordinate_map[derived_base_start] <= start);
	assert(m_coordinate_map[derived_base_end] >= end);

	start = m_coordinate_map[derived_base_start];
	end = m_coordinate_map[derived_base_end];
	
	Haplotype ret(m_ref_name,
				  start,
				  m_reference.substr(start - m_ref_position, end - start + 1));

	ret.m_sequence = m_sequence.substr(derived_base_start, derived_base_end - derived_base_start + 1);
	ret.m_coordinate_map = std::vector<size_t>(m_coordinate_map.begin() + derived_base_start,
											   m_coordinate_map.begin() + derived_base_end + 1);

	assert(ret.m_coordinate_map.front() == start);
	assert(ret.m_coordinate_map.back() == end);
	assert(ret.m_coordinate_map.size() == ret.m_sequence.size());

	return ret;
}

size_t Haplotype::get_reference_position_for_haplotype_base(size_t i) const
{
	assert(i < m_coordinate_map.size());
	return m_coordinate_map[i] == INSERTED_POSITION ? std::string::npos : m_coordinate_map[i];
}

void Haplotype::get_enclosing_reference_range_for_haplotype_range(size_t& hap_lower, size_t& hap_upper,
																  size_t& ref_lower, size_t& ref_upper) const
{
	while(hap_lower > 0 && m_coordinate_map[hap_lower] == INSERTED_POSITION) {
		hap_lower--;
	}

	while(hap_upper < m_coordinate_map.size() && m_coordinate_map[hap_upper] == INSERTED_POSITION) {
		hap_upper++;
	}

	if(hap_lower == 0 || hap_upper >= m_coordinate_map.size()) {
		hap_lower = hap_upper = ref_lower = ref_upper = std::string::npos;
	} else {
		ref_lower = m_coordinate_map[hap_lower];
		ref_upper = m_coordinate_map[hap_upper];
	}
}

size_t Haplotype::_find_derived_index_by_ref_lower_bound(size_t ref_index) const
{
	for(size_t i = 0; i < m_coordinate_map.size(); ++i) {
		if(m_coordinate_map[i] != INSERTED_POSITION && m_coordinate_map[i] >= ref_index) {
			return i;
		}
	}
	return m_coordinate_map.size();
}

void Haplotype::print_debug_info() const
{
	fprintf(stderr, "[haplotype-debug] ctg: %s, position: %lu\n", m_ref_name.c_str(), m_ref_position);
	fprintf(stderr, "[haplotype-debug] r-sequence: %s\n", m_reference.c_str());
	fprintf(stderr, "[haplotype-debug] h-sequence: %s\n", m_sequence.c_str());
	for(size_t i = 0; i < m_variants.size(); ++i) {
		fprintf(stderr, "[haplotype-debug] variant[%zu]: ", i); 
		m_variants[i].write_vcf(stderr);
	}
}
