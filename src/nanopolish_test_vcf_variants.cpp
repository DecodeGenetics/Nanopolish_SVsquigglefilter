//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_call_variants -- find variants wrt a reference
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <inttypes.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>
#include <algorithm>
#include <queue>
#include <sstream>
#include <fstream>
#include <set>
#include <omp.h>
#include <getopt.h>
#include <iterator>
#include "htslib/faidx.h"
#include "nanopolish_poremodel.h"
#include "nanopolish_transition_parameters.h"
#include "nanopolish_matrix.h"
#include "nanopolish_klcs.h"
#include "nanopolish_profile_hmm.h"
#include "nanopolish_alignment_db.h"
#include "nanopolish_anchor.h"
#include "nanopolish_variant.h"
#include "nanopolish_haplotype.h"
#include "nanopolish_pore_model_set.h"
#include "nanopolish_duration_model.h"
#include "nanopolish_variant_db.h"
#include "profiler.h"
#include "progress.h"
#include "stdaln.h"
#include "likelihoods_core.h"

#include <random>
#include <cstdlib>
#include <climits>
#include <limits>       // std::numeric_limits


// Macros
#define max3(x,y,z) std::max(std::max(x,y), z)

// Flags to turn on/off debugging information

//#define DEBUG_HMM_UPDATE 1
//#define DEBUG_HMM_EMISSION 1
//#define DEBUG_TRANSITION 1
//#define DEBUG_PATH_SELECTION 1
//#define DEBUG_SINGLE_SEGMENT 1
//#define DEBUG_SHOW_TOP_TWO 1
//#define DEBUG_SEGMENT_ID 193
//#define DEBUG_BENCHMARK 1

// Hack hack hack
float g_p_skip3, g_p_skip_self3, g_p_bad3, g_p_bad_self3;

//
// Getopt
//
#define SUBPROGRAM "variants"

static const char *CONSENSUS_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2015 Ontario Institute for Cancer Research\n";

static const char *CONSENSUS_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] --reads reads.fa --bam alignments.bam --genome genome.fa\n"
"Find SNPs using a signal-level HMM\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --version                        display version\n"
"      --help                           display this help and exit\n"
"      --snps                           only call SNPs\n"
"      --consensus=FILE                 run in consensus calling mode and write polished sequence to FILE\n"
"      --fix-homopolymers               run the experimental homopolymer caller\n"
"      --faster                         minimize compute time while slightly reducing consensus accuracy\n"
"  -w, --window=STR                     find variants in window STR (format: <chromsome_name>:<start>-<end>)\n"
"  -r, --reads=FILE                     the 2D ONT reads are in fasta FILE\n"
"  -b, --bam=FILE                       the reads aligned to the reference genome are in bam FILE\n"
"  -e, --event-bam=FILE                 the events aligned to the reference genome are in bam FILE\n"
"  -g, --genome=FILE                    the reference genome is in FILE\n"
"  -p, --ploidy=NUM                     the ploidy level of the sequenced genome\n"
"  -q  --methylation-aware=STR          turn on methylation aware polishing and test motifs given in STR (example: -m dcm,dam)\n"
"      --genotype=FILE                  call genotypes for the variants in the vcf FILE\n"
"  -o, --outfile=FILE                   write result to FILE [default: stdout]\n"
"  -t, --threads=NUM                    use NUM threads (default: 1)\n"
"  -m, --min-candidate-frequency=F      extract candidate variants from the aligned reads when the variant frequency is at least F (default 0.2)\n"
"  -d, --min-candidate-depth=D          extract candidate variants from the aligned reads when the depth is at least D (default: 20)\n"
"  -x, --max-haplotypes=N               consider at most N haplotype combinations (default: 1000)\n"
"      --max-rounds=N                   perform N rounds of consensus sequence improvement (default: 50)\n"
"  -c, --candidates=VCF                 read variant candidates from VCF, rather than discovering them from aligned reads\n"
"  -a, --alternative-basecalls-bam=FILE if an alternative basecaller was used that does not output event annotations\n"
"                                       then use basecalled sequences from FILE. The signal-level events will still be taken from the -b bam.\n"
"      --calculate-all-support          when making a call, also calculate the support of the 3 other possible bases\n"
"      --models-fofn=FILE               read alternative k-mer models from FILE\n"
"  -f  --flanking-sequence-len=INT      number of bases to look before and after the SV.\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
	static unsigned int verbose = 0;
	static std::string reads_file;
	static std::string bam_file;
	static std::string event_bam_file;
	static std::string genome_file;
	static std::string output_file;
	static std::string candidates_file;
	static std::string models_fofn;
	static std::string window;
	static std::string consensus_output;
	static std::string alternative_model_type = DEFAULT_MODEL_TYPE;
	static std::string alternative_basecalls_bam;
	static double min_candidate_frequency = 0.2f;
	static int min_candidate_depth = 20;
	static int calculate_all_support = false;
	static int snps_only = 0;
	static int show_progress = 0;
	static int num_threads = 1;
	static int consensus_mode = 0;
	static int fix_homopolymers = 0;
	static int genotype_only = 0;
	static int ploidy = 0;
	static int min_distance_between_variants = 10;
	static int min_flanking_sequence = 30;
	static int max_haplotypes = 1000;
	static int max_rounds = 50;
	static int screen_score_threshold = 100;
	static int screen_flanking_sequence = 50;
	static int debug_alignments = 0;
	static std::vector<std::string> methylation_types;
}

static const char* shortopts = "r:b:g:t:w:o:e:m:f:c:d:a:x:q:v";

enum { OPT_HELP = 1,
	   OPT_VERSION,
	   OPT_VCF,
	   OPT_PROGRESS,
	   OPT_SNPS_ONLY,
	   OPT_CALC_ALL_SUPPORT,
	   OPT_CONSENSUS,
	   OPT_FIX_HOMOPOLYMERS,
	   OPT_GENOTYPE,
	   OPT_MODELS_FOFN,
	   OPT_MAX_ROUNDS,
	   OPT_EFFORT,
	   OPT_FASTER,
	   OPT_P_SKIP,
	   OPT_P_SKIP_SELF,
	   OPT_P_BAD,
	   OPT_P_BAD_SELF };

static const struct option longopts[] = {
	{ "verbose",                   no_argument,       NULL, 'v' },
	{ "reads",                     required_argument, NULL, 'r' },
	{ "bam",                       required_argument, NULL, 'b' },
	{ "event-bam",                 required_argument, NULL, 'e' },
	{ "genome",                    required_argument, NULL, 'g' },
	{ "window",                    required_argument, NULL, 'w' },
	{ "outfile",                   required_argument, NULL, 'o' },
	{ "threads",                   required_argument, NULL, 't' },
	{ "min-candidate-frequency",   required_argument, NULL, 'm' },
	{ "min-candidate-depth",       required_argument, NULL, 'd' },
	{ "max-haplotypes",            required_argument, NULL, 'x' },
	{ "candidates",                required_argument, NULL, 'c' },
	{ "ploidy",                    required_argument, NULL, 'p' },
	{ "alternative-basecalls-bam", required_argument, NULL, 'a' },
	{ "methylation-aware",         required_argument, NULL, 'q' },
	{ "flanking-sequence-len",     required_argument, NULL, 'f' },
	{ "effort",                    required_argument, NULL, OPT_EFFORT },
	{ "max-rounds",                required_argument, NULL, OPT_MAX_ROUNDS },
	{ "genotype",                  required_argument, NULL, OPT_GENOTYPE },
	{ "models-fofn",               required_argument, NULL, OPT_MODELS_FOFN },
	{ "p-skip",                    required_argument, NULL, OPT_P_SKIP },
	{ "p-skip-self",               required_argument, NULL, OPT_P_SKIP_SELF },
	{ "p-bad",                     required_argument, NULL, OPT_P_BAD },
	{ "p-bad-self",                required_argument, NULL, OPT_P_BAD_SELF },
	{ "consensus",                 required_argument, NULL, OPT_CONSENSUS },
	{ "faster",                    no_argument,       NULL, OPT_FASTER },
	{ "fix-homopolymers",          no_argument,       NULL, OPT_FIX_HOMOPOLYMERS },
	{ "calculate-all-support",     no_argument,       NULL, OPT_CALC_ALL_SUPPORT },
	{ "snps",                      no_argument,       NULL, OPT_SNPS_ONLY },
	{ "progress",                  no_argument,       NULL, OPT_PROGRESS },
	{ "help",                      no_argument,       NULL, OPT_HELP },
	{ "version",                   no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

// If there is a single contig in the .fai file, return its name
// otherwise print an error message and exit
std::string get_single_contig_or_fail3()
{
	faidx_t *fai = fai_load(opt::genome_file.c_str());
	size_t n_contigs = faidx_nseq(fai);
	if(n_contigs > 1) {
		fprintf(stderr, "Error: genome has multiple contigs, please use -w to specify input region\n");
		exit(EXIT_FAILURE);
	}

	const char* name = faidx_iseq(fai, 0);
	return std::string(name);
}

int get_contig_length3(const std::string& contig)
{ 
	faidx_t *fai = fai_load(opt::genome_file.c_str());
	int len = faidx_seq_len(fai, contig.c_str());
	if(len == -1) {
		fprintf(stderr, "error: faidx could not get contig length for contig %s\n", contig.c_str());
		exit(EXIT_FAILURE);
	}
	fai_destroy(fai);
	return len;
}

void gather_reads_spanning_variant(FILE* out_fp)
{

	fprintf(out_fp, "#chr\tSV_index\tflank\tpos\tbreakpoint\tSV_type\tvariant_length\tread_id\tread_length\tevent_length\tmap_qual\tflanksize\tbase_score\talt_score\tdiff\n");

	uint32_t alignment_flags = HAF_ALLOW_PRE_CLIP | HAF_ALLOW_POST_CLIP;
	std::vector<Variant> all_variants = read_variants_from_file(opt::candidates_file, true);
	for (int vv = 0; vv < all_variants.size(); vv++)
	{
		Variant test_variant = all_variants[vv];
		
		std::string contig = test_variant.ref_name;
		
		int contig_length = get_contig_length3(contig);
		unsigned int start_base = 0;
		unsigned int end_base = contig_length - 1;

		int variant_start = test_variant.ref_position;
		int variant_end = test_variant.ref_position + test_variant.ref_length;
		int variant_len = std::max(test_variant.ref_length, test_variant.alt_length);

		const int BUFFER = 500;
 
		int start_a = test_variant.ref_position - opt::screen_flanking_sequence - BUFFER;
		int start_pos = std::max(start_a, 0);

		int stop_a = test_variant.ref_position + test_variant.ref_length + opt::screen_flanking_sequence + BUFFER;
		int stop_pos  = std::min(stop_a, contig_length);
		
		AlignmentDB alignments(opt::reads_file, opt::genome_file, opt::bam_file, opt::event_bam_file);
		if(!opt::alternative_basecalls_bam.empty()) {
			alignments.set_alternative_basecalls_bam(opt::alternative_basecalls_bam);
		}

		// Separate the begin and end sites of the SV, if they're larger than this bp.
		const int SV_LEN_CUTOFF_TO_SEPARATE = 50;
		const int TOO_SHORT_SV_LEN = 30;
		const int TOO_LONG_SV_LEN = 1000000;
		if ((variant_len < TOO_SHORT_SV_LEN) || (variant_len > TOO_LONG_SV_LEN))
		{
			// test the usual way, perhaps choose the begin/end coords better.
			std::cout << variant_len << ": variant_len" << std::endl;
			std::cout << test_variant.ref_position << ": test_variant.ref_position" << std::endl;
			std::cout << "SV is too short or too long." << std::endl;
		}
		else
		{
			if (test_variant.alt_seq == "<DEL>" || test_variant.alt_seq == "<INS>")
			{	
				// Doing LEFT FLANK.
				std::cout << "Left Flank Analysis " << std::endl;
				int min_required_ref_coord;
				auto seq_records = get_min_leftflankrefpos_and_seqrecords(alignments, contig, variant_start, SV_LEN_CUTOFF_TO_SEPARATE, opt::bam_file, min_required_ref_coord);

				// min required ref coord is found. No need for max on left flank because that is hand given.
				alignments.load_region_select_reads(contig, std::max(min_required_ref_coord - BUFFER, 0), stop_pos, test_variant.supporting_reads);
				left_flank_analysis(alignments, alignment_flags, contig, contig_length, seq_records, test_variant, opt::screen_flanking_sequence, opt::methylation_types, out_fp);

				// Doing RIGHT FLANK
				std::cout << "Right Flank Analysis " << std::endl;
				int max_required_ref_coord;
				seq_records = get_max_rightflankrefpos_and_seqrecords(alignments, contig, variant_start, test_variant.ref_length, SV_LEN_CUTOFF_TO_SEPARATE, opt::bam_file, max_required_ref_coord);

				// min required ref coord is found. No need for max on left flank because that is hand given.
				alignments.load_region_select_reads(contig, start_pos, std::min(max_required_ref_coord + BUFFER, contig_length), test_variant.supporting_reads);
				right_flank_analysis(alignments, alignment_flags, contig, contig_length, seq_records, test_variant, opt::screen_flanking_sequence, opt::methylation_types, out_fp);
			}
			else
			{
				std::cout << "SVs other than DEL or INS not supported yet." << std::endl;
			}
		}
	}

}

   
void parse_test_vcf_variants_options(int argc, char** argv)
{
	std::string methylation_motifs_str;
	bool die = false;
	for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
		std::istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case 'r': arg >> opt::reads_file; break;
			case 'g': arg >> opt::genome_file; break;
			case 'b': arg >> opt::bam_file; break;
			case 'e': arg >> opt::event_bam_file; break;
			case 'w': arg >> opt::window; break;
			case 'o': arg >> opt::output_file; break;
			case 'm': arg >> opt::min_candidate_frequency; break;
			case 'f': arg >> opt::screen_flanking_sequence; break;
			case 'd': arg >> opt::min_candidate_depth; break;
			case 'x': arg >> opt::max_haplotypes; break;
			case 'c': arg >> opt::candidates_file; break;
			case 'p': arg >> opt::ploidy; break;
			case 'q': arg >> methylation_motifs_str; break;
			case 'a': arg >> opt::alternative_basecalls_bam; break;
			case '?': die = true; break;
			case 't': arg >> opt::num_threads; break;
			case 'v': opt::verbose++; break;
			case OPT_CONSENSUS: arg >> opt::consensus_output; opt::consensus_mode = 1; break;
			case OPT_FIX_HOMOPOLYMERS: opt::fix_homopolymers = 1; break;
			case OPT_EFFORT: arg >> opt::screen_score_threshold; break;
			case OPT_FASTER: opt::screen_score_threshold = 25; break;
			case OPT_MAX_ROUNDS: arg >> opt::max_rounds; break;
			case OPT_GENOTYPE: opt::genotype_only = 1; arg >> opt::candidates_file; break;
			case OPT_MODELS_FOFN: arg >> opt::models_fofn; break;
			case OPT_CALC_ALL_SUPPORT: opt::calculate_all_support = 1; break;
			case OPT_SNPS_ONLY: opt::snps_only = 1; break;
			case OPT_PROGRESS: opt::show_progress = 1; break;
			case OPT_P_SKIP: arg >> g_p_skip3; break;
			case OPT_P_SKIP_SELF: arg >> g_p_skip_self3; break;
			case OPT_P_BAD: arg >> g_p_bad3; break;
			case OPT_P_BAD_SELF: arg >> g_p_bad_self3; break;
			case OPT_HELP:
				std::cout << CONSENSUS_USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				std::cout << CONSENSUS_VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
		}
	}

	if (argc - optind < 0) {
		std::cerr << SUBPROGRAM ": missing arguments\n";
		die = true;
	} else if (argc - optind > 0) {
		std::cerr << SUBPROGRAM ": too many arguments\n";
		die = true;
	}

	if(opt::num_threads <= 0) {
		std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::num_threads << "\n";
		die = true;
	}

	if(opt::reads_file.empty()) {
		std::cerr << SUBPROGRAM ": a --reads file must be provided\n";
		die = true;
	}

	if(opt::genome_file.empty()) {
		std::cerr << SUBPROGRAM ": a --genome file must be provided\n";
		die = true;
	}

	if(!opt::consensus_mode && opt::ploidy == 0) {
		std::cerr << SUBPROGRAM ": --ploidy parameter must be provided\n";
		die = true;
	} else if(opt::consensus_mode) {
		opt::ploidy = 1;
	}

	if(opt::bam_file.empty()) {
		std::cerr << SUBPROGRAM ": a --bam file must be provided\n";
		die = true;
	}

	if(!opt::models_fofn.empty()) {
		// initialize the model set from the fofn
		PoreModelSet::initialize(opt::models_fofn);
	}

	if(!methylation_motifs_str.empty()) {
		opt::methylation_types = split(methylation_motifs_str, ',');
		for(const std::string& mtype : opt::methylation_types) {
			// this call will abort if the alphabet does not exist
			const Alphabet* alphabet = get_alphabet_by_name(mtype);
			assert(alphabet != NULL);
		}
	}

	if (die)
	{
		std::cout << "\n" << CONSENSUS_USAGE_MESSAGE;
		exit(EXIT_FAILURE);
	}
}

int test_vcf_variants_main(int argc, char** argv)
{
	parse_test_vcf_variants_options(argc, argv);
	omp_set_num_threads(opt::num_threads);

	std::string contig;
	int start_base;
	int end_base;
	int contig_length = -1;

	//If a window has been specified, only call variants/polish in that range
	if(!opt::window.empty()) {
		// Parse the window string
		std::cout << "window is not empty" << std::endl;
		parse_region_string(opt::window, contig, start_base, end_base);
		contig_length = get_contig_length3(contig);
		if (end_base == 0) // not altered by parse region
			end_base = contig_length - 1;
		end_base = std::min(end_base, contig_length - 1);
		std::cout << "end is: " << end_base << std::endl;
		std::cout << "start is: " << start_base << std::endl;
	} else {
		// otherwise, run on the whole genome
		std::cout << "window is empty" << std::endl;

		contig = get_single_contig_or_fail3();
		contig_length = get_contig_length3(contig);
		start_base = 0;
		end_base = contig_length - 1;
	}

	int MIN_DISTANCE_TO_END = 40;
	if(contig_length - start_base < MIN_DISTANCE_TO_END) {
		fprintf(stderr, "Invalid polishing window: [%d %d] - please adjust -w parameter.\n", start_base, end_base);
		fprintf(stderr, "The starting coordinate of the polishing window must be at least %dbp from the contig end\n", MIN_DISTANCE_TO_END);
		exit(EXIT_FAILURE);
	}

	FILE* out_fp;
	if(!opt::output_file.empty()) {
		out_fp = fopen(opt::output_file.c_str(), "w");
	} else {
		out_fp = stdout;
	}

	// Build the VCF header
	std::vector<std::string> tag_fields;

	//
	tag_fields.push_back(
		Variant::make_vcf_tag_string("INFO", "TotalReads", 1, "Integer",
									  "The number of event-space reads used to call the variant"));

	tag_fields.push_back(
		Variant::make_vcf_tag_string("INFO", "SupportFraction", 1, "Float",
									  "The fraction of event-space reads that support the variant"));

	tag_fields.push_back(
		Variant::make_vcf_tag_string("INFO", "BaseCalledReadsWithVariant", 1, "Integer",
									  "The number of base-space reads that support the variant"));

	tag_fields.push_back(
		Variant::make_vcf_tag_string("INFO", "BaseCalledFraction", 1, "Float",
									  "The fraction of base-space reads that support the variant"));

	tag_fields.push_back(
			Variant::make_vcf_tag_string("INFO", "AlleleCount", 1, "Integer",
				"The inferred number of copies of the allele"));

	if(opt::calculate_all_support) {
		tag_fields.push_back(
				Variant::make_vcf_tag_string("INFO", "SupportFractionByBase", 4, "Integer",
					"The fraction of reads supporting A,C,G,T at this position"));

	}
	tag_fields.push_back(
			Variant::make_vcf_tag_string("FORMAT", "GT", 1, "String",
				"Genotype"));

	gather_reads_spanning_variant(out_fp);


	return 0;
}
