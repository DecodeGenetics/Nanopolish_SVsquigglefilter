//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_common -- Data structures and definitions
// shared across files
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nanopolish_common.h"
#include "nanopolish_squiggle_read.h"

// Split a string into parts based on the delimiter
std::vector<std::string> split(std::string in, char delimiter)
{
    std::vector<std::string> out;
    size_t lastPos = 0;
    size_t pos = in.find_first_of(delimiter);

    while(pos != std::string::npos)
    {
        out.push_back(in.substr(lastPos, pos - lastPos));
        lastPos = pos + 1;
        pos = in.find_first_of(delimiter, lastPos);
    }
    out.push_back(in.substr(lastPos));
    return out;
}

void parse_region_string(const std::string& region, std::string& contig, int& start, int& end)
{
    std::string region_copy = region;

    // Parse the window string

    // Delete commas
    region_copy.erase(std::remove(region_copy.begin(), region_copy.end(), ','), region_copy.end());
    
    // Replace ":" and "-" with spaces to make it parseable with stringstream
    std::replace(region_copy.begin(), region_copy.end(), ':', ' ');
    std::replace(region_copy.begin(), region_copy.end(), '-', ' ');

    std::stringstream parser(region_copy);
    
    //modified by dorukb such that we could use the contig name only.
    // TODO must be checked for compatibility with original nanopolish
    if (region_copy.find(' ') == std::string::npos) // only contig name
        parser >> contig;
    
    start = 0;
    end=0;
    parser >> contig >> start >> end;

    std::cout << "parser start: " << start << std::endl;
    std::cout << "parser end: " << end << std::endl;
}

SemVer parse_semver_string(const std::string& semver_str)
{
    SemVer out;
    const auto& vec = split(semver_str, '.');
    if(vec.size() == 3) {
        out = { atoi(vec[0].c_str()),
                atoi(vec[1].c_str()),
                atoi(vec[2].c_str()) 
              };
    } else {
        out = { 0, 0, 0 };
    }
    return out;
}

bool ends_with(const std::string& str, const std::string& suffix)
{
    if(suffix.empty()) {
        return true;
    }

    size_t pos = str.find(suffix);
    if(pos == std::string::npos) {
        return false;
    }
    return pos + suffix.size() == str.length();
}
// from: http://stackoverflow.com/questions/9330915/number-of-combinations-n-choose-r-in-c
size_t nChoosek(size_t n, size_t k)
{
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    int result = n;
    for( int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}

void bam_index_error_exit(const std::string& bam_filename)
{
    fprintf(stderr, "Error: could not load the .bai index file for %s\n", bam_filename.c_str());
    fprintf(stderr, "Please run 'samtools index %s' before nanopolish\n", bam_filename.c_str());
    exit(EXIT_FAILURE);
}
