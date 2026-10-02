#pragma once
#ifndef DASHING2_BEDSKETCH_H__
#define DASHING2_BEDSKETCH_H__
#include "d2.h"

namespace dashing2 {
std::pair<std::vector<RegT>, double> bed2sketch(const std::string &path, const Dashing2Options &opts);
// Cache file name of a BED or BigWig input: the input path followed by suffix, which names the
// options that change its sketch; under --outprefix, the file name tagged with a hash of its absolute path.
std::string interval_cache_path(const std::string &path, const Dashing2Options &opts, const std::string &suffix);
// Whether cache_path exists and path was not modified after it was written
bool interval_cache_is_fresh(const std::string &cache_path, const std::string &path);
}

#endif
