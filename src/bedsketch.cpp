#include "bedsketch.h"
#include <sys/stat.h>

namespace dashing2 {

// Every option that changes a BED sketch is part of its cache file name, followed by the
// sketch suffix (.ss, .opss, .bmh or .pmh), so a cache is only reused under the same options.
static std::string bed_cache_suffix(const Dashing2Options &opts) {
    std::string ret = ".sketchsize" + std::to_string(opts.sketchsize_);
    if(opts.count_threshold_ > 0) {
        ret += ".ct_threshold";
        if(std::fmod(opts.count_threshold_, 1.)) ret += std::to_string(opts.count_threshold_);
        else ret += std::to_string(int(opts.count_threshold_));
    }
    if(!opts.trim_chr_)
        ret += ".notrimchr";
    if(opts.bed_parse_normalize_intervals_)
        ret += ".normalize_intervals";
    if(opts.sspace_ != SPACE_SET) {
        ret += '.';
        ret += to_string(opts.ct());
        if(opts.ct() != EXACT_COUNTING)
            ret += std::to_string(opts.cssize_);
    }
    return ret + to_suffix(opts);
}

// Modification time in nanoseconds, or -1 if the file cannot be stat'ed
static int64_t bed_mtime_ns(const std::string &path) {
    struct stat st;
    if(::stat(path.data(), &st)) return -1;
#ifdef __APPLE__
    const auto &ts = st.st_mtimespec;
#else
    const auto &ts = st.st_mtim;
#endif
    return int64_t(ts.tv_sec) * 1000000000 + ts.tv_nsec;
}

std::pair<std::vector<RegT>, double> bed2sketch(const std::string &path, const Dashing2Options &opts) {
    if(opts.sspace_ > SPACE_PSET) throw std::invalid_argument("Can't do edit distance for BED files");
    if(opts.bed_parse_normalize_intervals_ && opts.sspace_ == SPACE_SET)
        throw std::invalid_argument("Can't normalize BED rows in set space. Use SPACE_MULTISET or SPACE_PSET");
    std::ifstream ifs(path);
    const bool op = opts.one_perm();
    FullSetSketch ss(opts.count_threshold_, opts.sketchsize_);
    OPSetSketch opss(opts.sketchsize_);
    Counter ctr(opts.cssize_);
    std::pair<std::vector<RegT>, double> ret({std::vector<RegT>(opts.sketchsize_), 0.});
    auto &retvec(ret.first);
    const std::string suffix = bed_cache_suffix(opts);
    std::string cache_path = path + suffix;
    DBG_ONLY(std::fprintf(stderr, "Using %s\n", op ? "oneperm": "fullsetsketch");)

    if(opts.trim_folder_paths()) {
        // The cache file keeps the sketch suffix, so it is never named like an input BED file.
        // Files with the same name in different directories must not share a cache file,
        // so the file name is tagged with a hash of its absolute path.
        char *const abspath = ::realpath(path.data(), nullptr);
        const std::string key = abspath ? abspath: path;
        std::free(abspath);
        char buf[24];
        cache_path = trim_folder(path) + std::string(buf, std::snprintf(buf, sizeof(buf), ".%016llx", static_cast<unsigned long long>(XXH3_64bits(key.data(), key.size())))) + suffix;
        if(opts.outprefix_.size())
            cache_path = opts.outprefix_ + '/' + cache_path;
    }
    // A cache file older than its input is not reused, so an input that was overwritten is sketched again.
    if(opts.cache_sketches_ && bns::isfile(cache_path) && bed_mtime_ns(path) <= bed_mtime_ns(cache_path)) {
        auto [ifp, ispopen] = xopen(cache_path);
        // The cache holds the cardinality followed by exactly sketchsize_ registers.
        const bool ok = std::fread(&ret.second, sizeof(ret.second), 1, ifp) == 1
                     && std::fread(retvec.data(), sizeof(RegT), retvec.size(), ifp) == retvec.size();
        if(ispopen) ::pclose(ifp); else std::fclose(ifp);
        if(!ok) throw std::runtime_error("Failed to read cached BED sketch from " + cache_path);
        return ret;
    }
    for(std::string s;std::getline(ifs, s);) {
        if(s.empty() || s.front() == '#') continue;
        char *p = s.data(), *p2;
        if((p2 = std::strchr(p, '\t')) == nullptr)
            throw std::invalid_argument(std::string("Malformed line: ") + s);
        if(opts.trim_chr_ && ((*p == 'c' || *p == 'C') && p[1] == 'h' && p[2] == 'r'))
            p += 3;
        const uint64_t chrhash = XXH3_64bits(p, p2 - p);
        p = p2 + 1;
        const unsigned long start = std::strtoul(p, &p, 10), stop =  std::strtoul(p, &p, 10);
        const double inc = opts.bed_parse_normalize_intervals_ ? 1. / (stop - start): 1;
        // Consider SIMDifying packing these before adding?
        // If set space, sketch directly.
        if(opts.sspace_ == SPACE_SET) {
            if(op) for(auto i = start; i < stop; opss.update(chrhash ^ i++));
            else   for(auto i = start; i < stop; ss.update(chrhash ^ i++));
        } else {
            // else, we need to compute counts before we sketch
            for(auto i = start; i < stop; ctr.add(chrhash ^ i++, inc));
        }
    }
    if(opts.sspace_ > SPACE_SET) {
        if(opts.ct() == EXACT_COUNTING) {
            if(opts.sspace_ == SPACE_MULTISET) {
                BagMinHash bmh(opts.sketchsize_);
                ctr.finalize(bmh);
                std::copy(bmh.data(), bmh.data() + opts.sketchsize_, retvec.data());
                ret.second = bmh.total_weight();
            } else {
                sketch::pmh2_t pmh(opts.sketchsize_);
                ctr.finalize(pmh);
                std::copy(pmh.data(), pmh.data() + opts.sketchsize_, retvec.data());
                ret.second = pmh.total_weight();
            }
        } else {
#define __FS() do {\
    for(size_t i = 0; i < csz; ++i) sketcher.update(i, ctr.count_sketch_[i]);\
    auto p = sketcher.data();\
    ret.second = sketcher.total_weight();\
    std::copy(p, p + opts.sketchsize_, retvec.data());\
    } while(0)
            const size_t csz = ctr.count_sketch_.size();
            if(opts.sspace_ == SPACE_MULTISET) {
                BagMinHash sketcher(opts.sketchsize_);
                __FS();
            } else {
                sketch::pmh2_t sketcher(opts.sketchsize_);
                __FS();
            }
#undef __FS
        }
    } else {
        ret.second = op ? opss.getcard(): ss.getcard();
        RegT *sptr = op ? opss.data(): ss.data();
        std::copy(sptr, sptr + opts.sketchsize_, retvec.data());
    }
    // Like sketches of sequence files, BED sketches are written to disk only with --cache.
    if(opts.cache_sketches_) {
        std::FILE *ofp = bfopen(cache_path.data(), "wb");
        if(ofp == nullptr) THROW_EXCEPTION(std::runtime_error(std::string("Failed to open file ") + cache_path + " for writing sketch."));
        std::fwrite(&ret.second, 1, sizeof(ret.second), ofp);
        std::fwrite(retvec.data(), opts.sketchsize_, sizeof(RegT), ofp);
        std::fclose(ofp);
    }
    return ret;
}

} // namespace dashing2
