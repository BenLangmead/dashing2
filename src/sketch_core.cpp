#include "sketch_core.h"
#include "cmp_main.h"
#include <cinttypes>

namespace dashing2 {
extern size_t MEMSIGTHRESH;

// .xz, .bz2 and .zst inputs are read through a decompression pipe, which starts even when the
// file is missing, so every sequence file is checked before sketching.
static void require_readable_inputs(const std::vector<std::string> &paths) {
    for(const auto &line: paths) {
        for_each_substr([](const std::string &path) {
            if(std::FILE *fp = std::fopen(path.data(), "rb")) {
                std::fclose(fp);
                return;
            }
            std::fprintf(stderr, "Could not open input file %s: %s\n", path.data(), std::strerror(errno));
            std::exit(EXIT_FAILURE);
        }, line);
    }
}

INLINE size_t nbytes_from_line(const std::string &line) {
    size_t ret = 0;
    for_each_substr([&ret](const std::string &s) {ret += bns::filesize(s.data());}, line);
    return ret;
}

// Sketch and k-mer files written by dashing2 are binary; parsing one as FASTA/FASTQ yields meaningless k-mers.
static void reject_sketch_files(const std::vector<std::string> &paths) {
    static constexpr const char *sketch_exts[] = {".ss", ".opss", ".bmh", ".pmh", ".kmerset64", ".kmerset128", ".mmerseq64", ".mmerseq128"};
    for(const auto &line: paths) {
        for_each_substr([](const std::string &path) {
            if(std::none_of(std::begin(sketch_exts), std::end(sketch_exts), [&](const char *ext) {return endswith(path, ext);}))
                return;
            std::FILE *fp = std::fopen(path.data(), "rb");
            if(!fp) return;
            const int c = std::fgetc(fp);
            std::fclose(fp);
            if(c != '>' && c != '@') {
                std::fprintf(stderr, "%s looks like a dashing2 sketch file rather than FASTA/FASTQ; use cmp --presketched to compare sketch files\n", path.data());
                std::exit(EXIT_FAILURE);
            }
        }, line);
    }
}

SketchingResult &sketch_core(SketchingResult &result, Dashing2DistOptions &opts, const std::vector<std::string> &paths, std::string &outfile) {
    if(opts.kmer_result() == FULL_MMER_SEQUENCE && outfile.empty()) {
        THROW_EXCEPTION(std::runtime_error("outfile must be specified for --seq mode."));
    }
    if(opts.dtype_ == DataType::FASTX) require_readable_inputs(paths);
    result.signatures_.memthreshold(MEMSIGTHRESH);
    result.kmers_.memthreshold(MEMSIGTHRESH);
    const size_t npaths = paths.size();
    std::string tmpfile;
    if(opts.dtype_ == DataType::FASTX) {
        reject_sketch_files(paths);
        if(opts.parse_by_seq_) {
            if(paths.size() != 1) {
                result.nperfile_.resize(paths.size());
                THROW_EXCEPTION(std::runtime_error("parse-by-seq currently only handles one file at a time. To process multiple files, simply concatenate them into one file, and run dashing2 on that."));
            }
            KSeqHolder kseqs(std::max(opts.nthreads(), 1u));
            fastx2sketch_byseq(result, opts, paths.front(), kseqs.kseqs_, outfile, true, 512);
        } else {
            fastx2sketch(result, opts, paths, outfile);
        }
    } else if(opts.dtype_ == DataType::LEAFCUTTER) {
        if(outfile.empty()) THROW_EXCEPTION(std::runtime_error("Outfile required for LeafCutter input."));
        auto res = lf2sketch(paths, opts);
        result.names_ = std::move(res.sample_names());
        result.nperfile_.resize(res.nsamples_per_file().size());
        std::copy(res.nsamples_per_file().begin(), res.nsamples_per_file().end(), result.nperfile_.begin());
        std::FILE *of = bfopen(outfile.data(), "wb");
        uint64_t ns = result.names_.size();
        checked_fwrite(&ns, sizeof(ns), 1, of);
        ns = opts.sketchsize_;
        checked_fwrite(&ns, sizeof(ns), 1, of);
        std::fclose(of);
        result.signatures_.assign(outfile);
        result.signatures_.resize(res.registers().size());
        std::copy(res.registers().begin(), res.registers().end(), result.signatures_.begin());
    } else if(opts.dtype_ == DataType::BED || opts.dtype_ == DataType::BIGWIG) {
        // Interval inputs are only sketched; there are no k-mer set or sequence files to compare exactly
        if(opts.kmer_result_ >= FULL_MMER_SET) {
            std::fprintf(stderr, "--bed and --bigwig inputs support sketches only, not --set, --countdict or --seq\n");
            std::exit(EXIT_FAILURE);
        }
        std::vector<std::pair<size_t, uint64_t>> filesizes = get_filesizes(paths);
        result.signatures_.resize(npaths * opts.sketchsize_);
        result.names_.resize(npaths);
        result.cardinalities_.resize(npaths);
        if(opts.dtype_ == DataType::BED) {
            OMP_PFOR_DYN
            for(size_t i = 0; i < npaths; ++i) {
                auto myind = filesizes.size() ? filesizes[i].second: uint64_t(i);
                auto &p(paths[myind]);
                result.names_[myind] = p;
                auto [sig, card] = bed2sketch(p, opts);
                result.cardinalities_[myind] = card;
                std::copy(sig.begin(), sig.end(), &result.signatures_[myind * opts.sketchsize_]);
            }
        } else {
            // BigWig sketching is parallelized within files
            if(opts.by_chrom_) {
                std::vector<flat_hash_map<std::string, std::vector<RegT>>> bc(npaths);
                std::vector<flat_hash_map<std::string, double>> dbc(npaths);
                result.nperfile_.resize(npaths);
                for(size_t i = 0; i < npaths; ++i) {
                    auto &p = paths[i];
                    auto res = bw2sketch(p, opts, /*parallel_process=*/true);
                    result.names_[i] = p;
                    result.cardinalities_[i] = res.card_;
                    result.nperfile_[i] = res.chrmap_->size();
                    bc[i] = std::move(*res.chrmap_.get());
                    dbc[i] = std::move(res.cardmap_);
                    DBG_ONLY(std::fprintf(stderr, "Cardinality %g found from path %s/%zu\n", res.card_, p.data(), i);)
                    //std::copy(sigs.begin(), sigs.end(), &result.signatures_[opts.sketchsize_* i]);
                }
                const size_t total_n = std::accumulate(result.nperfile_.begin(), result.nperfile_.end(), size_t(0));
                size_t offset = 0;
                result.signatures_.resize(total_n * opts.sketchsize_);
                result.names_.resize(total_n);
                result.cardinalities_.resize(total_n);
                for(size_t i = 0; i < npaths; ++i) {
                    size_t myi = 0;
                    for(const auto &pair: bc[i]) {
                        const auto om = offset + myi++;
                        std::copy(pair.second.begin(), pair.second.end(), &result.signatures_[opts.sketchsize_ * om]);
                        result.names_[om] = paths[i] + ":" + pair.first;
                        result.cardinalities_[om] = dbc[i][pair.first];
                    }
                    offset += bc[i].size();
                }
            } else {
                OMP_PFOR_DYN
                for(size_t i = 0; i < npaths; ++i) {
                    auto myind = filesizes.size() ? filesizes[i].second: uint64_t(i);
                    auto &p(paths[myind]);
                    result.names_[myind] = p;
                    auto res = bw2sketch(p, opts, /*parallel_process=*/false);
                    std::copy(res.global_->begin(), res.global_->end(), &result.signatures_[myind * opts.sketchsize_]);
                    result.cardinalities_[myind] = res.card_;
                }
            }
        }
    }
    std::FILE *ofp;
    if(opts.kmer_result_ == FULL_MMER_SEQUENCE) {
        if((ofp = bfopen(outfile.data(), "r+")) == nullptr) THROW_EXCEPTION(std::runtime_error("Failed to open output file for mmer sequence results."));
        size_t offset = result.names_.size();
        checked_fwrite(&offset, sizeof(offset), 1, ofp);
        {
            const uint32_t k = opts.k_, w = opts.w_;
            checked_fwrite(&k, sizeof(k), 1, ofp);
            checked_fwrite(&w, sizeof(w), 1, ofp);
            uint32_t dtype = (uint32_t)opts.input_mode() | (int(opts.canonicalize()) << 8);
            checked_fwrite(&dtype, sizeof(dtype), 1, ofp);
        }
        checked_fwrite(result.cardinalities_.data(), sizeof(double), result.cardinalities_.size(), ofp);
        offset = 0;
        for(size_t i = 0; i < result.nperfile_.size(); ++i) {
            if(result.nperfile_[i]) {
                checked_fwrite(&result.signatures_.at(offset), sizeof(RegT), result.nperfile_[i], ofp);
                offset += result.nperfile_.at(i);
            }
        }
        std::fclose(ofp);
    } else {
        if(outfile.size() && outfile != "/dev/stdout" && outfile != "-" && (opts.dtype_ == DataType::BED || opts.dtype_ == DataType::BIGWIG)) {
            // Interval sketches are held in memory rather than mapped to outfile, so write the whole stacked file
            const uint64_t t = result.cardinalities_.size(), sketchsize = opts.sketchsize_;
            std::FILE *fp = bfopen(outfile.data(), "wb");
            if(!fp) THROW_EXCEPTION(std::runtime_error("Failed to open file "s + outfile + " for writing"));
            checked_fwrite(fp, &t, sizeof(t));
            checked_fwrite(fp, &sketchsize, sizeof(sketchsize));
            checked_fwrite(fp, result.cardinalities_.data(), result.cardinalities_.size() * sizeof(double));
            checked_fwrite(fp, result.signatures_.data(), result.signatures_.size() * sizeof(RegT));
            std::fclose(fp);
        } else if(outfile.size() && outfile != "/dev/stdout" && outfile != "-") {
            // This should not overlap with the memory mapped for result.signatures_
            const uint64_t t = result.cardinalities_.size();
            std::FILE *fp = bfopen(outfile.data(), "r+");
            if(!fp) THROW_EXCEPTION(std::runtime_error("Failed to open file "s + outfile + " for in-place modification"));
            const uint64_t sketchsize = opts.sketchsize_;
            checked_fwrite(fp, &t, sizeof(t));
            checked_fwrite(fp, &sketchsize, sizeof(sketchsize));
            checked_fwrite(fp, result.cardinalities_.data(), result.cardinalities_.size() * sizeof(double));
            std::fclose(fp);
        } else {
            if(outfile.size() && (outfile == "-" || outfile == "/dev/stdout")) {
                THROW_EXCEPTION(std::runtime_error("Not yet supported: writing stacked sketches to file streams. This may change."));
            }
        }
    }
    if(!outfile.empty() && result.names_.size()) {
        if((ofp = bfopen((outfile + ".names.txt").data(), "wb")) == nullptr)
            THROW_EXCEPTION(std::runtime_error(std::string("Failed to open outfile at ") + outfile + ".names.txt"));
        std::fputs("#Name\tCardinality\n", ofp);
        for(size_t i = 0; i < result.names_.size(); ++i) {
            const auto &n(result.names_[i]);
            checked_fwrite(n.data(), 1, n.size(), ofp);
            if(result.cardinalities_.size()) {
                std::fprintf(ofp, tfmt<double>, result.cardinalities_[i]);
            }
            if(auto &kf = result.kmercountfiles_; !kf.empty())
                std::fputc('\t', ofp), checked_fwrite(kf[i].data(), 1, kf[i].size(), ofp);
            std::fputc('\n', ofp);
        }
        std::fclose(ofp);
    }
    if(!outfile.empty() && result.kmercounts_.size()) {
        const size_t nb = result.kmercounts_.size() * sizeof(decltype(result.kmercounts_)::value_type);
        DBG_ONLY(std::fprintf(stderr, "Writing kmercounts of size %zu\n", result.kmercounts_.size());)
        if((ofp = bfopen((outfile + ".kmercounts.f64").data(), "wb"))) {
            checked_fwrite(result.kmercounts_.data(), nb, 1, ofp);
            std::fclose(ofp);
        } else {
            DBG_ONLY(std::fprintf(stderr, "Failed to open file at %s to write k-mer counts, failing silently.\n", (outfile + ".kmercounts.f64").data());)
        }
    }
    return result;
}

std::vector<std::pair<size_t, uint64_t>> get_filesizes(const std::vector<std::string> &paths) {
    const size_t npaths = paths.size();
    std::vector<std::pair<size_t, uint64_t>> filesizes(npaths);
    OMP_PFOR
    for(size_t i = 0; i < npaths; ++i) {
        filesizes[i] = {nbytes_from_line(paths[i]), uint64_t(i)};
    }
    std::sort(filesizes.begin(), filesizes.end(), std::greater<>());
    return filesizes;
}

}
