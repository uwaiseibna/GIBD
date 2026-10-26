/*
 * gfa2gbwt.cpp  —  build a bidirectional GBWT from the P- and W-lines of a GFA
 *
 * Reads the path (P) AND walk (W) lines of a GFA, encodes each as an oriented
 * node sequence, and builds a bidirectional DynamicGBWT (interleaved fwd/rev,
 * the form gibd expects). Each path becomes one GBWT sequence with
 * full PanSN metadata so gibd can reconstruct sample#hap#contig names:
 *   - P-line  <name>           -> sample=<name>, contig="", hap=0, offset=0
 *                                 (each haplotype thread keeps a DISTINCT sample
 *                                  name, so merge mode does not collapse them).
 *   - W-line  <sample> <hap> <contig> <start> <end> <walk>
 *                              -> sample,hap,contig as given, offset=<start>
 *                                 (a chunked HPRC reference appears as W-lines;
 *                                  --ref-prefix GRCh38#0#chr22 then matches and
 *                                  the subrange start gives the bp offset).
 *
 * vg's `vg gbwt -E` collapses every haplotype sample name to _gbwt_ref; building
 * the GBWT here instead preserves per-thread identity AND the reference path, so
 * the same chunk can be swept in merge (cM) mode with correct thread-thread
 * labels and a working reference coordinate system.
 *
 * Node ids are taken verbatim from the GFA segment ids, so the resulting GBWT
 * lines up with the same GFA passed to gibd --gfa for node lengths.
 *
 * Usage:  gfa2gbwt  graph.gfa  out_basename
 *   -> writes out_basename.gbwt
 *
 * Build (see Makefile gfa2gbwt rule), or:
 *   g++ -std=c++17 -O2 -I deps/gbwt/include -I $SDSL_PREFIX/include \
 *       -L deps/gbwt/lib -L $SDSL_PREFIX/lib \
 *       gfa2gbwt.cpp -o gfa2gbwt -lgbwt -lsdsl -lzstd
 */

#include <gbwt/dynamic_gbwt.h>
#include <gbwt/gbwt.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " graph.gfa out_basename\n";
        return 1;
    }
    const std::string gfa_path = argv[1];
    const std::string out_base = argv[2];

    std::ifstream gfa(gfa_path);
    if (!gfa) throw std::runtime_error("Cannot open GFA: " + gfa_path);

    using vt = gbwt::vector_type::value_type;

    struct PathRec {
        gbwt::vector_type nodes;
        std::string sample;
        std::string contig;   // "" for P-lines
        size_t      phase = 0; // haplotype index
        size_t      offset = 0; // subrange start bp (W-lines)
    };
    std::vector<PathRec> recs;

    std::string line;
    while (std::getline(gfa, line)) {
        if (line.empty()) continue;

        if (line[0] == 'P') {
            // P <name> <seg1+,seg2-,...> <overlaps>
            std::stringstream ss(line);
            std::string tag, name, seglist;
            std::getline(ss, tag, '\t');
            std::getline(ss, name, '\t');
            std::getline(ss, seglist, '\t');
            if (seglist.empty()) continue;

            PathRec r;
            r.sample = name;
            std::stringstream seg_ss(seglist);
            std::string tok;
            while (std::getline(seg_ss, tok, ',')) {
                if (tok.size() < 2) continue;
                char orient = tok.back();
                int id = std::atoi(tok.substr(0, tok.size() - 1).c_str());
                bool rev = (orient == '-');
                r.nodes.push_back(static_cast<vt>(gbwt::Node::encode(id, rev)));
            }
            if (!r.nodes.empty()) recs.push_back(std::move(r));

        } else if (line[0] == 'W') {
            // W <sample> <hapidx> <contig> <start> <end> <walk(>id<id...)>
            std::stringstream ss(line);
            std::string tag, sample, hapstr, contig, startstr, endstr, walk;
            std::getline(ss, tag, '\t');
            std::getline(ss, sample, '\t');
            std::getline(ss, hapstr, '\t');
            std::getline(ss, contig, '\t');
            std::getline(ss, startstr, '\t');
            std::getline(ss, endstr, '\t');
            std::getline(ss, walk, '\t');
            if (walk.empty()) continue;

            PathRec r;
            r.sample = sample;
            r.contig = contig;
            r.phase  = static_cast<size_t>(std::atol(hapstr.c_str()));
            r.offset = static_cast<size_t>(std::atol(startstr.c_str()));
            // walk = sequence of (orient char)(segment id), orient in {'>','<'}
            size_t i = 0;
            while (i < walk.size()) {
                char o = walk[i];
                if (o != '>' && o != '<') { ++i; continue; }
                ++i;
                size_t j = i;
                while (j < walk.size() && walk[j] != '>' && walk[j] != '<') ++j;
                if (j > i) {
                    int id = std::atoi(walk.substr(i, j - i).c_str());
                    bool rev = (o == '<');
                    r.nodes.push_back(static_cast<vt>(gbwt::Node::encode(id, rev)));
                }
                i = j;
            }
            if (!r.nodes.empty()) recs.push_back(std::move(r));
        }
    }
    if (recs.empty())
        throw std::runtime_error("No P- or W-lines found in " + gfa_path);

    // ---- Build bidirectional DynamicGBWT (interleave fwd $ rev $ per path) ----
    gbwt::vector_type text;
    for (const auto& r : recs) {
        text.insert(text.end(), r.nodes.begin(), r.nodes.end());
        text.push_back(static_cast<vt>(gbwt::ENDMARKER));
        gbwt::vector_type rev;
        gbwt::reversePath(r.nodes, rev);
        text.insert(text.end(), rev.begin(), rev.end());
        text.push_back(static_cast<vt>(gbwt::ENDMARKER));
    }
    gbwt::DynamicGBWT dynamic;
    dynamic.insert(text, /*has_both_orientations=*/true);

    // ---- Metadata: full PanSN (distinct sample names + contig names + offset) ----
    std::vector<std::string> sample_names, contig_names;
    std::map<std::string, gbwt::size_type> sample_idx, contig_idx;
    auto intern = [](const std::string& s, std::vector<std::string>& names,
                     std::map<std::string, gbwt::size_type>& idx) -> gbwt::size_type {
        auto it = idx.find(s);
        if (it != idx.end()) return it->second;
        gbwt::size_type id = static_cast<gbwt::size_type>(names.size());
        names.push_back(s);
        idx[s] = id;
        return id;
    };

    size_t max_hap = 0;
    for (const auto& r : recs) {
        gbwt::size_type s = intern(r.sample, sample_names, sample_idx);
        gbwt::size_type c = intern(r.contig, contig_names, contig_idx);
        dynamic.metadata.addPath(s, c,
                                 static_cast<gbwt::size_type>(r.phase),
                                 static_cast<gbwt::size_type>(r.offset));
        if (r.phase + 1 > max_hap) max_hap = r.phase + 1;
    }
    dynamic.metadata.setSamples(sample_names);
    dynamic.metadata.setContigs(contig_names);
    dynamic.metadata.setHaplotypes(max_hap > 0 ? max_hap : 1);
    dynamic.header.set(gbwt::GBWTHeader::FLAG_METADATA);

    // ---- Compress and save ----
    gbwt::GBWT index(dynamic);
    const std::string out_file = out_base + ".gbwt";
    std::ofstream out(out_file, std::ios::binary);
    if (!out) throw std::runtime_error("Cannot write " + out_file);
    index.serialize(out);

    std::cerr << "Wrote " << out_file << ": " << recs.size() << " paths ("
              << sample_names.size() << " samples, " << contig_names.size()
              << " contigs), " << index.sequences() << " sequences, sigma="
              << index.sigma() << "\n";
    return 0;
}
