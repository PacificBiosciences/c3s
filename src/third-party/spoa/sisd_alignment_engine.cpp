/*!
 * @file sisd_alignment_engine.cpp
 *
 * @brief SisdAlignmentEngine class source file
 */

#include <algorithm>
#include <limits>

#include "graph.hpp"
#include "sisd_alignment_engine.hpp"

namespace spoa {

constexpr int32_t kNegativeInfinity = std::numeric_limits<int32_t>::min() + 1024;

std::unique_ptr<AlignmentEngine> createSisdAlignmentEngine(AlignmentType alignment_type,
                                                           int8_t match, int8_t mismatch,
                                                           int8_t gap)
{

    return std::unique_ptr<AlignmentEngine>(
        new SisdAlignmentEngine(alignment_type, match, mismatch, gap));
}

struct SisdAlignmentEngine::Implementation
{
    std::vector<uint32_t> node_id_to_rank;
    std::vector<int32_t> sequence_profile;
    std::vector<int32_t> M;

    Implementation() : node_id_to_rank(), sequence_profile(), M() {}
};

SisdAlignmentEngine::SisdAlignmentEngine(AlignmentType alignment_type, int8_t match,
                                         int8_t mismatch, int8_t gap)
    : AlignmentEngine(alignment_type, match, mismatch, gap), pimpl_(new Implementation())
{}

SisdAlignmentEngine::~SisdAlignmentEngine() {}

void SisdAlignmentEngine::prealloc(uint32_t max_sequence_size, uint32_t alphabet_size)
{

    this->realloc(max_sequence_size, alphabet_size * max_sequence_size, alphabet_size);
}

void SisdAlignmentEngine::realloc(uint32_t matrix_width, uint32_t matrix_height, uint32_t num_codes)
{

    if (pimpl_->node_id_to_rank.size() < matrix_height - 1) {
        pimpl_->node_id_to_rank.resize(matrix_height - 1, 0);
    }
    if (pimpl_->sequence_profile.size() < num_codes * matrix_width) {
        pimpl_->sequence_profile.resize(num_codes * matrix_width, 0);
    }
    if (pimpl_->M.size() < matrix_height * matrix_width) {
        pimpl_->M.resize(matrix_width * matrix_height, 0);
    }
}

void SisdAlignmentEngine::initialize(const char* sequence, uint32_t sequence_size,
                                     const std::unique_ptr<Graph>& graph) noexcept
{

    uint32_t matrix_width = sequence_size + 1;
    uint32_t matrix_height = graph->nodes().size() + 1;

    for (uint32_t i = 0; i < graph->num_codes(); ++i) {
        char c = graph->decoder(i);
        pimpl_->sequence_profile[i * matrix_width] = 0;
        for (uint32_t j = 0; j < sequence_size; ++j) {
            pimpl_->sequence_profile[i * matrix_width + (j + 1)] =
                (c == sequence[j] ? match_ : mismatch_);
        }
    }

    const auto& rank_to_node_id = graph->rank_to_node_id();

    for (uint32_t i = 0; i < rank_to_node_id.size(); ++i) {
        pimpl_->node_id_to_rank[rank_to_node_id[i]] = i;
    }

    // vertical conditions
    if (alignment_type_ == AlignmentType::kSW || alignment_type_ == AlignmentType::kOV) {
        for (uint32_t i = 0; i < matrix_height; ++i) {
            pimpl_->M[i * matrix_width] = 0;
        }
    } else if (alignment_type_ == AlignmentType::kNW) {
        pimpl_->M[0] = 0;
        for (const auto& node_id : rank_to_node_id) {
            uint32_t i = pimpl_->node_id_to_rank[node_id] + 1;
            const auto& node = graph->nodes()[node_id];
            if (node->in_edges().empty()) {
                pimpl_->M[i * matrix_width] = gap_;
            } else {
                int32_t penalty = kNegativeInfinity;
                for (const auto& edge : node->in_edges()) {
                    uint32_t pred_i = pimpl_->node_id_to_rank[edge->begin_node_id()] + 1;
                    penalty = std::max(penalty, pimpl_->M[pred_i * matrix_width]);
                }
                pimpl_->M[i * matrix_width] = penalty + gap_;
            }
        }
    }

    // horizontal conditions
    if (alignment_type_ == AlignmentType::kSW) {
        for (uint32_t j = 0; j < matrix_width; ++j) {
            pimpl_->M[j] = 0;
        }
    } else if (alignment_type_ == AlignmentType::kOV || alignment_type_ == AlignmentType::kNW) {
        for (uint32_t j = 1; j < matrix_width; ++j) {
            pimpl_->M[j] = j * gap_;
        }
    }
}

Alignment SisdAlignmentEngine::align_sequence_with_graph(const char* sequence,
                                                         uint32_t sequence_size,
                                                         const std::unique_ptr<Graph>& graph)
{

    if (graph->nodes().empty() || sequence_size == 0) {
        return Alignment();
    }

    return align(sequence, sequence_size, graph);
}

Alignment SisdAlignmentEngine::align(const char* sequence, uint32_t sequence_size,
                                     const std::unique_ptr<Graph>& graph) noexcept
{

    uint32_t matrix_width = sequence_size + 1;
    uint32_t matrix_height = graph->nodes().size() + 1;
    const auto& rank_to_node_id = graph->rank_to_node_id();

    // realloc
    this->realloc(matrix_width, matrix_height, graph->num_codes());

    // initialize
    this->initialize(sequence, sequence_size, graph);

    int32_t max_score = alignment_type_ == AlignmentType::kSW ? 0 : kNegativeInfinity;
    int32_t max_i = -1;
    int32_t max_j = -1;
    auto update_max_score = [&max_score, &max_i, &max_j](int32_t* M_row, uint32_t i,
                                                         uint32_t j) -> void {
        if (max_score < M_row[j]) {
            max_score = M_row[j];
            max_i = i;
            max_j = j;
        }
        return;
    };

    // alignment
    for (uint32_t node_id : rank_to_node_id) {
        const auto& node = graph->nodes()[node_id];
        const auto& char_profile = &(pimpl_->sequence_profile[node->code() * matrix_width]);

        uint32_t i = pimpl_->node_id_to_rank[node_id] + 1;

        int32_t* M_row = &(pimpl_->M[i * matrix_width]);

        uint32_t pred_i = node->in_edges().empty()
                              ? 0
                              : pimpl_->node_id_to_rank[node->in_edges()[0]->begin_node_id()] + 1;

        int32_t* M_pred_row = &(pimpl_->M[pred_i * matrix_width]);

        for (uint32_t j = 1; j < matrix_width; ++j) {
            // update M
            M_row[j] = std::max(M_pred_row[j - 1] + char_profile[j], M_pred_row[j] + gap_);
        }

        // check other predeccessors
        for (uint32_t p = 1; p < node->in_edges().size(); ++p) {
            pred_i = pimpl_->node_id_to_rank[node->in_edges()[p]->begin_node_id()] + 1;

            M_pred_row = &(pimpl_->M[pred_i * matrix_width]);

            for (uint32_t j = 1; j < matrix_width; ++j) {
                // update M
                M_row[j] = std::max(M_pred_row[j - 1] + char_profile[j],
                                    std::max(M_row[j], M_pred_row[j] + gap_));
            }
        }

        for (uint32_t j = 1; j < matrix_width; ++j) {
            // update M
            M_row[j] = std::max(M_row[j - 1] + gap_, M_row[j]);

            if (alignment_type_ == AlignmentType::kSW) {
                M_row[j] = std::max(M_row[j], 0);
                update_max_score(M_row, i, j);

            } else if (alignment_type_ == AlignmentType::kNW &&
                       (j == matrix_width - 1 && node->out_edges().empty())) {
                update_max_score(M_row, i, j);

            } else if (alignment_type_ == AlignmentType::kOV && (node->out_edges().empty())) {
                update_max_score(M_row, i, j);
            }
        }
    }

    // backtrack
    Alignment alignment;

    uint32_t i = max_i;
    uint32_t j = max_j;

    auto sw_condition = [this, &i, &j, &matrix_width]() {
        return (pimpl_->M[i * matrix_width + j] == 0) ? false : true;
    };
    auto nw_condition = [&i, &j]() { return (i == 0 && j == 0) ? false : true; };
    auto ov_condition = [&i, &j]() { return (i == 0 || j == 0) ? false : true; };

    uint32_t prev_i = 0;
    uint32_t prev_j = 0;

    while ((alignment_type_ == AlignmentType::kSW && sw_condition()) ||
           (alignment_type_ == AlignmentType::kNW && nw_condition()) ||
           (alignment_type_ == AlignmentType::kOV && ov_condition())) {

        auto M_ij = pimpl_->M[i * matrix_width + j];
        bool predecessor_found = false;

        if (i != 0 && j != 0) {
            const auto& node = graph->nodes()[rank_to_node_id[i - 1]];
            int32_t match_cost = pimpl_->sequence_profile[node->code() * matrix_width + j];

            uint32_t pred_i =
                node->in_edges().empty()
                    ? 0
                    : pimpl_->node_id_to_rank[node->in_edges()[0]->begin_node_id()] + 1;

            if (M_ij == pimpl_->M[pred_i * matrix_width + (j - 1)] + match_cost) {
                prev_i = pred_i;
                prev_j = j - 1;
                predecessor_found = true;
            }

            if (!predecessor_found) {
                const auto& edges = node->in_edges();
                for (uint32_t p = 1; p < edges.size(); ++p) {
                    uint32_t pred_i = pimpl_->node_id_to_rank[edges[p]->begin_node_id()] + 1;

                    if (M_ij == pimpl_->M[pred_i * matrix_width + (j - 1)] + match_cost) {
                        prev_i = pred_i;
                        prev_j = j - 1;
                        predecessor_found = true;
                        break;
                    }
                }
            }
        }

        if (!predecessor_found && i != 0) {
            const auto& node = graph->nodes()[rank_to_node_id[i - 1]];

            uint32_t pred_i =
                node->in_edges().empty()
                    ? 0
                    : pimpl_->node_id_to_rank[node->in_edges()[0]->begin_node_id()] + 1;

            if (M_ij == pimpl_->M[pred_i * matrix_width + j] + gap_) {
                prev_i = pred_i;
                prev_j = j;
                predecessor_found = true;
            }

            if (!predecessor_found) {
                const auto& edges = node->in_edges();
                for (uint32_t p = 1; p < edges.size(); ++p) {
                    uint32_t pred_i = pimpl_->node_id_to_rank[edges[p]->begin_node_id()] + 1;

                    if (M_ij == pimpl_->M[pred_i * matrix_width + j] + gap_) {
                        prev_i = pred_i;
                        prev_j = j;
                        predecessor_found = true;
                        break;
                    }
                }
            }
        }

        if (!predecessor_found && M_ij == pimpl_->M[i * matrix_width + j - 1] + gap_) {
            prev_i = i;
            prev_j = j - 1;
            predecessor_found = true;
        }

        alignment.emplace_back(i == prev_i ? -1 : rank_to_node_id[i - 1], j == prev_j ? -1 : j - 1);

        i = prev_i;
        j = prev_j;
    }

    std::reverse(alignment.begin(), alignment.end());
    return alignment;
}

}  // namespace spoa
