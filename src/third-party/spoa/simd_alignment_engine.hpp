/*!
 * @file simd_alignment_engine.hpp
 *
 * @brief SimdAlignmentEngine class header file
 */

#pragma once

#include <stdint.h>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "alignment_engine.hpp"

namespace spoa {

class Graph;

class SimdAlignmentEngine;
std::unique_ptr<AlignmentEngine> createSimdAlignmentEngine(AlignmentType alignment_type,
                                                           int8_t match, int8_t mismatch,
                                                           int8_t gap);

class SimdAlignmentEngine : public AlignmentEngine
{
public:
    ~SimdAlignmentEngine();

    void prealloc(uint32_t max_sequence_size, uint32_t alphabet_size) override;

    Alignment align_sequence_with_graph(const char* sequence, uint32_t sequence_size,
                                        const std::unique_ptr<Graph>& graph) override;

    friend std::unique_ptr<AlignmentEngine> createSimdAlignmentEngine(AlignmentType alignment_type,
                                                                      int8_t match, int8_t mismatch,
                                                                      int8_t gap);

private:
    SimdAlignmentEngine(AlignmentType alignment_type, int8_t match, int8_t mismatch, int8_t gap);
    SimdAlignmentEngine(const SimdAlignmentEngine&) = delete;
    const SimdAlignmentEngine& operator=(const SimdAlignmentEngine&) = delete;

    template <typename T>
    Alignment align(const char* sequence, uint32_t sequence_size,
                    const std::unique_ptr<Graph>& graph) noexcept;

    void realloc(uint32_t matrix_width, uint32_t matrix_height, uint32_t num_codes);

    template <typename T>
    void initialize(const char* sequence, const std::unique_ptr<Graph>& graph,
                    uint32_t normal_matrix_width, uint32_t matrix_width,
                    uint32_t matrix_height) noexcept;

    struct Implementation;
    std::unique_ptr<Implementation> pimpl_;
};

}  // namespace spoa
