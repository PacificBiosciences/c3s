/*!
 * @file sisd_alignment_engine.hpp
 *
 * @brief SisdAlignmentEngine class header file
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

class SisdAlignmentEngine;
std::unique_ptr<AlignmentEngine> createSisdAlignmentEngine(AlignmentType alignment_type,
                                                           int8_t match, int8_t mismatch,
                                                           int8_t gap);

class SisdAlignmentEngine : public AlignmentEngine
{
public:
    ~SisdAlignmentEngine();

    void prealloc(uint32_t max_sequence_size, uint32_t alphabet_size) override;

    Alignment align_sequence_with_graph(const char* sequence, uint32_t sequence_size,
                                        const std::unique_ptr<Graph>& graph) override;

    friend std::unique_ptr<AlignmentEngine> createSisdAlignmentEngine(AlignmentType alignment_type,
                                                                      int8_t match, int8_t mismatch,
                                                                      int8_t gap);

private:
    SisdAlignmentEngine(AlignmentType alignment_type, int8_t match, int8_t mismatch, int8_t gap);
    SisdAlignmentEngine(const SisdAlignmentEngine&) = delete;
    const SisdAlignmentEngine& operator=(const SisdAlignmentEngine&) = delete;

    Alignment align(const char* sequence, uint32_t sequence_size,
                    const std::unique_ptr<Graph>& graph) noexcept;

    void realloc(uint32_t matrix_width, uint32_t matrix_height, uint32_t num_codes);

    void initialize(const char* sequence, uint32_t sequence_size,
                    const std::unique_ptr<Graph>& graph) noexcept;

    struct Implementation;
    std::unique_ptr<Implementation> pimpl_;
};

}  // namespace spoa
