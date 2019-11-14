/*!
 * @file simd_alignment_engine.cpp
 *
 * @brief SimdAlignmentEngine class source file
 */

#include <algorithm>
#include <limits>

extern "C"
{
#include <immintrin.h>  // AVX2 and lower
}

#include "graph.hpp"
#include "simd_alignment_engine.hpp"

namespace spoa {

// Taken from https://gcc.gnu.org/viewcvs/gcc?view=revision&revision=216149
inline void* align(size_t __align, size_t __size, void*& __ptr, size_t& __space) noexcept
{

    const auto __intptr = reinterpret_cast<uintptr_t>(__ptr);
    const auto __aligned = (__intptr - 1u + __align) & -__align;
    const auto __diff = __aligned - __intptr;
    if ((__size + __diff) > __space)
        return nullptr;
    else {
        __space -= __diff;
        return __ptr = reinterpret_cast<void*>(__aligned);
    }
}

template <typename T>
T* allocateAlignedMemory(T** storage, uint32_t size, uint32_t alignment)
{
    *storage = new T[size + alignment - 1];
    void* ptr = static_cast<void*>(*storage);
    size_t storage_size = (size + alignment - 1) * sizeof(T);
    return static_cast<T*>(align(alignment, size * sizeof(T), ptr, storage_size));
}

template <typename T>
struct InstructionSet;

#if defined(__AVX2__)

constexpr uint32_t kRegisterSize = 256;
using __mxxxi = __m256i;

inline __mxxxi _mmxxx_load_si(__mxxxi const* mem_addr) { return _mm256_load_si256(mem_addr); }

inline void _mmxxx_store_si(__mxxxi* mem_addr, const __mxxxi& a)
{
    _mm256_store_si256(mem_addr, a);
}

inline __mxxxi _mmxxx_or_si(const __mxxxi& a, const __mxxxi& b) { return _mm256_or_si256(a, b); }

#define _mmxxx_slli_si(a, n)                                                                      \
    n < 16                                                                                        \
        ? _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0, 0, 2, 0)), 16 - n) \
        : _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0, 0, 2, 0))

#define _mmxxx_srli_si(a, n) \
    _mm256_srli_si256(_mm256_permute2x128_si256(a, a, _MM_SHUFFLE(2, 0, 0, 1)), n - 16)

template <>
struct InstructionSet<int16_t>
{
    using type = int16_t;
    static constexpr uint32_t kNumVar = kRegisterSize / (8 * sizeof(type));
    static constexpr uint32_t kLogNumVar = 4;
    static constexpr uint32_t kLSS = 2;   // Left Shift Size
    static constexpr uint32_t kRSS = 30;  // Right Shift Size
    static inline __mxxxi _mmxxx_add_epi(const __mxxxi& a, const __mxxxi& b)
    {
        return _mm256_add_epi16(a, b);
    }
    static inline __mxxxi _mmxxx_sub_epi(const __mxxxi& a, const __mxxxi& b)
    {
        return _mm256_sub_epi16(a, b);
    }
    static inline __mxxxi _mmxxx_min_epi(const __mxxxi& a, const __mxxxi& b)
    {
        return _mm256_min_epi16(a, b);
    }
    static inline __mxxxi _mmxxx_max_epi(const __mxxxi& a, const __mxxxi& b)
    {
        return _mm256_max_epi16(a, b);
    }
    static inline __mxxxi _mmxxx_set1_epi(type a) { return _mm256_set1_epi16(a); }
    static inline void _mmxxx_prefix_max(__mxxxi& a, const __mxxxi* masks, const __mxxxi* penalties)
    {

        a = _mmxxx_max_epi(
            a, _mmxxx_or_si(masks[0], _mmxxx_slli_si(_mmxxx_add_epi(a, penalties[0]), 2)));
        a = _mmxxx_max_epi(
            a, _mmxxx_or_si(masks[1], _mmxxx_slli_si(_mmxxx_add_epi(a, penalties[1]), 4)));
        a = _mmxxx_max_epi(
            a, _mmxxx_or_si(masks[2], _mmxxx_slli_si(_mmxxx_add_epi(a, penalties[2]), 8)));
        a = _mmxxx_max_epi(
            a, _mmxxx_or_si(masks[3], _mmxxx_slli_si(_mmxxx_add_epi(a, penalties[3]), 16)));
    }
};

template <>
struct InstructionSet<int32_t>
{
    using type = int32_t;
    static constexpr uint32_t kNumVar = kRegisterSize / (8 * sizeof(type));
    static constexpr uint32_t kLogNumVar = 3;
    static constexpr uint32_t kLSS = 4;
    static constexpr uint32_t kRSS = 28;
    static inline __mxxxi _mmxxx_add_epi(const __mxxxi& a, const __mxxxi& b)
    {
        return _mm256_add_epi32(a, b);
    }
    static inline __mxxxi _mmxxx_sub_epi(const __mxxxi& a, const __mxxxi& b)
    {
        return _mm256_sub_epi32(a, b);
    }
    static inline __mxxxi _mmxxx_min_epi(const __mxxxi& a, const __mxxxi& b)
    {
        return _mm256_min_epi32(a, b);
    }
    static inline __mxxxi _mmxxx_max_epi(const __mxxxi& a, const __mxxxi& b)
    {
        return _mm256_max_epi32(a, b);
    }
    static inline __mxxxi _mmxxx_set1_epi(type a) { return _mm256_set1_epi32(a); }
    static inline void _mmxxx_prefix_max(__mxxxi& a, const __mxxxi* masks, const __mxxxi* penalties)
    {

        a = _mmxxx_max_epi(
            a, _mmxxx_or_si(masks[0], _mmxxx_slli_si(_mmxxx_add_epi(a, penalties[0]), 4)));
        a = _mmxxx_max_epi(
            a, _mmxxx_or_si(masks[1], _mmxxx_slli_si(_mmxxx_add_epi(a, penalties[1]), 8)));
        a = _mmxxx_max_epi(
            a, _mmxxx_or_si(masks[2], _mmxxx_slli_si(_mmxxx_add_epi(a, penalties[2]), 16)));
    }
};

#elif defined(__SSE4_1__)

constexpr uint32_t kRegisterSize = 128;
using __mxxxi = __m128i;

inline __mxxxi _mmxxx_load_si(__mxxxi const* mem_addr) { return _mm_load_si128(mem_addr); }

inline void _mmxxx_store_si(__mxxxi* mem_addr, const __mxxxi& a) { _mm_store_si128(mem_addr, a); }

inline __mxxxi _mmxxx_or_si(const __mxxxi& a, const __mxxxi& b) { return _mm_or_si128(a, b); }

#define _mmxxx_slli_si(a, n) _mm_slli_si128(a, n)

#define _mmxxx_srli_si(a, n) _mm_srli_si128(a, n)

template <>
struct InstructionSet<int16_t>
{
    using type = int16_t;
    static constexpr uint32_t kNumVar = kRegisterSize / (8 * sizeof(type));
    static constexpr uint32_t kLogNumVar = 3;
    static constexpr uint32_t kLSS = 2;
    static constexpr uint32_t kRSS = 14;
    static inline __mxxxi _mmxxx_add_epi(const __mxxxi& a, const __mxxxi& b)
    {
        return _mm_add_epi16(a, b);
    }
    static inline __mxxxi _mmxxx_sub_epi(const __mxxxi& a, const __mxxxi& b)
    {
        return _mm_sub_epi16(a, b);
    }
    static inline __mxxxi _mmxxx_min_epi(const __mxxxi& a, const __mxxxi& b)
    {
        return _mm_min_epi16(a, b);
    }
    static inline __mxxxi _mmxxx_max_epi(const __mxxxi& a, const __mxxxi& b)
    {
        return _mm_max_epi16(a, b);
    }
    static inline __mxxxi _mmxxx_set1_epi(type a) { return _mm_set1_epi16(a); }
    static inline void _mmxxx_prefix_max(__mxxxi& a, const __mxxxi* masks, const __mxxxi* penalties)
    {

        a = _mmxxx_max_epi(
            a, _mmxxx_or_si(masks[0], _mmxxx_slli_si(_mmxxx_add_epi(a, penalties[0]), 2)));
        a = _mmxxx_max_epi(
            a, _mmxxx_or_si(masks[1], _mmxxx_slli_si(_mmxxx_add_epi(a, penalties[1]), 4)));
        a = _mmxxx_max_epi(
            a, _mmxxx_or_si(masks[2], _mmxxx_slli_si(_mmxxx_add_epi(a, penalties[2]), 8)));
    }
};

template <>
struct InstructionSet<int32_t>
{
    using type = int32_t;
    static constexpr uint32_t kNumVar = kRegisterSize / (8 * sizeof(type));
    static constexpr uint32_t kLogNumVar = 2;
    static constexpr uint32_t kLSS = 4;
    static constexpr uint32_t kRSS = 12;
    static inline __mxxxi _mmxxx_add_epi(const __mxxxi& a, const __mxxxi& b)
    {
        return _mm_add_epi32(a, b);
    }
    static inline __mxxxi _mmxxx_sub_epi(const __mxxxi& a, const __mxxxi& b)
    {
        return _mm_sub_epi32(a, b);
    }
    static inline __mxxxi _mmxxx_min_epi(const __mxxxi& a, const __mxxxi& b)
    {
        return _mm_min_epi32(a, b);
    }
    static inline __mxxxi _mmxxx_max_epi(const __mxxxi& a, const __mxxxi& b)
    {
        return _mm_max_epi32(a, b);
    }
    static inline __mxxxi _mmxxx_set1_epi(type a) { return _mm_set1_epi32(a); }
    static inline void _mmxxx_prefix_max(__mxxxi& a, const __mxxxi* masks, const __mxxxi* penalties)
    {

        a = _mmxxx_max_epi(
            a, _mmxxx_or_si(masks[0], _mmxxx_slli_si(_mmxxx_add_epi(a, penalties[0]), 4)));
        a = _mmxxx_max_epi(
            a, _mmxxx_or_si(masks[1], _mmxxx_slli_si(_mmxxx_add_epi(a, penalties[1]), 8)));
    }
};

#endif

#if defined(__AVX2__) || defined(__SSE4_1__)

template <typename T>
void _mmxxx_print(const __mxxxi& a)
{

    __attribute__((aligned(kRegisterSize / 8))) typename T::type unpacked[T::kNumVar];
    _mmxxx_store_si(reinterpret_cast<__mxxxi*>(unpacked), a);

    for (uint32_t i = 0; i < T::kNumVar; i++) {
        printf("%d ", unpacked[i]);
    }
}

template <typename T>
typename T::type _mmxxx_max_value(const __mxxxi& a)
{

    typename T::type max_score = 0;
    __attribute__((aligned(kRegisterSize / 8))) typename T::type unpacked[T::kNumVar];
    _mmxxx_store_si(reinterpret_cast<__mxxxi*>(unpacked), a);

    for (uint32_t i = 0; i < T::kNumVar; i++) {
        max_score = std::max(max_score, unpacked[i]);
    }

    return max_score;
}

template <typename T>
typename T::type _mmxxx_value_at(const __mxxxi& a, uint32_t i)
{

    __attribute__((aligned(kRegisterSize / 8))) typename T::type unpacked[T::kNumVar];
    _mmxxx_store_si(reinterpret_cast<__mxxxi*>(unpacked), a);

    return unpacked[i];
}

template <typename T>
int32_t _mmxxx_index_of(const __mxxxi* row, uint32_t row_width, typename T::type value)
{

    for (uint32_t i = 0; i < row_width; ++i) {
        __attribute__((aligned(kRegisterSize / 8))) typename T::type unpacked[T::kNumVar];
        _mmxxx_store_si(reinterpret_cast<__mxxxi*>(unpacked), row[i]);

        for (uint32_t j = 0; j < T::kNumVar; j++) {
            if (unpacked[j] == value) {
                return i * T::kNumVar + j;
            }
        }
    }

    return -1;
}

#endif

std::unique_ptr<AlignmentEngine> createSimdAlignmentEngine(AlignmentType alignment_type,
                                                           int8_t match, int8_t mismatch,
                                                           int8_t gap)
{

#if defined(__AVX2__) || defined(__SSE4_1__)

    return std::unique_ptr<AlignmentEngine>(
        new SimdAlignmentEngine(alignment_type, match, mismatch, gap));

#else

    return nullptr;

#endif
}

struct SimdAlignmentEngine::Implementation
{

#if defined(__AVX2__) || defined(__SSE4_1__)

    std::vector<uint32_t> node_id_to_rank;

    std::unique_ptr<__mxxxi[]> sequence_profile_storage;
    uint32_t sequence_profile_size;
    __mxxxi* sequence_profile;

    std::vector<int32_t> first_column;
    std::unique_ptr<__mxxxi[]> M_storage;
    uint32_t M_size;
    __mxxxi* M;

    std::unique_ptr<__mxxxi[]> masks_storage;
    uint32_t masks_size;
    __mxxxi* masks;

    std::unique_ptr<__mxxxi[]> penalties_storage;
    uint32_t penalties_size;
    __mxxxi* penalties;

    Implementation()
        : node_id_to_rank()
        , sequence_profile_storage(nullptr)
        , sequence_profile_size(0)
        , sequence_profile(nullptr)
        , first_column()
        , M_storage(nullptr)
        , M_size(0)
        , M(nullptr)
        , masks_storage(nullptr)
        , masks_size(0)
        , masks(nullptr)
        , penalties_storage(nullptr)
        , penalties_size(0)
        , penalties(nullptr)
    {}

#endif
};

SimdAlignmentEngine::SimdAlignmentEngine(AlignmentType alignment_type, int8_t match,
                                         int8_t mismatch, int8_t gap)
    : AlignmentEngine(alignment_type, match, mismatch, gap), pimpl_(new Implementation())
{}

SimdAlignmentEngine::~SimdAlignmentEngine() {}

void SimdAlignmentEngine::prealloc(uint32_t max_sequence_size, uint32_t alphabet_size)
{

#if defined(__AVX2__) || defined(__SSE4_1__)

    uint32_t longest_path =
        max_sequence_size * (alphabet_size + 1) + 1 + InstructionSet<int16_t>::kNumVar;

    uint32_t max_penalty = std::max(std::max(abs(match_), abs(mismatch_)), abs(gap_));

    if (max_penalty * longest_path < std::numeric_limits<int16_t>::max()) {
        this->realloc((max_sequence_size / InstructionSet<int16_t>::kNumVar) + 1,
                      alphabet_size * max_sequence_size, alphabet_size);
    } else {
        this->realloc((max_sequence_size / InstructionSet<int32_t>::kNumVar) + 1,
                      alphabet_size * max_sequence_size, alphabet_size);
    }

#endif
}

void SimdAlignmentEngine::realloc(uint32_t matrix_width, uint32_t matrix_height, uint32_t num_codes)
{

#if defined(__AVX2__) || defined(__SSE4_1__)

    if (pimpl_->node_id_to_rank.size() < matrix_height - 1) {
        pimpl_->node_id_to_rank.resize(matrix_height - 1, 0);
    }
    if (pimpl_->sequence_profile_size < num_codes * matrix_width) {
        __mxxxi* storage = nullptr;
        pimpl_->sequence_profile_size = num_codes * matrix_width;
        pimpl_->sequence_profile =
            allocateAlignedMemory(&storage, pimpl_->sequence_profile_size, kRegisterSize / 8);
        pimpl_->sequence_profile_storage.reset();
        pimpl_->sequence_profile_storage = std::unique_ptr<__mxxxi[]>(storage);
    }
    if (pimpl_->first_column.size() < matrix_height) {
        pimpl_->first_column.resize(matrix_height, 0);
    }
    if (pimpl_->M_size < matrix_height * matrix_width) {
        __mxxxi* storage = nullptr;
        pimpl_->M_size = matrix_height * matrix_width;
        pimpl_->M = allocateAlignedMemory(&storage, pimpl_->M_size, kRegisterSize / 8);
        pimpl_->M_storage.reset();
        pimpl_->M_storage = std::unique_ptr<__mxxxi[]>(storage);
    }
    if (pimpl_->masks_size < InstructionSet<int16_t>::kLogNumVar + 1) {
        __mxxxi* storage = nullptr;
        pimpl_->masks_size = InstructionSet<int16_t>::kLogNumVar + 1;
        pimpl_->masks = allocateAlignedMemory(&storage, pimpl_->masks_size, kRegisterSize / 8);
        pimpl_->masks_storage.reset();
        pimpl_->masks_storage = std::unique_ptr<__mxxxi[]>(storage);
    }
    if (pimpl_->penalties_size < InstructionSet<int16_t>::kLogNumVar) {
        __mxxxi* storage = nullptr;
        pimpl_->penalties_size = InstructionSet<int16_t>::kLogNumVar;
        pimpl_->penalties =
            allocateAlignedMemory(&storage, pimpl_->penalties_size, kRegisterSize / 8);
        pimpl_->penalties_storage.reset();
        pimpl_->penalties_storage = std::unique_ptr<__mxxxi[]>(storage);
    }

#endif
}

template <typename T>
void SimdAlignmentEngine::initialize(const char* sequence, const std::unique_ptr<Graph>& graph,
                                     uint32_t normal_matrix_width, uint32_t matrix_width,
                                     uint32_t matrix_height) noexcept
{

#if defined(__AVX2__) || defined(__SSE4_1__)

    int32_t padding_penatly = -1 * std::max(std::max(abs(match_), abs(mismatch_)), abs(gap_));

    __attribute__((aligned(kRegisterSize / 8))) typename T::type unpacked[T::kNumVar] = {};

    for (uint32_t i = 0; i < graph->num_codes(); ++i) {
        char c = graph->decoder(i);
        for (uint32_t j = 0; j < matrix_width; ++j) {
            for (uint32_t k = 0; k < T::kNumVar; ++k) {
                unpacked[k] = (j * T::kNumVar + k) < normal_matrix_width
                                  ? (c == sequence[j * T::kNumVar + k] ? match_ : mismatch_)
                                  : padding_penatly;
            }
            pimpl_->sequence_profile[i * matrix_width + j] =
                _mmxxx_load_si(reinterpret_cast<const __mxxxi*>(unpacked));
        }
    }

    const auto& rank_to_node_id = graph->rank_to_node_id();

    for (uint32_t i = 0; i < rank_to_node_id.size(); ++i) {
        pimpl_->node_id_to_rank[rank_to_node_id[i]] = i;
    }

    typename T::type negative_infinity = std::numeric_limits<typename T::type>::min() + 1024;

    __mxxxi zeroes = T::_mmxxx_set1_epi(0);

    // vertical conditions
    if (alignment_type_ == AlignmentType::kSW || alignment_type_ == AlignmentType::kOV) {
        for (uint32_t i = 0; i < matrix_height; ++i) {
            pimpl_->first_column[i] = 0;
        }
    } else if (alignment_type_ == AlignmentType::kNW) {
        pimpl_->first_column[0] = 0;
        for (const auto& node_id : rank_to_node_id) {
            uint32_t i = pimpl_->node_id_to_rank[node_id] + 1;
            const auto& node = graph->nodes()[node_id];
            if (node->in_edges().empty()) {
                pimpl_->first_column[i] = gap_;
            } else {
                int32_t penalty = negative_infinity;
                for (const auto& edge : node->in_edges()) {
                    uint32_t pred_i = pimpl_->node_id_to_rank[edge->begin_node_id()] + 1;
                    penalty = std::max(penalty, pimpl_->first_column[pred_i]);
                }
                pimpl_->first_column[i] = penalty + gap_;
            }
        }
    }

    // horizontal conditions
    if (alignment_type_ == AlignmentType::kSW) {
        for (uint32_t j = 0; j < matrix_width; ++j) {
            pimpl_->M[j] = zeroes;
        }
    } else if (alignment_type_ == AlignmentType::kOV || alignment_type_ == AlignmentType::kNW) {
        for (uint32_t j = 0; j < matrix_width; ++j) {
            pimpl_->M[j] = T::_mmxxx_set1_epi(gap_ + j * T::kNumVar * gap_);
            __mxxxi penalty = T::_mmxxx_set1_epi(gap_);

            for (uint32_t k = 1; k < T::kNumVar; ++k) {
                penalty = _mmxxx_slli_si(penalty, T::kLSS);
                pimpl_->M[j] = T::_mmxxx_add_epi(pimpl_->M[j], penalty);
            }
        }
    }

#endif
}

Alignment SimdAlignmentEngine::align_sequence_with_graph(const char* sequence,
                                                         uint32_t sequence_size,
                                                         const std::unique_ptr<Graph>& graph)
{

    if (graph->nodes().empty() || sequence_size == 0) {
        return Alignment();
    }

#if defined(__AVX2__) || defined(__SSE4_1__)

    uint32_t longest_path =
        graph->nodes().size() + 1 + sequence_size + InstructionSet<int16_t>::kNumVar;

    uint32_t max_penalty = std::max(std::max(abs(match_), abs(mismatch_)), abs(gap_));

    if (max_penalty * longest_path < std::numeric_limits<int16_t>::max()) {
        return align<InstructionSet<int16_t>>(sequence, sequence_size, graph);
    } else {
        return align<InstructionSet<int32_t>>(sequence, sequence_size, graph);
    }

#else

    return Alignment();

#endif
}

template <typename T>
Alignment SimdAlignmentEngine::align(const char* sequence, uint32_t sequence_size,
                                     const std::unique_ptr<Graph>& graph) noexcept
{

#if defined(__AVX2__) || defined(__SSE4_1__)

    uint32_t normal_matrix_width = sequence_size;
    uint32_t matrix_width =
        (sequence_size +
         (sequence_size % T::kNumVar == 0 ? 0 : T::kNumVar - sequence_size % T::kNumVar)) /
        T::kNumVar;
    uint32_t matrix_height = graph->nodes().size() + 1;
    const auto& rank_to_node_id = graph->rank_to_node_id();

    // realloc
    this->realloc(matrix_width, matrix_height, graph->num_codes());

    // initialize
    this->initialize<T>(sequence, graph, normal_matrix_width, matrix_width, matrix_height);

    typename T::type negative_infinity = std::numeric_limits<typename T::type>::min() + 1024;

    __attribute__((aligned(kRegisterSize / 8))) typename T::type unpacked[T::kNumVar] = {0};

    for (uint32_t i = 0, j = 0; i < T::kNumVar && j < T::kLogNumVar; ++i) {
        unpacked[i] = negative_infinity;
        if ((i & (i + 1)) == 0) {
            pimpl_->masks[j++] = _mmxxx_load_si(reinterpret_cast<const __mxxxi*>(unpacked));
        }
    }
    pimpl_->masks[T::kLogNumVar] = _mmxxx_slli_si(T::_mmxxx_set1_epi(negative_infinity), T::kLSS);

    pimpl_->penalties[0] = T::_mmxxx_set1_epi(gap_);
    for (uint32_t i = 1; i < T::kLogNumVar; ++i) {
        pimpl_->penalties[i] =
            T::_mmxxx_add_epi(pimpl_->penalties[i - 1], pimpl_->penalties[i - 1]);
    }

    typename T::type max_score = alignment_type_ == AlignmentType::kSW ? 0 : negative_infinity;
    int32_t max_i = -1;
    int32_t max_j = -1;
    uint32_t last_column_id = (normal_matrix_width - 1) % T::kNumVar;
    __mxxxi zeroes = T::_mmxxx_set1_epi(0);
    __mxxxi penalty = T::_mmxxx_set1_epi(gap_);

    // alignment
    for (uint32_t node_id : rank_to_node_id) {
        const auto& node = graph->nodes()[node_id];
        __mxxxi* char_profile = &(pimpl_->sequence_profile[node->code() * matrix_width]);

        uint32_t i = pimpl_->node_id_to_rank[node_id] + 1;
        uint32_t pred_i = node->in_edges().empty()
                              ? 0
                              : pimpl_->node_id_to_rank[node->in_edges()[0]->begin_node_id()] + 1;

        __mxxxi* M_row = &(pimpl_->M[i * matrix_width]);
        __mxxxi* M_pred_row = &(pimpl_->M[pred_i * matrix_width]);

        __mxxxi x = _mmxxx_srli_si(T::_mmxxx_set1_epi(pimpl_->first_column[pred_i]), T::kRSS);

        for (uint32_t j = 0; j < matrix_width; ++j) {
            // get diagonal
            __mxxxi t1 = _mmxxx_srli_si(M_pred_row[j], T::kRSS);
            M_row[j] = _mmxxx_or_si(_mmxxx_slli_si(M_pred_row[j], T::kLSS), x);
            x = t1;

            // update M
            M_row[j] = T::_mmxxx_max_epi(T::_mmxxx_add_epi(M_row[j], char_profile[j]),
                                         T::_mmxxx_add_epi(M_pred_row[j], penalty));
        }

        // check other predecessors
        for (uint32_t p = 1; p < node->in_edges().size(); ++p) {
            pred_i = pimpl_->node_id_to_rank[node->in_edges()[p]->begin_node_id()] + 1;

            M_pred_row = &(pimpl_->M[pred_i * matrix_width]);

            x = _mmxxx_srli_si(T::_mmxxx_set1_epi(pimpl_->first_column[pred_i]), T::kRSS);

            for (uint32_t j = 0; j < matrix_width; ++j) {
                // get diagonal
                __mxxxi t1 = _mmxxx_srli_si(M_pred_row[j], T::kRSS);
                __mxxxi m = _mmxxx_or_si(_mmxxx_slli_si(M_pred_row[j], T::kLSS), x);
                x = t1;

                // updage M
                M_row[j] = T::_mmxxx_max_epi(
                    M_row[j], T::_mmxxx_max_epi(T::_mmxxx_add_epi(m, char_profile[j]),
                                                T::_mmxxx_add_epi(M_pred_row[j], penalty)));
            }
        }

        __mxxxi score = T::_mmxxx_set1_epi(negative_infinity);
        x = _mmxxx_srli_si(T::_mmxxx_add_epi(T::_mmxxx_set1_epi(pimpl_->first_column[i]), penalty),
                           T::kRSS);

        for (uint32_t j = 0; j < matrix_width; ++j) {

            // add last element of previous vector into this one
            M_row[j] = T::_mmxxx_max_epi(M_row[j], _mmxxx_or_si(x, pimpl_->masks[T::kLogNumVar]));

            T::_mmxxx_prefix_max(M_row[j], pimpl_->masks, pimpl_->penalties);

            x = _mmxxx_srli_si(T::_mmxxx_add_epi(M_row[j], penalty), T::kRSS);

            if (alignment_type_ == AlignmentType::kSW) {
                M_row[j] = T::_mmxxx_max_epi(M_row[j], zeroes);
            }
            score = T::_mmxxx_max_epi(score, M_row[j]);
        }

        if (alignment_type_ == AlignmentType::kSW) {
            int32_t max_row_score = _mmxxx_max_value<T>(score);
            if (max_score < max_row_score) {
                max_score = max_row_score;
                max_i = i;
            }

        } else if (alignment_type_ == AlignmentType::kOV) {
            if (node->out_edges().empty()) {
                int32_t max_row_score = _mmxxx_max_value<T>(score);
                if (max_score < max_row_score) {
                    max_score = max_row_score;
                    max_i = i;
                }
            }

        } else if (alignment_type_ == AlignmentType::kNW) {
            if (node->out_edges().empty()) {
                int32_t max_row_score = _mmxxx_value_at<T>(M_row[matrix_width - 1], last_column_id);
                if (max_score < max_row_score) {
                    max_score = max_row_score;
                    max_i = i;
                }
            }
        }
    }

    if (max_i == -1 && max_j == -1) {  // no alignment found
        return Alignment();
    }

    if (alignment_type_ == AlignmentType::kSW) {
        max_j = _mmxxx_index_of<T>(&(pimpl_->M[max_i * matrix_width]), matrix_width, max_score);

    } else if (alignment_type_ == AlignmentType::kOV) {
        if (graph->nodes()[rank_to_node_id[max_i - 1]]->out_edges().empty()) {
            max_j = _mmxxx_index_of<T>(&(pimpl_->M[max_i * matrix_width]), matrix_width, max_score);
        } else {
            max_j = normal_matrix_width - 1;
        }

    } else if (alignment_type_ == AlignmentType::kNW) {
        max_j = normal_matrix_width - 1;
    }

    // backtrack
    uint32_t max_num_predecessors = 1;
    for (uint32_t i = 0; i < (uint32_t)max_i; ++i) {
        max_num_predecessors = std::max(
            max_num_predecessors, (uint32_t)graph->nodes()[rank_to_node_id[i]]->in_edges().size());
    }

    typename T::type* backtrack_storage = nullptr;
    typename T::type* M = allocateAlignedMemory(
        &backtrack_storage, 3 * T::kNumVar + 2 * T::kNumVar * max_num_predecessors,
        kRegisterSize / 8);
    typename T::type* M_pred = &(M[T::kNumVar]);
    typename T::type* M_diag_pred = &(M_pred[T::kNumVar * max_num_predecessors]);
    typename T::type* M_left_pred = &(M_diag_pred[T::kNumVar * max_num_predecessors]);
    typename T::type* profile = &(M_left_pred[T::kNumVar]);

    std::vector<uint32_t> predecessors;

    int32_t i = max_i;
    int32_t j = max_j;
    int32_t prev_i = 0, prev_j = 0;

    uint32_t j_div = j / T::kNumVar;
    uint32_t j_mod = j % T::kNumVar;

    bool load_next_segment = true;

    Alignment alignment;

    do {
        // check stop condition
        if (j == -1 || i == 0) {
            break;
        }

        const auto& node = graph->nodes()[rank_to_node_id[i - 1]];
        // load everything
        if (load_next_segment) {
            predecessors.clear();

            // load current cells
            _mmxxx_store_si(reinterpret_cast<__mxxxi*>(M), pimpl_->M[i * matrix_width + j_div]);

            // load predecessors cells
            if (node->in_edges().empty()) {
                predecessors.emplace_back(0);
                _mmxxx_store_si(reinterpret_cast<__mxxxi*>(M_pred), pimpl_->M[j_div]);

            } else {
                uint32_t store_pos = 0;
                for (const auto& edge : node->in_edges()) {
                    predecessors.emplace_back(pimpl_->node_id_to_rank[edge->begin_node_id()] + 1);
                    _mmxxx_store_si(reinterpret_cast<__mxxxi*>(&M_pred[store_pos * T::kNumVar]),
                                    pimpl_->M[predecessors.back() * matrix_width + j_div]);
                    ++store_pos;
                }
            }

            // load query profile cells
            _mmxxx_store_si(reinterpret_cast<__mxxxi*>(profile),
                            pimpl_->sequence_profile[node->code() * matrix_width + j_div]);
        }

        // check stop condition
        if (alignment_type_ == AlignmentType::kSW && M[j_mod] == 0) {
            break;
        }

        if (j_mod == 0) {
            // border case
            if (j_div > 0) {
                _mmxxx_store_si(reinterpret_cast<__mxxxi*>(M_left_pred),
                                pimpl_->M[i * matrix_width + j_div - 1]);

                for (uint32_t p = 0; p < predecessors.size(); ++p) {
                    _mmxxx_store_si(reinterpret_cast<__mxxxi*>(&M_diag_pred[p * T::kNumVar]),
                                    pimpl_->M[predecessors[p] * matrix_width + (j_div - 1)]);
                }
            } else {
                M_left_pred[T::kNumVar - 1] = pimpl_->first_column[i];

                for (uint32_t p = 0; p < predecessors.size(); ++p) {
                    M_diag_pred[(p + 1) * T::kNumVar - 1] = pimpl_->first_column[predecessors[p]];
                }
            }
        }

        // find best predecessor cell
        bool predecessor_found = false;

        if (i != 0) {
            for (uint32_t p = 0; p < predecessors.size(); ++p) {
                if ((j_mod == 0 &&
                     M[j_mod] == M_diag_pred[(p + 1) * T::kNumVar - 1] + profile[j_mod]) ||
                    (j_mod != 0 &&
                     M[j_mod] == M_pred[p * T::kNumVar + j_mod - 1] + profile[j_mod])) {

                    prev_i = predecessors[p];
                    prev_j = j - 1;
                    predecessor_found = true;
                    break;
                }
            }
        }

        if (!predecessor_found && i != 0) {
            for (uint32_t p = 0; p < predecessors.size(); ++p) {
                if (M[j_mod] == M_pred[p * T::kNumVar + j_mod] + gap_) {
                    prev_i = predecessors[p];
                    prev_j = j;
                    predecessor_found = true;
                    break;
                }
            }
        }

        if (!predecessor_found) {
            if ((j_mod == 0 && M[j_mod] == M_left_pred[T::kNumVar - 1] + gap_) ||
                (j_mod != 0 && M[j_mod] == M[j_mod - 1] + gap_)) {
                prev_i = i;
                prev_j = j - 1;
                predecessor_found = true;
            }
        }

        alignment.emplace_back(i == prev_i ? -1 : rank_to_node_id[i - 1], j == prev_j ? -1 : j);

        // update for next round
        load_next_segment = (i == prev_i ? false : true) ||
                            (j != prev_j && prev_j % T::kNumVar == T::kNumVar - 1 ? true : false);

        i = prev_i;
        j = prev_j;
        j_div = j / T::kNumVar;
        j_mod = j % T::kNumVar;

    } while (true);

    delete[] backtrack_storage;

    // update alignment for NW (backtrack stops on first row or column)
    if (alignment_type_ == AlignmentType::kNW) {
        while (i == 0 && j != -1) {
            alignment.emplace_back(-1, j);
            --j;
        }
        while (i != 0 && j == -1) {
            alignment.emplace_back(rank_to_node_id[i - 1], -1);

            const auto& node = graph->nodes()[rank_to_node_id[i - 1]];
            if (node->in_edges().empty()) {
                i = 0;
            } else {
                for (const auto& edge : node->in_edges()) {
                    uint32_t pred_i = pimpl_->node_id_to_rank[edge->begin_node_id()] + 1;
                    if (pimpl_->first_column[i] == pimpl_->first_column[pred_i] + gap_) {
                        i = pred_i;
                        break;
                    }
                }
            }
        }
    }

    std::reverse(alignment.begin(), alignment.end());
    return alignment;

#else

    return Alignment();

#endif
}
}  // namespace spoa
