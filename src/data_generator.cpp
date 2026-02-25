/**
 * @file data_generator.cpp
 * @brief Implementation of synthetic-data generator.
 */

#include "motif/data_generator.hpp"
#include <algorithm>

namespace motif {

DataGenerator::DataGenerator(unsigned int seed) : rng_(seed) {}

char DataGenerator::random_base() {
    constexpr char bases[] = {'A', 'C', 'G', 'T'};
    return bases[base_dist_(rng_)];
}

Sequence DataGenerator::generate_random_sequence(size_t length) {
    Sequence seq;
    seq.reserve(length);
    for (size_t i = 0; i < length; ++i) {
        seq.push_back(random_base());
    }
    return seq;
}

SequenceList DataGenerator::generate_random_sequences(size_t count, size_t length) {
    SequenceList sequences;
    sequences.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        sequences.push_back(generate_random_sequence(length));
    }
    return sequences;
}

SequenceList DataGenerator::inject_motif(
    const SequenceList& sequences,
    const std::string& motif,
    double injection_rate
) {
    double rate = injection_rate;
    if (rate < 0.0) {
        rate = 0.0;
    }
    if (rate > 1.0) {
        rate = 1.0;
    }
    if (motif.empty() || rate <= 0.0) {
        return sequences;
    }

    SequenceList result;
    result.reserve(sequences.size());

    for (const auto& seq : sequences) {
        if (seq.size() < motif.size()) {
            result.push_back(seq);
            continue;
        }

        if (prob_dist_(rng_) < rate) {
            std::uniform_int_distribution<size_t> pos_dist(0, seq.size() - motif.size());
            const size_t pos = pos_dist(rng_);
            std::string updated = seq;
            updated.replace(pos, motif.size(), motif);
            result.push_back(std::move(updated));
        } else {
            result.push_back(seq);
        }
    }

    return result;
}

std::pair<SequenceList, SequenceList> DataGenerator::generate_test_data(
    size_t num_sequences,
    size_t seq_length,
    const std::string& motif,
    double injection_rate
) {
    SequenceList primary = generate_random_sequences(num_sequences, seq_length);
    primary = inject_motif(primary, motif, injection_rate);

    SequenceList control = generate_random_sequences(num_sequences, seq_length);
    return {primary, control};
}

} // namespace motif
