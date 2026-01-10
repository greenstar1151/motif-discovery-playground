/**
 * @file data_generator.cpp
 * @brief Implementation of Data Generator
 */

#include "motif/data_generator.hpp"

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
        seq += random_base();
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
    SequenceList result;
    result.reserve(sequences.size());
    
    for (const auto& seq : sequences) {
        if (seq.length() < motif.length()) {
            result.push_back(seq);
            continue;
        }
        
        if (prob_dist_(rng_) < injection_rate) {
            // Inject motif at random position
            std::uniform_int_distribution<size_t> pos_dist(
                0, seq.length() - motif.length()
            );
            size_t pos = pos_dist(rng_);
            
            std::string new_seq = seq;
            new_seq.replace(pos, motif.length(), motif);
            result.push_back(new_seq);
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
    // Generate primary set with injected motif
    SequenceList primary = generate_random_sequences(num_sequences, seq_length);
    primary = inject_motif(primary, motif, injection_rate);
    
    // Generate control set (pure random)
    SequenceList control = generate_random_sequences(num_sequences, seq_length);
    
    return {primary, control};
}

} // namespace motif
