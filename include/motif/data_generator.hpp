/**
 * @file data_generator.hpp
 * @brief Synthetic Data Generation for Testing
 * 
 * Utilities for generating test sequences with known embedded motifs.
 * Essential for validating motif discovery algorithms.
 */

#ifndef MOTIF_DATA_GENERATOR_HPP
#define MOTIF_DATA_GENERATOR_HPP

#include "motif/types.hpp"
#include <random>

namespace motif {

/**
 * @brief Random data generator for motif discovery testing
 */
class DataGenerator {
public:
    /**
     * @brief Construct generator with optional seed
     * @param seed Random seed (default uses random device)
     */
    explicit DataGenerator(unsigned int seed = std::random_device{}());
    
    /**
     * @brief Generate a random DNA sequence
     * @param length Sequence length
     * @return Random DNA string
     */
    Sequence generate_random_sequence(size_t length);
    
    /**
     * @brief Generate multiple random sequences
     * @param count Number of sequences
     * @param length Length of each sequence
     * @return List of random DNA sequences
     */
    SequenceList generate_random_sequences(size_t count, size_t length);
    
    /**
     * @brief Inject a motif into sequences
     * 
     * @param sequences Input sequences
     * @param motif Motif to inject
     * @param injection_rate Probability of injection per sequence (0.0-1.0)
     * @return New sequence list with motif injected
     */
    SequenceList inject_motif(
        const SequenceList& sequences,
        const std::string& motif,
        double injection_rate = 0.5
    );
    
    /**
     * @brief Generate primary and control sets for testing
     * 
     * @param num_sequences Number of sequences per set
     * @param seq_length Length of each sequence
     * @param motif Motif to inject in primary set
     * @param injection_rate Injection rate for primary set
     * @return Pair of (primary_sequences, control_sequences)
     */
    std::pair<SequenceList, SequenceList> generate_test_data(
        size_t num_sequences,
        size_t seq_length,
        const std::string& motif,
        double injection_rate = 0.5
    );

private:
    std::mt19937 rng_;
    std::uniform_int_distribution<int> base_dist_{0, 3};
    std::uniform_real_distribution<double> prob_dist_{0.0, 1.0};
    
    char random_base();
};

} // namespace motif

#endif // MOTIF_DATA_GENERATOR_HPP
