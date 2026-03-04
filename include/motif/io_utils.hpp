/**
 * @file io_utils.hpp
 * @brief I/O and data-preparation utilities.
 *
 * Provides:
 *   - k-mer-preserving Euler-tour sequence shuffling
 *   - Control set generation via shuffling
 *   - JASPAR PWM file loading
 */

#ifndef MOTIF_IO_UTILS_HPP
#define MOTIF_IO_UTILS_HPP

#include "motif/types.hpp"

#include <random>
#include <string>
#include <vector>

namespace motif {

/**
 * @brief Shuffle a sequence preserving k-mer composition (Euler-tour).
 *
 * Uses an Euler path on the (k-1)-mer de Bruijn graph to produce a
 * random permutation that preserves the exact dinucleotide (or higher)
 * frequencies.
 *
 * @param seq   Input DNA sequence.
 * @param rng   Mersenne-Twister RNG instance.
 * @param k     Preserved k-mer order (default 2 = dinucleotide preserving).
 * @return Shuffled sequence.
 */
std::string shuffle_sequence_kmer_preserving(
    const std::string& seq,
    std::mt19937& rng,
    size_t k = 2
);

/**
 * @brief Generate control sequences by k-mer-preserving shuffling.
 *
 * @param input  Sequences to shuffle.
 * @param seed   RNG seed.
 * @param k      Preserved k-mer order (default 2).
 * @return Copy of @p input with each sequence shuffled.
 */
SequenceList generate_control_sequences(
    const SequenceList& input,
    unsigned int seed,
    size_t k = 2
);

/**
 * @brief Load a PWM from a JASPAR-format file.
 *
 * Expected format:
 * @code
 *   > MA0007.2  ...
 *   A  [ 5  1  3 ... ]
 *   C  [ 2  0  7 ... ]
 *   G  [ 3  9  0 ... ]
 *   T  [ 0  0  0 ... ]
 * @endcode
 *
 * @param filepath Path to .jaspar file.
 * @return Normalised PWM.
 * @throws std::runtime_error on I/O or format errors.
 */
PWM load_jaspar_pwm(const std::string& filepath);

} // namespace motif

#endif // MOTIF_IO_UTILS_HPP
