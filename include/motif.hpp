/**
 * @file motif.hpp
 * @brief Main header - includes all motif discovery components
 * 
 * Single include for all motif discovery functionality.
 */

#ifndef MOTIF_HPP
#define MOTIF_HPP

#include "motif/types.hpp"
#include "motif/kmer_counting.hpp"
#include "motif/enrichment_score.hpp"
#include "motif/pwm_construction.hpp"
#include "motif/pwm_scoring.hpp"
#include "motif/sequence_masking.hpp"
#include "motif/data_generator.hpp"
#include "motif/fasta_reader.hpp"
#include "motif/markov.hpp"
#include "motif/seed_discovery.hpp"
#include "motif/em_refinement.hpp"
#include "motif/io_utils.hpp"

#endif // MOTIF_HPP
