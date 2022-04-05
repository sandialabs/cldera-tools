#ifndef CLDERA_PROFILING_INTERFACE_HPP
#define CLDERA_PROFILING_INTERFACE_HPP

#include <mpi.h>

#include "cldera_profiling_types.hpp"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * This file contains the declaration of cldera interface routines
 *
 * These routines are meant to be called from Fortran code,
 * which will then not need this header file. However, we still
 * provide this file for the sake of testing, and possibly
 * other C/C++ codes that might use the profiling library,
 * so that they do not need to declare all the external functions.
 */

// Initialize/finalize the cldera profiling session
void cldera_profiling_init (const MPI_Fint& fcomm);
void cldera_profiling_clean_up ();

// Set a field ptr/metadata in the profiling archive
void cldera_add_partitioned_field (
    const char* name,
    const int& rank,
    const int*& dims,
    const int& num_parts,
    const int& part_dim);

// Shortcut in case of single partition
void cldera_add_field (
    const char* name,
    const int& rank,
    const int*& dims);

void cldera_set_field_partition (
    const char* name,
    const int&  part,
    const int&  part_beg,
    const cldera::Real*& data);

// Shortcut in case of single partition
void cldera_set_field (
    const char* name,
    const cldera::Real*& data);

void cldera_commit_fields ();
void cldera_commit_field (const char* name);

// Read from yaml file the list of fields to profile,
// as well as what statistics to be profiled
void cldera_init_requests ();

// Compute all the reuqested stats and store results
// in the profiling archive
void cldera_compute_stats (const cldera::Real& time);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // CLDERA_PROFILING_INTERFACE_HPP
