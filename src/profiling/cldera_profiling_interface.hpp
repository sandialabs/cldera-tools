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
void cldera_init_c (const MPI_Fint fcomm);
void cldera_clean_up_c ();

// Set a field ptr/metadata in the profiling archive
void cldera_add_partitioned_field_c (
    const char*& name,
    const int   rank,
    const int*  dims,
    const int   num_parts,
    const int   part_dim);

// Shortcut in case of single partition
void cldera_add_field_c (
    const char*& name,
    const int   rank,
    const int*  dims);

void cldera_set_field_partition_c (
    const char*& name,
    const int   part,
    const int   part_size,
    const cldera::Real*&  data);

// Shortcut in case of single partition
void cldera_set_field_c (
    const char*& name,
    const cldera::Real*&  data);

void cldera_commit_fields_c ();
void cldera_commit_field_c (const char*& name);

// Read from yaml file the list of fields to profile,
// as well as what statistics to be profiled
void cldera_init_requests_c ();

// Compute all the reuqested stats and store results
// in the profiling archive
void cldera_compute_stats_c (const int ymd, const int tod);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // CLDERA_PROFILING_INTERFACE_HPP
