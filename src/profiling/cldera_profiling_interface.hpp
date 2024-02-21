#ifndef CLDERA_PROFILING_INTERFACE_HPP
#define CLDERA_PROFILING_INTERFACE_HPP

#include <mpi.h>

extern "C" {

namespace cldera {

void cldera_init_c (const char*& session_name,
                    const MPI_Fint fcomm,
                    const int case_t0_ymd, const int case_t0_tod,
                    const int run_t0_ymd, const int run_t0_tod,
                    const int stop_ymd, const int stop_tod);

void cldera_clean_up_c ();

void cldera_add_field_c (const char*& name,
                         const int    rank,
                         const int*   dims,
                         const char** dimnames,
                         const bool   is_view,
                         const char*& data_type);

void cldera_add_partitioned_field_c (
    const char*&  name,
    const int     rank,
    const int*    dims,
    const char**  dimnames,
    const int     num_parts,
    const int     part_dim,
    const int     part_dim_alloc_size,
    const bool    is_view,
    const char*&  dtype);

void cldera_set_field_part_size_c (
    const char*& name,
    const int   part,
    const int   part_size);

void cldera_set_field_part_data_c (
    const char*& name,
    const int   part,
    const void*& data,
    const char*& dtype);

void cldera_commit_field_c (const char*& name);

void cldera_commit_all_fields_c ();

void cldera_compute_stats_c (const int ymd, const int tod);

} // namespace cldera

} // extern "C"

#endif // CLDERA_PROFILING_INTERFACE_HPP
