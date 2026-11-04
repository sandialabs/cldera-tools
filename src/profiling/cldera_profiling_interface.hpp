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

void cldera_set_field_part_extent_c (
    const char*& name,
    const int   part,
    const int   part_extent);

void cldera_set_field_part_data_c (
    const char*& name,
    const int   part,
    const void*& data,
    const char*& dtype);

// Like the above, but cldera may WRITE into the host-app memory.
// The field must have been made writable via cldera_set_field_writable_c
// (after add_partitioned_field, before setting any data).
void cldera_set_field_part_data_nonconst_c (
    const char*& name,
    const int   part,
    void*& data,
    const char*& dtype);

void cldera_set_field_writable_c (const char*& name);

// Tell the SAI writer (if any) whether the target field is stored bottom-up
// (EAM alt_data, flipped later by extfrc_set) or top-down (model order).
void cldera_set_sai_storage_bottom_up_c (const bool bottom_up);

// Returns 1 if the context yaml requested an SAI writer, 0 otherwise.
int cldera_has_sai_writer_c ();

void cldera_commit_field_c (const char*& name);

void cldera_commit_all_fields_c ();

void cldera_compute_stats_c (const int ymd, const int tod);

// Run the controller if it's time to (also called internally by compute_stats).
// Result is broadcast to all ranks and stored in the context.
void cldera_compute_controller_c (const int ymd, const int tod);

// Copy the latest controller output (Tg/yr per injection site) into vals.
// On input, *nvals is the capacity of vals; on output, the number of values.
// *has_value is 0 if the controller has not produced any output yet.
void cldera_get_controller_output_c (double* vals, int* nvals, int* has_value);

} // namespace cldera

} // extern "C"

#endif // CLDERA_PROFILING_INTERFACE_HPP
