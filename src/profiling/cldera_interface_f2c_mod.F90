module cldera_interface_f2c_mod

  interface

    ! Initialize cldera session and main structures
    subroutine cldera_init_c (fcomm,ymd,tod) bind(C)
      use iso_c_binding, only: c_int
      integer (kind=c_int), value, intent(in) :: fcomm,ymd,tod
    end subroutine cldera_init_c

    ! Add a partitioned field to cldera data base
    subroutine cldera_add_partitioned_field_c (fname, rank, dims, dimnames, nparts, part_dim,view) bind(c)
      use iso_c_binding, only: c_int, c_char, c_bool, c_ptr
      type(c_ptr), intent(in) :: fname
      integer (kind=c_int), value, intent(in) :: rank,nparts,part_dim
      integer (kind=c_int), intent(in) :: dims(rank)
      logical (kind=c_bool), value, intent(in) :: view
      type(c_ptr), intent(in) :: dimnames(rank)
    end subroutine cldera_add_partitioned_field_c

    ! Set size of a particular field partition in the cldera data base
    subroutine cldera_set_field_part_size_c (fname, part, part_size) bind(c)
      use iso_c_binding, only: c_int, c_char, c_double, c_ptr
      integer (kind=c_int), value, intent(in) :: part, part_size
      type(c_ptr), intent(in) :: fname
    end subroutine cldera_set_field_part_size_c

    ! Set pointer to host app data for a particular field partition
    subroutine cldera_set_field_part_data_c (fname, part, data) bind(c)
      use iso_c_binding, only: c_int, c_char, c_double, c_ptr
      integer (kind=c_int), value, intent(in) :: part
      type(c_ptr), intent(in) :: fname
      type(c_ptr), intent(in) :: data
    end subroutine cldera_set_field_part_data_c

    ! Check all parts have been set in each field
    subroutine cldera_commit_all_fields_c () bind(c)
    end subroutine cldera_commit_all_fields_c

    ! Compute requested stats
    subroutine cldera_compute_stats_c (ymd,tod) bind(c)
      use iso_c_binding, only: c_int
      integer(kind=c_int), intent(in), value :: ymd, tod
    end subroutine cldera_compute_stats_c

    ! Finalize any pending op (e.g., I/O) and clean up the cldera session
    subroutine cldera_clean_up_c () bind(c)
    end subroutine cldera_clean_up_c
  end interface

end module cldera_interface_f2c_mod

