module cldera_interface_f2c_mod

  implicit none

  interface

    ! Initialize cldera session and main structures
    subroutine cldera_init_c (fcomm) bind(C)
      use iso_c_binding, only: c_int
      integer (kind=c_int), value, intent(in) :: fcomm
    end subroutine cldera_init_c

    ! Add a partitioned field to cldera data base
    subroutine cldera_add_partitioned_field_c (name, rank, dims, nparts, part_dim) bind(c)
      use iso_c_binding, only: c_int, c_char

      ! type(c_ptr), intent(in) :: name
      character(kind=c_char), intent(in) :: name(*)
      integer (kind=c_int), intent(in) :: rank,nparts,part_dim
      integer (kind=c_int), intent(in), dimension(rank) :: dims
    end subroutine cldera_add_partitioned_field_c

    ! Set data of a particular field partition in the cldera data base
    subroutine cldera_set_field_partition_c (name, part, beg, data) bind(c)
      use iso_c_binding, only: c_int, c_char, c_double
      character (kind=c_char), intent(in) :: name(*)
      integer (kind=c_int),  intent(in) :: part, beg
      real(kind=c_double), intent(in), target :: data(:)
    end subroutine cldera_set_field_partition_c

    ! Check all parts have been set in each field
    subroutine cldera_commit_all_fields_c () bind(c)
    end subroutine cldera_commit_all_fields_c

    ! Finalize any pending op (e.g., I/O) and clean up the cldera session
    subroutine cldera_clean_up_c () bind(c)
    end subroutine cldera_clean_up_c
  end interface

end module cldera_interface_f2c_mod

