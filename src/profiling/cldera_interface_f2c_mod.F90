module cldera_interface_f2c_mod

  interface

    ! Initialize cldera session and main structures
    subroutine cldera_init_c (fcomm) bind(C)
      use iso_c_binding, only: c_int
      integer (kind=c_int), value, intent(in) :: fcomm
    end subroutine cldera_init_c

    ! Add a partitioned field to cldera data base
    subroutine cldera_add_partitioned_field_c (fname, rank, dims, nparts, part_dim) bind(c)
      use iso_c_binding, only: c_int, c_char, c_ptr
      type(c_ptr), intent(in) :: fname
      integer (kind=c_int), value, intent(in) :: rank,nparts,part_dim
      integer (kind=c_int), intent(in), dimension(rank) :: dims
    end subroutine cldera_add_partitioned_field_c

    ! Set data of a particular field partition in the cldera data base
    subroutine cldera_set_field_partition_c (fname, part, part_size, data) bind(c)
      use iso_c_binding, only: c_int, c_char, c_double, c_ptr
      integer (kind=c_int), value, intent(in) :: part, part_size
      type(c_ptr), intent(in) :: fname
      type(c_ptr), intent(in) :: data
    end subroutine cldera_set_field_partition_c

    ! Check all parts have been set in each field
    subroutine cldera_commit_all_fields_c () bind(c)
    end subroutine cldera_commit_all_fields_c

    ! Get number of fields to be tracked by clders
    function cldera_get_num_fields_c () bind(c)
      use iso_c_binding, only: c_int
      integer(kind=c_int) :: cldera_get_num_fields_c
    end function cldera_get_num_fields_c

    ! Get name of i-th field to be tracked by cldera
    subroutine cldera_get_field_name_c (i,fname) bind(c)
      use iso_c_binding, only: c_int, c_char, c_ptr
      integer(kind=c_int), intent(in), value :: i
      type(c_ptr), intent(in) :: fname
    end subroutine cldera_get_field_name_c

    ! Finalize any pending op (e.g., I/O) and clean up the cldera session
    subroutine cldera_clean_up_c () bind(c)
    end subroutine cldera_clean_up_c
  end interface

end module cldera_interface_f2c_mod

