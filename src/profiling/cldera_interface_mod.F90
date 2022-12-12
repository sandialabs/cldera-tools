#include "cldera_config.f"

module cldera_interface_mod

  use iso_c_binding,   only: c_loc, r8 => c_double
  use iso_fortran_env, only: output_unit
  implicit none

  integer, parameter, public :: max_str_len = CLDERA_MAX_NAME_LEN

  integer :: iulog = output_unit
  logical :: masterproc = .false.

  interface f2c
    module procedure string_f2c, int_f2c
  end interface f2c
  interface c2f
    module procedure string_c2f
  end interface c2f

  interface cldera_set_field_part_data
    module procedure cldera_set_field_part_data_real_1d, &
                     cldera_set_field_part_data_real_2d, &
                     cldera_set_field_part_data_real_3d, &
                     cldera_set_field_part_data_int_1d
  end interface cldera_set_field_part_data
contains

  ! Initialize cldera session and main structures
  subroutine cldera_init (comm, ymd, tod)
    use cldera_interface_f2c_mod, only: cldera_init_c
    integer, intent(in) :: comm, ymd, tod

    call cldera_init_c(f2c(comm),ymd,tod)
  end subroutine cldera_init

  subroutine cldera_set_log_unit (log_unit)
    integer, intent(in) :: log_unit

    iulog = log_unit
  end subroutine cldera_set_log_unit

  subroutine cldera_set_masterproc (am_i_master)
    logical, intent(in) :: am_i_master

    masterproc = am_i_master
  end subroutine cldera_set_masterproc

  ! Add a partitioned field to cldera data base
  subroutine cldera_add_partitioned_field(fname,rank,dims,dimnames,nparts,part_dim,view,dtype)
    use iso_c_binding, only: c_char, c_int, c_bool, c_loc, c_ptr
    use cldera_interface_f2c_mod, only: capf_c => cldera_add_partitioned_field_c
    character (len=*), intent(in) :: fname
    character (len=*), intent(in) :: dimnames(:)
    integer, intent(in) :: nparts, part_dim
    integer, intent(in) :: dims(:)
    integer, intent(in) :: rank
    logical, intent(in), optional :: view
    character (len=*), intent(in), optional :: dtype

    type(c_ptr),allocatable :: c_dimnames_ptrs(:)

    integer :: part_dim_c, i
    logical (kind=c_bool) view_c
    integer(kind=c_int), allocatable :: c_dims(:)
    character (kind=c_char, len=max_str_len), target :: fname_c
    character (kind=c_char, len=max_str_len), allocatable, target :: c_dimnames(:)
    character (kind=c_char, len=max_str_len), target :: dtype_c

    if (present(dtype)) then
      dtype_c = f2c(dtype)
    else
      dtype_c = f2c("real")
    endif

    if (masterproc) then
      write(iulog,fmt='(a)',advance='no') "[cldera profiling] Adding field "//trim(fname)//"("
      do i=1,rank-1
        write(iulog,fmt='(a)',advance='no') trim(dimnames(i))//","
      enddo
      write(iulog,fmt='(a)') trim(dimnames(rank))//")"
    endif

    fname_c = f2c(fname)

    allocate(c_dims(rank))
    allocate(c_dimnames(rank))
    allocate(c_dimnames_ptrs(rank))
    ! Flip dims,dimnames arrays, since in C the first is the slowest striding
    do i=1,rank
      c_dims(i) = f2c(dims(rank-i+1))
      c_dimnames(i) = f2c(dimnames(i))
      c_dimnames_ptrs(rank-i+1) = c_loc(c_dimnames(i))
    enddo

    ! We default to fields being views of Model data, but the user can
    ! override this, and register them as hard copies
    if (present(view)) then
      view_c = LOGICAL(view,kind=c_bool)
    else
      view_c = LOGICAL(.true.,kind=c_bool)
    endif
    ! Flip part_dim, since we flipped the dimensions
    ! Also, flip dims, since in C the first is the slowest striding
    call capf_c(c_loc(fname_c),f2c(rank),c_dims,c_dimnames_ptrs, &
                f2c(nparts),f2c(rank-part_dim),view_c,c_loc(dtype_c))
  end subroutine cldera_add_partitioned_field

  ! Set data of a particular field partition in the cldera data base
  subroutine cldera_set_field_part_size (fname,part,part_size)
    use iso_c_binding, only: c_char, c_loc
    use cldera_interface_f2c_mod, only: csfps_c => cldera_set_field_part_size_c
    character (len=*), intent(in) :: fname
    integer,  intent(in) :: part, part_size

    character (kind=c_char, len=max_str_len), target :: fname_c

    fname_c = f2c(fname)
    ! Subtract 1 to part, b/c of C-vs-Fortran index base
    call csfps_c(c_loc(fname_c),f2c(part-1),f2c(part_size))
  end subroutine cldera_set_field_part_size
  subroutine cldera_set_field_part_data_real_1d (fname,part,data)
    use iso_c_binding, only: c_char, c_loc
    use cldera_interface_f2c_mod, only: csfpd_c => cldera_set_field_part_data_c
    character (len=*), intent(in) :: fname
    integer,  intent(in) :: part
    real(r8), intent(in), target :: data(:)

    character (kind=c_char, len=max_str_len), target :: fname_c, dtype_c

    fname_c = f2c(fname)
    dtype_c = f2c("real")

    ! Subtract 1 to part, b/c of C-vs-Fortran index base
    call csfpd_c(c_loc(fname_c),f2c(part-1),c_loc(data(1)),c_loc(dtype_c))
  end subroutine cldera_set_field_part_data_real_1d
  subroutine cldera_set_field_part_data_real_2d (fname,part,data)
    use iso_c_binding, only: c_char, c_loc
    use cldera_interface_f2c_mod, only: csfpd_c => cldera_set_field_part_data_c
    character (len=*), intent(in) :: fname
    integer,  intent(in) :: part
    real(r8), intent(in), target :: data(:,:)

    character (kind=c_char, len=max_str_len), target :: fname_c, dtype_c

    fname_c = f2c(fname)
    dtype_c = f2c("real")
    ! Subtract 1 to part, b/c of C-vs-Fortran index base
    call csfpd_c(c_loc(fname_c),f2c(part-1),c_loc(data(1,1)),c_loc(dtype_c))
  end subroutine cldera_set_field_part_data_real_2d
  subroutine cldera_set_field_part_data_real_3d (fname,part,data)
    use iso_c_binding, only: c_char, c_loc
    use cldera_interface_f2c_mod, only: csfpd_c => cldera_set_field_part_data_c
    character (len=*), intent(in) :: fname
    integer,  intent(in) :: part
    real(r8), intent(in), target :: data(:,:,:)

    character (kind=c_char, len=max_str_len), target :: fname_c, dtype_c

    fname_c = f2c(fname)
    dtype_c = f2c("real")

    ! Subtract 1 to part, b/c of C-vs-Fortran index base
    call csfpd_c(c_loc(fname_c),f2c(part-1),c_loc(data(1,1,1)),c_loc(dtype_c))
  end subroutine cldera_set_field_part_data_real_3d
  subroutine cldera_set_field_part_data_int_1d (fname,part,data)
    use iso_c_binding, only: c_char, c_loc
    use cldera_interface_f2c_mod, only: csfpd_c => cldera_set_field_part_data_c
    character (len=*), intent(in) :: fname
    integer,  intent(in) :: part
    integer, intent(in), target :: data(:)

    character (kind=c_char, len=max_str_len), target :: fname_c, dtype_c

    fname_c = f2c(fname)
    dtype_c = f2c("int")

    ! Subtract 1 to part, b/c of C-vs-Fortran index base
    call csfpd_c(c_loc(fname_c),f2c(part-1),c_loc(data(1)),c_loc(dtype_c))
  end subroutine cldera_set_field_part_data_int_1d

  ! Check all parts have been set in each field
  subroutine cldera_commit_all_fields ()
    use cldera_interface_f2c_mod, only: cldera_commit_all_fields_c

    call cldera_commit_all_fields_c()
  end subroutine cldera_commit_all_fields

  ! Compute all stats
  subroutine cldera_compute_stats (ymd, tod)
    use cldera_interface_f2c_mod, only: cldera_compute_stats_c
    integer, intent(in) :: ymd, tod

    call cldera_compute_stats_c(f2c(ymd),f2c(tod))
  end subroutine cldera_compute_stats

  ! Finalize any pending op (e.g., I/O) and clean up the cldera session
  subroutine cldera_clean_up ()
    use cldera_interface_f2c_mod, only: cldera_clean_up_c

    call cldera_clean_up_c()
  end subroutine cldera_clean_up

!--------------------------------------------------------!
!                   INTENRAL FUNCTIONS                   !
!--------------------------------------------------------!

  function string_f2c (f_string) result(c_string)
    use iso_c_binding, only: c_char, c_null_char
    character (len=*), intent(in) :: f_string

    character (kind=c_char,len=max_str_len) :: c_string

    c_string = trim(f_string) // c_null_char
  end function string_f2c

  function string_c2f (c_string) result(f_string)
    use iso_c_binding, only: c_char, c_null_char
    character (kind=c_char,len=*), intent(in) :: c_string

    character (len=max_str_len) :: f_string
    integer :: len

    len = index(c_string,c_null_char)-1
    f_string = trim(c_string(1:len))
  end function string_c2f

  function int_f2c (int_f) result(int_c)
    use iso_c_binding, only: c_int
    integer, intent(in) :: int_f
    integer(kind=c_int) :: int_c

    int_c = int_f
  end function int_f2c

end module cldera_interface_mod
