#include "cldera_config.f"

module cldera_interface_mod

  use iso_c_binding, only: c_loc, r8 => c_double
  implicit none

  integer, parameter, public :: max_str_len = CLDERA_MAX_NAME_LEN

  interface f2c
    module procedure string_f2c, int_f2c, real_f2c, aint_f2c
  end interface f2c
  interface c2f
    module procedure string_c2f
  end interface c2f

  interface cldera_set_field_part_data
    module procedure cldera_set_field_part_data_1d, &
                     cldera_set_field_part_data_2d, &
                     cldera_set_field_part_data_3d
  end interface cldera_set_field_part_data
contains

  ! Initialize cldera session and main structures
  subroutine cldera_init (comm)
    use cldera_interface_f2c_mod, only: cldera_init_c
    integer, intent(in) :: comm

    call cldera_init_c(f2c(comm))
  end subroutine cldera_init

  function cldera_get_field_names () result(names)
    use iso_c_binding, only: c_char, c_loc
    use cldera_interface_f2c_mod, only: cgnf_c=>cldera_get_num_fields_c, &
                                        cgfn_c=>cldera_get_field_name_c

    character (len=max_str_len), allocatable :: names (:)
    character (kind=c_char, len=max_str_len), target :: fname_c
    integer :: nfields, ifield

    nfields = cgnf_c()
    allocate(names(nfields))

    do ifield=1,nfields
      call cgfn_c (ifield-1, c_loc(fname_c))
      names(ifield) = c2f(fname_c)
    enddo
  end function cldera_get_field_names

  ! Add a partitioned field to cldera data base
  subroutine cldera_add_partitioned_field(fname,dims,dimnames,nparts,part_dim)
    use iso_c_binding, only: c_char, c_loc, c_ptr
    use cldera_interface_f2c_mod, only: capf_c => cldera_add_partitioned_field_c
    character (len=*), intent(in) :: fname
    character (len=*), intent(in) :: dimnames(:)
    integer, intent(in) :: nparts, part_dim
    integer, intent(in) :: dims(:)

    type(c_ptr),allocatable :: c_dimnames_ptrs(:)

    integer :: part_dim_c, rank, i
    character (kind=c_char, len=max_str_len), target :: fname_c
    character (kind=c_char, len=max_str_len), allocatable, target :: c_dimnames(:)

    ! The index of partition in (dim0,...,dimN) needs to be flipped,
    ! since C will read the dims array backwards
    rank = size(dims)
    part_dim_c = rank - part_dim
    fname_c = f2c(fname)

    allocate(c_dimnames(rank))
    allocate(c_dimnames_ptrs(rank))
    do i=1,rank
      c_dimnames(i) = f2c(dimnames(i))
      c_dimnames_ptrs(i) = c_loc(c_dimnames(i))
    enddo

    call capf_c(c_loc(fname_c),f2c(rank),f2c(dims),c_dimnames_ptrs,f2c(nparts),f2c(part_dim_c))
  end subroutine cldera_add_partitioned_field

  ! Set data of a particular field partition in the cldera data base
  subroutine cldera_set_field_part_size (fname,part,part_size)
    use iso_c_binding, only: c_char, c_loc
    use cldera_interface_f2c_mod, only: csfps_c => cldera_set_field_part_size_c
    character (len=*), intent(in) :: fname
    integer,  intent(in) :: part, part_size

    character (kind=c_char, len=max_str_len), target :: fname_c

    fname_c = f2c(fname)
    call csfps_c(c_loc(fname_c),f2c(part-1),f2c(part_size))

  end subroutine cldera_set_field_part_size
  subroutine cldera_set_field_part_data_1d (fname,part,data)
    use iso_c_binding, only: c_char, c_loc
    use cldera_interface_f2c_mod, only: csfpd_c => cldera_set_field_part_data_c
    character (len=*), intent(in) :: fname
    integer,  intent(in) :: part
    real(r8), intent(in), target :: data(:)

    character (kind=c_char, len=max_str_len), target :: fname_c

    fname_c = f2c(fname)
    call csfpd_c(c_loc(fname_c),f2c(part-1),c_loc(data(1)))

  end subroutine cldera_set_field_part_data_1d
  subroutine cldera_set_field_part_data_2d (fname,part,data)
    use iso_c_binding, only: c_char, c_loc
    use cldera_interface_f2c_mod, only: csfpd_c => cldera_set_field_part_data_c
    character (len=*), intent(in) :: fname
    integer,  intent(in) :: part
    real(r8), intent(in), target :: data(:,:)

    character (kind=c_char, len=max_str_len), target :: fname_c

    fname_c = f2c(fname)
    call csfpd_c(c_loc(fname_c),f2c(part-1),c_loc(data(1,1)))
  end subroutine cldera_set_field_part_data_2d
  subroutine cldera_set_field_part_data_3d (fname,part,data)
    use iso_c_binding, only: c_char, c_loc
    use cldera_interface_f2c_mod, only: csfpd_c => cldera_set_field_part_data_c
    character (len=*), intent(in) :: fname
    integer,  intent(in) :: part
    real(r8), intent(in), target :: data(:,:,:)

    character (kind=c_char, len=max_str_len), target :: fname_c

    fname_c = f2c(fname)
    call csfpd_c(c_loc(fname_c),f2c(part-1),c_loc(data(1,1,1)))
  end subroutine cldera_set_field_part_data_3d

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

  function real_f2c (real_f) result(real_c)
    use iso_c_binding, only: c_double
    real, intent(in) :: real_f
    real(kind=c_double) :: real_c

    real_c = real_f
  end function real_f2c

  function aint_f2c (aint_f) result(aint_c)
    use iso_c_binding, only: c_int
    integer, intent(in) :: aint_f(:)
    integer(kind=c_int), allocatable :: aint_c(:)
    integer :: i,n

    n = size(aint_f)
    allocate (aint_c(n))
    do i=1,n
      aint_c(i) = aint_f(n-i+1)
    enddo
  end function aint_f2c

end module cldera_interface_mod
