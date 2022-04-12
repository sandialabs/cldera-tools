module cldera_interface_mod

  use iso_c_binding, only: r8 => c_double
  implicit none

  integer, parameter :: max_str_len = 256

  interface f2c
    module procedure string_f2c, int_f2c, aint_f2c
  end interface f2c

contains

  ! Initialize cldera session and main structures
  subroutine cldera_init (comm)
    use cldera_interface_f2c_mod, only: cldera_init_c
    integer, intent(in) :: comm

    call cldera_init_c(f2c(comm))
  end subroutine cldera_init

  ! Add a partitioned field to cldera data base
  subroutine cldera_add_partitioned_field(name,rank,dims,nparts,part_dim)
    use cldera_interface_f2c_mod, only: capf_c => cldera_add_partitioned_field_c
    character (len=*), intent(in) :: name
    integer, intent(in) :: rank, nparts, part_dim
    integer, intent(in) :: dims(:)

    call capf_c(f2c(name),f2c(rank),f2c(dims),f2c(nparts),f2c(part_dim))
  end subroutine cldera_add_partitioned_field

  ! Set data of a particular field partition in the cldera data base
  subroutine cldera_set_field_partition (name,part,beg,data)
    use cldera_interface_f2c_mod, only: csfp_c => cldera_set_field_partition_c
    character (len=*), intent(in) :: name
    integer,  intent(in) :: part, beg
    real(r8), intent(in), target :: data(:)

    call csfp_c(f2c(name),f2c(part),f2c(beg),data)

  end subroutine cldera_set_field_partition

  ! Check all parts have been set in each field
  subroutine cldera_commit_all_fields ()
    use cldera_interface_f2c_mod, only: cldera_commit_all_fields_c

    call cldera_commit_all_fields_c()
  end subroutine cldera_commit_all_fields

  ! Finalize any pending op (e.g., I/O) and clean up the cldera session
  subroutine cldera_clean_up ()
    use cldera_interface_f2c_mod, only: cldera_clean_up_c

    call cldera_profilng_clean_up_c()
  end subroutine cldera_clean_up

  function string_f2c (f_string) result(c_string)
    use iso_c_binding, only: c_char, c_null_char
    character (len=*), intent(in) :: f_string

    character (kind=c_char,len=max_str_len) :: c_string

    c_string = trim(f_string) // c_null_char
  end function string_f2c

  function int_f2c (int_f) result(int_c)
    use iso_c_binding, only: c_int
    integer, intent(in) :: int_f
    integer(kind=c_int) :: int_c

    int_c = int_f
  end function int_f2c

  function aint_f2c (aint_f) result(aint_c)
    use iso_c_binding, only: c_int
    integer, intent(in) :: aint_f(:)
    integer(kind=c_int), allocatable :: aint_c(:)
    integer :: i,n

    n = size(aint_f)
    allocate (aint_c(n))
    do i=1,n
      aint_c(i) = aint_f(i)
    enddo
  end function aint_f2c

end module cldera_interface_mod
