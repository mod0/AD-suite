
! taping --------------------------------------------


      subroutine push_s0(x)
! $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      double precision :: x
! $OpenAD$ END DECLS
        double_tape(double_tape_pointer)=x
        double_tape_pointer=double_tape_pointer+1
      end subroutine 

      subroutine pop_s0(x)
! $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      double precision :: x
! $OpenAD$ END DECLS
        double_tape_pointer=double_tape_pointer-1
        x=double_tape(double_tape_pointer)
      end subroutine

      subroutine push_s1(x)
! $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      double precision :: x(:)
! $OpenAD$ END DECLS
        write(*,*) "In push_s1(x) size(x,1) = ", size(x,1)
        double_tape(double_tape_pointer:double_tape_pointer+size(x,1)-1)=x(:)
        double_tape_pointer=double_tape_pointer+size(x,1)
      end subroutine 

      subroutine pop_s1(x)
! $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      double precision :: x(:)
! $OpenAD$ END DECLS
        double_tape_pointer=double_tape_pointer-size(x,1)
        x(:)=double_tape(double_tape_pointer:double_tape_pointer+size(x,1)-1)
        write(*,*) "In pop_s1(x) doubel_tape_pointer is", double_tape_pointer
      end subroutine

      subroutine push_s2(x)
! $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      double precision :: x(:,:)
! $OpenAD$ END DECLS
        double_tape(double_tape_pointer:double_tape_pointer+(size(x,1)*size(x,2))-1)=reshape(x,(/size(x,1)*size(x,2)/))
        double_tape_pointer=double_tape_pointer+(size(x,1)*size(x,2))
      end subroutine 

      subroutine pop_s2(x)
! $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      double precision :: x(:,:)
! $OpenAD$ END DECLS
        double_tape_pointer=double_tape_pointer-(size(x,1)*size(x,2))
        x(:,:)=reshape(double_tape(double_tape_pointer:),shape(x))
        write(*,*) "In pop_s2(x) doubel_tape_pointer is", double_tape_pointer
      end subroutine

      subroutine apush(x)
! $OpenAD$ INLINE DECLS
      use OAD_tape
      use OAD_active
      implicit none
      type(active) :: x
! $OpenAD$ END DECLS
        double_tape(double_tape_pointer)=x%v
        double_tape_pointer=double_tape_pointer+1
      end subroutine 

      subroutine apop(x)
! $OpenAD$ INLINE DECLS
      use OAD_tape
      use OAD_active
      implicit none
      type(active) :: x
! $OpenAD$ END DECLS
        double_tape_pointer=double_tape_pointer-1
        x%v=double_tape(double_tape_pointer)
      end subroutine

      subroutine push_i_s0(x)
! $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      integer :: x
! $OpenAD$ END DECLS
        integer_tape(integer_tape_pointer)=x
        integer_tape_pointer=integer_tape_pointer+1
      end subroutine 

      subroutine pop_i_s0(x)
! $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      integer :: x
! $OpenAD$ END DECLS
        integer_tape_pointer=integer_tape_pointer-1
        x=integer_tape(integer_tape_pointer)
      end subroutine

      subroutine push_i_s1(x)
! $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      integer :: x(:)
! $OpenAD$ END DECLS
      integer_tape(integer_tape_pointer:integer_tape_pointer+size(x)-1)=x(:)
      integer_tape_pointer=integer_tape_pointer+size(x)
      end subroutine 

      subroutine pop_i_s1(x)
! $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      integer :: x(:)
! $OpenAD$ END DECLS
        integer_tape_pointer=integer_tape_pointer-size(x,1)
        x(:)=integer_tape(integer_tape_pointer:integer_tape_pointer+size(x,1)-1)
      end subroutine

      subroutine push_i_s2(x)
! $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      integer :: x(:,:)
! $OpenAD$ END DECLS
        integer_tape(integer_tape_pointer:integer_tape_pointer+(size(x,1)*size(x,2))-1)=reshape(x,(/size(x,1)*size(x,2)/))
        integer_tape_pointer=integer_tape_pointer+(size(x,1)*size(x,2))
      end subroutine 

      subroutine pop_i_s2(x)
! $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      integer :: x(:,:)
! $OpenAD$ END DECLS
        integer_tape_pointer=integer_tape_pointer-(size(x,1)*size(x,2))
        x(:,:)=reshape(integer_tape(integer_tape_pointer:),shape(x))
      end subroutine

      subroutine push_b(x)
! $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      logical :: x
! $OpenAD$ END DECLS
        logical_tape(logical_tape_pointer)=x
        logical_tape_pointer=logical_tape_pointer+1
      end subroutine 

      subroutine pop_b(x)
! $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      logical :: x
! $OpenAD$ END DECLS
        logical_tape_pointer=logical_tape_pointer-1
        x=logical_tape(logical_tape_pointer)
      end subroutine

      subroutine push_s(s)
! $OpenAD$ INLINE DECLS
        use OAD_tape
        implicit none
        character(*) :: s
! $OpenAD$ END DECLS
        stringlength_tape(stringlength_tape_pointer)=len(s)
        stringlength_tape_pointer=stringlength_tape_pointer+1
        character_tape(character_tape_pointer:character_tape_pointer+len(s))=s(1:len(s))
        character_tape_pointer=character_tape_pointer+len(s)
      end subroutine 

      subroutine pop_s(s)
! $OpenAD$ INLINE DECLS
        use OAD_tape
        implicit none
        character(*) :: s
! $OpenAD$ END DECLS
        stringlength_tape_pointer=stringlength_tape_pointer-1
        character_tape_pointer=character_tape_pointer-stringlength_tape(stringlength_tape_pointer)
        s(1:len(s))=character_tape(character_tape_pointer:character_tape_pointer+stringlength_tape(stringlength_tape_pointer))
      end subroutine

