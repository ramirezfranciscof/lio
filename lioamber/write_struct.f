      Subroutine write_struct() 
      use garcha_mod

       write(999,*) natom+nsol
       write(999,*)
       do i=1,natom
         write(999,'(I4,F10.6,F10.6,F10.6)') Iz(i),r(i,1:3)*0.529177D0
       enddo
        
       do j=1,nsol
         n=natom+j
         write(999,'(I4,F10.6,F10.6,F10.6)') Pc(n),r(n,1:3)*0.529177D0
       enddo
      
      return
      end subroutine
