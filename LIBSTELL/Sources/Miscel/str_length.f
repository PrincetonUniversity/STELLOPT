c-------------------------------------------------------------------/str_length
c This routine determines the length of the string, excluding
c trailing blanks.  (Dick Wieland)
C
C  DMC OCT 1991 -- PIRATED FROM DICK'S STRUTIL.FOR AND PLACED IN
C  THE TRANSP CMSMMS LIBRARY VAXONLY -- TO RESOLVE THE REFERENCES
C  FROM UFLIB ON THE DECSTATION WITHOUT HAVING TO PULL IN ALL OF
C  STRUTIL
 
      integer function str_length(mainstr)
 
      character*(*) mainstr
      integer   mlen
      integer   null
      integer   blank
 
      data null/0/
      data blank/32/
 
c-      fortran 90 can zero length strings
      if (len(mainstr)==0) then
         str_length=0
         return
      end if
c-  if the first character is NUL, assume an empty string;
      if ( ichar(mainstr(1:1)) .eq. 0 )   then
         str_length = 0
         return
      end if
 
c-  start the leftward moving test for blanks or NULs
      mlen = len(mainstr)
      do while ((mlen.gt.1).and.(ichar(mainstr(mlen:mlen)).eq.blank
     >           .or. ichar(mainstr(mlen:mlen)).eq.null))
        mlen = mlen - 1
      end do
 
c-  MLEN=1 a special case:
      if ( mlen .eq. 1 .and. (ichar(mainstr(1:1)) .eq. blank .or.
     >           ichar(mainstr(mlen:mlen)).eq.null)) mlen = 0
 
      str_length = mlen
 
      return
      end
