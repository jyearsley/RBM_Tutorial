
      PROGRAM create_station_list
c
c  This program takes a flow direction file, a maskfile, and an input file 
c  containing paths and a list of station names, lat lons, and basin areas and writes
c  a station list for the routing program.  The program searches an adjustable range of cells
c  around the closest lat lon match and selects the closest basin area match.
c
c  The program will also write out a file with number of cells upstream for checking the 
c  locations. 

c  Note that the simple search only identifies an approximate location, which must be checked.
c
c  written by A. Wood mostly from bits of D. Lohmann and G. O'Donnell's routing code
c  modified by A. Hamlet 4/26/01


      IMPLICIT NONE

      
c     change dimensions here
c     nrow and ncol should be larger than the grid
      INTEGER NROW, NCOL, CNT,FLAG,NLOCS,Z,K,W,M,N
      PARAMETER (NROW = 200, NCOL = 200)
      CHARACTER CHARDUM*14
      REAL VOIDVAL,FACTOR_SUM,XC,YC,SIZE,LAT,LON
      REAL CLAT,CLON
      REAL ESTAREA,CHECK
      INTEGER CROW,CCOL
      INTEGER FPI,FPJ,FNCELLS,ST(200,200)
      INTEGER PMAX,I,J,L,IROW,ICOL,PI,PJ,NO_OF_BOX,DPREC

      PARAMETER (PMAX = 10000)  !Put something larger than your basin cellno here
      
      INTEGER DIREC(NCOL,NROW,2)
      REAL    XMASK(NCOL,NROW), FRACTION(NCOL,NROW)
      REAL    NEWMASK(NCOL,NROW),AREA
      INTEGER CATCHIJ(PMAX,2), H(NCOL,NROW),HH(NCOL,NROW)
      
      CHARACTER*21 NAME
      CHARACTER HEADER(5)*40
      CHARACTER*72 FILE_INPUT, FILENAME, INPATH
      CHARACTER*21 NAME2,NAME3

C***********************************************************
C     OPEN NECESSARY FILES
C***********************************************************
      W=-99
      Z=1


c     open input file
      OPEN(1, FILE='Willamette.inp',STATUS='OLD', ERR=9001)


c     get file name for output row and column list 
      READ(1,*) FILENAME
      OPEN(2,FILE=FILENAME,STATUS='UNKNOWN',ERR=9002)


c     get data from flow direction file
      READ(1,*) FILENAME
      CALL READ_DIREC(DIREC,NCOL,NROW,H,XC,YC,SIZE,
     $     FILENAME,IROW,ICOL)


c     read fraction file
c      READ(1,*) FILENAME
c      CALL READ_FRACTION(FRACTION,NCOL,NROW,FILENAME,
c     $        IROW,ICOL)


c     read number of locations from input file
      READ(1,*) NLOCS


c     loop over number of locations to process

      DO K=1,NLOCS

      READ(1,*) CROW, CCOL

c     search catchment expects row and column in reverse order
c     i.e. position 1,1  is IROW, ICOL
      PJ = IROW + 1 -CROW
      PI = ICOL +1 -CCOL

c     find cells routed to the cell
c     note that this version of search catchment writes a list of 
c     cells to a hard coded file

      CALL SEARCH_CATCHMENT
     &   (PI,PJ,DIREC,NCOL,NROW,
     &   NO_OF_BOX,CATCHIJ,PMAX,IROW,ICOL,NEWMASK,VOIDVAL)

      END DO

      do I = 1,30
        do J = 1,37
          if (newmask(i,j) .gt. 0.0) then
            hh(i,j) = h(i,j)
          else
            hh(i,j) = 0
          end if
        end do
      end do
      do j = 37,1,-1
        write(50,'(50i2)') (hh(i,j),i=30,1,-1)
      end do

      STOP
 9001 WRITE(*,*) 'CANNOT OPEN INPUT FILE: '
 9002 WRITE(*,*) 'CANNOT OPEN OUTPUT FILE ' 
      END

      SUBROUTINE SEARCH_CATCHMENT
     & (PI,PJ,DIREC,NCOL,NROW,NO_OF_BOX,CATCHIJ,PMAX,
     $  IROW,ICOL,NEWMASK,VOIDVAL)

      IMPLICIT NONE

      INTEGER PI,PJ,I,J,NCOL,NROW,PMAX,ICOL,IROW,N
      INTEGER II, JJ, III, JJJ,NO_OF_BOX
      INTEGER DIREC(NCOL,NROW,2)
      INTEGER CATCHIJ(PMAX,2)
      REAL NEWMASK(NCOL,NROW),VOIDVAL

C****** CATCHMENTS ***************************************

      NO_OF_BOX = 0

      DO I = 1, ICOL
         DO J = 1, IROW
            N=0
            NEWMASK(I,J) = VOIDVAL 
            write(*,*) 'i=',i,'j=',j,'no_of_box=', no_of_box
            II = I
            JJ = J


 300        CONTINUE
            N= N+1

            IF(N.GT.7000) THEN
               print*, 'sink encountered'
               print*, 'offending cell', II,JJ
               GOTO 310              
               END IF


            IF ((II .GT. ICOL) .OR. (II .LT.1) .OR. 
     &          (JJ .GT. IROW) .OR. (JJ .LT.1)) THEN
               GOTO 310
            END IF
            IF ((II .EQ. PI) .AND. (JJ .EQ. PJ)) THEN 
               NO_OF_BOX = NO_OF_BOX + 1
              CATCHIJ(NO_OF_BOX,1) = I
               CATCHIJ(NO_OF_BOX,2) = J
               NEWMASK(I,J) = 1.0                 !put current cell in
               GOTO 310                           !the mask
            ELSE IF ((DIREC(II,JJ,1).NE.0) .AND.    !check if the current
     &             (DIREC(II,JJ,2) .NE.0)) THEN   !ii,jj cell routes down
                     III = DIREC(II,JJ,1)         !to the subbasin outlet
                     JJJ = DIREC(II,JJ,2)         !point, following the

            IF(DIREC(III,JJJ,1).EQ.II .AND. DIREC(III,JJJ,2).EQ.JJ) THEN
               print*, 'cell ', 129-II,97-JJ,' directs to itself'
              GOTO 310
               END IF

                     II  = III                    !direction of direc(,)
                     JJ  = JJJ                    !from each cell
                     GOTO 300
             END IF                               !if you get there,
c                                                 !no_of_box increments  
 310        CONTINUE                              !and you try another   
         END DO                                   !cell.
      END DO
 
      WRITE(*,*) 'Number of grid cells upstream of present station',
     $     no_of_box

      RETURN
      END



      SUBROUTINE READ_DIREC(DIREC,NCOL,NROW,H,XC,YC,SIZE
     $     ,FILENAME,IROW,ICOL)

c  reads the flow direction file.

      IMPLICIT NONE

      INTEGER NCOL,NROW,I,J,IROW,ICOL,IMISS
      INTEGER DIREC(NCOL,NROW,2) 
      INTEGER H(NCOL,NROW)
      REAL XC, YC, SIZE
      CHARACTER*72 FILENAME
      CHARACTER*14 CDUM 

      OPEN(10, FILE = FILENAME, STATUS='OLD',ERR=9001)

      READ(10,*) CDUM, ICOL    !note: ARC/INFO style header
      READ(10,*) CDUM, IROW
      READ(10,*) CDUM, XC
      READ(10,*) CDUM, YC
      READ(10,*) CDUM, SIZE
      READ(10,*) CDUM, IMISS
      print*, 'Flow direction file dimensions:'
      print*, 'cols:', icol,' rows:',irow
      print*, 'x-corner:', xc,' ycorner:',yc
      print*, 'size:', size, ' void:', imiss
      IF(IROW.GT.NROW .OR. ICOL.GT.NCOL)THEN
         WRITE(*,*) 'Incorrect dimensions:'
         WRITE(*,*) 'Reset nrow and ncol in main to;',
     $        irow, icol
         STOP
      ENDIF
      
      DO J = IROW,1,-1
         READ(10,*) (H(I,J), I=ICOL,1,-1) 
      END DO      
      CLOSE(10)

      DO I = 1, ICOL
         DO J = 1,IROW
            IF (H(I,J) .EQ. 9) THEN     !note outlet is 9, not 0
               DIREC(I,J,1) = 0
               DIREC(I,J,2) = 0
            ELSE IF (H(I,J) .EQ. 1) THEN
               DIREC(I,J,1) = I
               DIREC(I,J,2) = J+1
            ELSE IF (H(I,J) .EQ. 2) THEN
               DIREC(I,J,1) = I-1
               DIREC(I,J,2) = J+1
            ELSE IF (H(I,J) .EQ. 3) THEN
               DIREC(I,J,1) = I-1
               DIREC(I,J,2) = J
            ELSE IF (H(I,J) .EQ. 4) THEN 
               DIREC(I,J,1) = I-1
               DIREC(I,J,2) = J-1
            ELSE IF (H(I,J) .EQ. 5) THEN
               DIREC(I,J,1) = I
               DIREC(I,J,2) = J-1
            ELSE IF (H(I,J) .EQ. 6) THEN
               DIREC(I,J,1) = I+1
               DIREC(I,J,2) = J-1
            ELSE IF (H(I,J) .EQ. 7) THEN
               DIREC(I,J,1) = I+1
               DIREC(I,J,2) = J
            ELSE IF (H(I,J) .EQ. 8) THEN
               DIREC(I,J,1) = I+1
               DIREC(I,J,2) = J+1
            END IF
         END DO
      END DO
      RETURN
 9001 WRITE(*,*) 'CANNOT OPEN INPUT FILE IN READ_DIREC',
     $  FILENAME
      STOP
      END



