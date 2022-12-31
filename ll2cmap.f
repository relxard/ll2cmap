PROGRAM ll2cmap 

!-------------------------------------------------------------------------------
! Program to convert a global {L}atitude-{L}ongitude HYSPLIT formatted 
! meteorology file {to} an HYSPLIT formatted file on a {C}onformal {MAP} 
! projection. This code, LL2CMAP, is distributed under the GNU General Public 
! license https://www.gnu.org/licenses/gpl-3.0.en.html
!-------------------------------------------------------------------------------
! Compilation requires the conformal mapping conversion library: MAP="libcmapf.a"   
! https://www.arl.noaa.gov/wp_arl/wp-content/uploads/utilities/cmap/cmapf.v1_0.tar.gz
! OR
! https://www.iamg.org/documents/oldftp/VOL23/v23-1-5.tar.Z
!-------------------------------------------------------------------------------
! Use the following statements to compile the code with gfortran:
! OPT="-ffree-form -fconvert=big-endian -frecord-marker=4"
! gfortran -oll2cmap ${OPT} ll2cmap.f ${MAP} 
!-------------------------------------------------------------------------------
! Author: roland.draxler@meteozone.com
! 12 Dec 2022 - Initial version
!-------------------------------------------------------------------------------

  IMPLICIT none

  INTEGER                    :: NARG          ! number of command line arguments
  INTEGER                    :: I,J,K         ! temporary index values
  CHARACTER(1)               :: MAP           ! (A)uto|{P}olar|{L}ambert|{M}ercator 
  CHARACTER(1)               :: PROCESS       ! {N}eighbor|(I)nterpolation
  CHARACTER(1)               :: VERTICAL      ! output {I}nput|{S}igma
  CHARACTER(2)               :: CGRID         ! indicator for >999 grid sizes
  CHARACTER(4)               :: KVAR, MODEL   ! variable and model identification
  CHARACTER(50)              :: LABEL         ! label field of data record
  CHARACTER(50)              :: LABELS        ! saved label field for index record
  CHARACTER(5000)            :: HEADER        ! header information for index record
  CHARACTER(256)             :: FINP, FOUT    ! input and output file names
  CHARACTER(256)             :: LABELA        ! label field for argument list
  LOGICAL                    :: FTEST         ! file exist test result
  LOGICAL                    :: DIAG          ! when true diagnostic output on
  LOGICAL                    :: ANAL          ! when true analyze input file content

  REAL                       :: CLAT1,CLON1   ! lower left corner of input grid
  REAL                       :: DLATLON       ! spacing of the input grid
  REAL                       :: CLAT,CLON     ! center of the output grid
  REAL                       :: GRIDS(12)     ! grid definition array

  INTEGER                    :: IYR,IMO,IDA   ! date of input/output and ...
  INTEGER                    :: IHR,IFH       ! hour and forecast hour (<=99)
  INTEGER                    :: ICX           ! extended forecast hour (<=999)

  INTEGER                    :: NX,NY,NZ      ! input grid dimensions
  INTEGER                    :: NXY,LEN,LENH  ! array size, record size, header size
  INTEGER                    :: NDX           ! number of index records
  INTEGER                    :: KSYS          ! 1:sigma 2:pressure 3:z-terrain 4:hybrid

  INTEGER                    :: KREC,KRET     ! input file record counter and EOF test
  INTEGER                    :: NREC          ! output file record counter
  INTEGER                    :: IREC          ! record number of the output index record

! input data structure
  INTEGER, PARAMETER         :: MAXL=256          ! maximum number of levels
  INTEGER, PARAMETER         :: MAXV=64           ! maximum number of variables per level
  REAL                       :: PLVL(MAXL)        ! level heights/pressure/sigma 
  REAL                       :: SLVL(MAXL)        ! optional output sigma level
  REAL                       :: PSFC=1013.0       ! nominal surface pressure for sigma
  REAL                       :: PTOP=10.0         ! nominal model top for sigma 
  INTEGER                    :: NUMVAR(MAXL)      ! number of variables per level
  CHARACTER(4)               :: VARLVL(MAXL,MAXV) ! variable ID by level
  INTEGER                    :: CHKSUM(MAXL,MAXV) ! rotating checksum
  INTEGER                    :: LEV               ! input data level 0=sfc >0=upper
  INTEGER                    :: NEXP              ! packing exponent
  REAL                       :: PREC,VAR1         ! precision and value at 1,1
  REAL,          ALLOCATABLE :: RVARB(:,:)        ! data array one variable and level
  INTEGER(1),    ALLOCATABLE :: KPACK(:)          ! 1-byte packed equivalent of RVARB
  REAL,          ALLOCATABLE :: SMET(:,:,  :)     ! surface data
  REAL,          ALLOCATABLE :: VMET(:,:,:,:)     ! upper level data

! maximum output data structure desired variables
  INTEGER                    :: NVSFC=11
  INTEGER                    :: NVLVL=7
  CHARACTER(4)               :: VCHAR0(100),VCHAR1(100)  ! desired variables
  CHARACTER(4)               :: FCHAR0(100),FCHAR1(100)  ! final matching variables

! output grid and map conversion variables
  INTEGER                    :: NXP,NYP,NXYP  ! output grid and array size
  REAL                       :: GSIZE,GTEMP   ! output grid size
  REAL                       :: PARMAP(9)     ! conformal mapping array
  REAL, ALLOCATABLE          :: TLAT(:,:)     ! latitude of each output grid point
  REAL, ALLOCATABLE          :: TLON(:,:)     ! longitude of each output grid point
  REAL, ALLOCATABLE          :: SDAT(:,:,  :) ! surface output variable 
  REAL, ALLOCATABLE          :: VDAT(:,:,:,:) ! upper level output variable
  INTEGER(1), ALLOCATABLE    :: NPACK(:)      ! 1-byte integer packed output

! variables for CPU internal time and clock functions
  REAL           :: SS
  INTEGER(KIND=8):: DSEC,JULSEC
  INTEGER        :: YR,MO,DA,HR,MN
  INTEGER        :: DATE_TIME(8)
  CHARACTER(12)  :: REAL_CLOCK(3)

! important constants
  REAL, PARAMETER :: REARTH = 6371.2    ! radius of earth in km
  REAL, PARAMETER :: PI     = 3.14159265358979
  REAL, PARAMETER :: DEGPRD = 180.0/PI  ! deg per radian

! set output structure in ARL character format
  DATA VCHAR0 /'PRSS','SHGT','PBLH','TPP6','RH2M','U10M','V10M','T02M',  &
               'UMOF','VMOF','SHTF',89*'    '/
  DATA VCHAR1 /'HGTS','PRES','TEMP','UWND','VWND','WWND','RELH',93*'    '/

!--------------------
  INTERFACE
  SUBROUTINE UNPACK(KPACK,RVARB,NEXP,VAR1,PREC)
  IMPLICIT none
  INTEGER(1),  INTENT(IN)  :: kpack(:)    ! packed input array
  REAL,        INTENT(OUT) :: rvarb(:,:)  ! unpacked output array
  INTEGER,     INTENT(IN)  :: nexp        ! packing exponent
  REAL,        INTENT(IN)  :: var1        ! real array value at 1,1
  REAL,        INTENT(IN)  :: prec        ! precision   
  END SUBROUTINE unpack
  SUBROUTINE MKGRID(MAP,GSIZE,CLAT,CLON,GRIDS,PARMAP,TLAT,TLON)
  IMPLICIT none
  CHARACTER(1), INTENT(INOUT) :: MAP  ! force map projection
  REAL, INTENT(IN)    :: GSIZE        ! grid spacing in km
  REAL, INTENT(IN)    :: CLAT,CLON    ! center of the output grid
  REAL, INTENT(OUT)   :: GRIDS(12)    ! lat-lon input grid values
  REAL, INTENT(OUT)   :: PARMAP(9)    ! conformal map parameters
  REAL, INTENT(OUT)   :: TLAT(:,:)    ! lat-lon of output grid pts
  REAL, INTENT(OUT)   :: TLON(:,:)
  END SUBROUTINE mkgrid
  SUBROUTINE REGRID(V1,V2,TLAT,TLON,CLAT1,CLON1,DLATLON,PROCESS)
  IMPLICIT none
  REAL,         INTENT(IN)  :: V1(:,:)     ! old array on lat-lon grid 
  REAL,         INTENT(OUT) :: V2(:,:)     ! new array on conformal grid
  REAL,         INTENT(IN)  :: TLAT(:,:)   ! lat-lon position of conformal points
  REAL,         INTENT(IN)  :: TLON(:,:)
  REAL,         INTENT(IN)  :: CLAT1,CLON1 ! lat-lon lower left corner point
  REAL,         INTENT(IN)  :: DLATLON     ! lat-lon grid spacing
  CHARACTER(1), INTENT(IN)  :: PROCESS     ! iterpolation process
  END SUBROUTINE regrid
  SUBROUTINE U2GRID(PARMAP,UU,VV)
  IMPLICIT none
  REAL, INTENT(IN)     :: PARMAP(9) ! map conversion constants
  REAL, INTENT(INOUT)  :: UU(:,:)   ! U wind component (S->N)
  REAL, INTENT(INOUT)  :: VV(:,:)   ! V wind component (W->E)
  END SUBROUTINE u2grid
  SUBROUTINE REPACK(RVAR,KPACK,PREC,NEXP,VAR1,KSUM)
  IMPLICIT none
  REAL,       INTENT(IN)  :: RVAR(:,:)   ! data array to be packed
  INTEGER(1), INTENT(OUT) :: KPACK(:)    ! packed int*1 output array
  REAL,       INTENT(OUT) :: PREC        ! precision of packed data array
  INTEGER,    INTENT(OUT) :: NEXP        ! packing scaling exponent
  REAL,       INTENT(OUT) :: VAR1        ! value of real array at position (1,1)
  INTEGER,    INTENT(OUT) :: KSUM        ! rotating checksum
  END SUBROUTINE repack
  SUBROUTINE WINDEX(GRIDS,NREC,LABEL,HEADER,NUMVAR,VARLVL,CHKSUM,PLVL,SLVL,  &
                    NX,NY,NZ,NXY,KSYS,LENH,NDX)
  IMPLICIT none
  REAL,           INTENT(IN)    :: GRIDS(12)        ! grid definition array
  INTEGER,        INTENT(IN)    :: NREC             ! output record number
  CHARACTER(50),  INTENT(IN)    :: LABEL            ! first 50 bytes
  CHARACTER(5000),INTENT(INOUT) :: HEADER           ! header portion of index
  INTEGER,        INTENT(IN)    :: NUMVAR(:)        ! number of variables
  CHARACTER(4),   INTENT(IN)    :: VARLVL(:,:)      ! variable list
  INTEGER,        INTENT(IN)    :: CHKSUM(:,:)      ! rotating checksum  
  REAL,           INTENT(IN)    :: PLVL(:)          ! pressure vertical coordinate
  REAL,           INTENT(IN)    :: SLVL(:)          ! sigma vertical coordinate
  INTEGER,        INTENT(IN)    :: NX,NY,NZ,NXY     ! grid dimensions
  INTEGER,        INTENT(IN)    :: KSYS             ! vertical coordinate system
  INTEGER,        INTENT(IN)    :: LENH             ! header length
  INTEGER,        INTENT(IN)    :: NDX              ! number of index records
  END SUBROUTINE windex
  SUBROUTINE VSIGMA(FCHAR0,FCHAR1,PTOP,PLVL,SLVL,SM,VM)
  IMPLICIT none
  CHARACTER(4), INTENT(IN)    :: FCHAR0(:)   ! surface variables
  CHARACTER(4), INTENT(IN)    :: FCHAR1(:)   ! upper variables 
  REAL,         INTENT(IN)    :: PTOP        ! nominal model top for sigma  
  REAL,         INTENT(IN)    :: PLVL(:)     ! level heights/pressure/sigma
  REAL,         INTENT(IN)    :: SLVL(:)     ! optional output sigma level
  REAL,         INTENT(IN)    :: SM(:,:,:)   ! surface meteorology  
  REAL,         INTENT(INOUT) :: VM(:,:,:,:) ! vertical meteorology  
  END SUBROUTINE vsigma
  END INTERFACE
!--------------------

  CALL DATE_AND_TIME(REAL_CLOCK(1),REAL_CLOCK(2),REAL_CLOCK(3),DATE_TIME)
  READ(REAL_CLOCK(1),'(I4,2I2)')YR,MO,DA   
  READ(REAL_CLOCK(2),'(2I2,F6.3)')HR,MN,SS
  DSEC=HR*3600+MN*60+NINT(SS)
! DSEC=JULSEC(YR,MO,DA,HR,MN,INT(SS))  {function for complete time}

  NARG=IARGC() 
  IF(NARG.EQ.0)THEN
     WRITE(*,*)'Program to convert a global latitude-longitude ARL formatted' 
     WRITE(*,*)'meteorology file to an ARL formatted file on a conformal map' 
     WRITE(*,*)'projection, either polar stereographic, Mercator, or Lambert' 
     WRITE(*,*)'Conformal. The converted file will be a geographic extract of'
     WRITE(*,*)'the orginal, containing the minimum number of variables required'
     WRITE(*,*)'to run HYSPLIT. Additional variables may be added to the output' 
     WRITE(*,*)'as a command line argument. The vertical coordinate of the output'
     WRITE(*,*)'data will be unchanged unless the pressure to sigma interpolation'
     WRITE(*,*)'option is selected, available for pressure level input data.' 
     WRITE(*,*)'The remapping uses bi-linear interpolation of the meteorology' 
     WRITE(*,*)'variables unless the nearest neighbor option is selected.'
     WRITE(*,*)' ' 
     WRITE(*,*)'USAGE: ll2cmap -[options (default)]'
     WRITE(*,*)'  -a[Analyze input file]'
     WRITE(*,*)'  -d[Diagnostic output turned on]'
     WRITE(*,*)'  -i[Input file name]'
     WRITE(*,*)'  -o[Output file name]'
     WRITE(*,*)'  -c[Center lat:lon of output (40.0:-90.0)]'
     WRITE(*,*)'  -g[Grid resolution in km (otherwise from the input file)]'
     WRITE(*,*)'  -n[Number of grid points in lat:lon (25:25)]'
     WRITE(*,*)'  -m[Map projection (A)uto|{P}olar|{L}ambert|{M}ercator]'
     WRITE(*,*)'  -p[Process for remapping {N}eighbor|(I)nterpolation]'
     WRITE(*,*)'  -v[Vertical coordinate for output (I)nput|{S}igma]'
     WRITE(*,*)'  -s[Surface variable to add e.g. PRSS; one for each -s]'
     WRITE(*,*)'  -u[Upper level variable to add e.g. WVEL; one for each -u]'
     STOP      
  END IF

! command line defaults
  CLAT= 40.0
  CLON=-90.0
  GSIZE=-1.0  ! autoset when <0
  NXP=25  
  NYP=25  
  FINP=''
  FOUT=''
  VERTICAL='I'
  PROCESS='I'
  MAP='A'
  DIAG=.FALSE.
  ANAL=.FALSE.

! go through each argument
  DO WHILE (NARG.GT.0)

     CALL GETARG(NARG,LABELA)
     SELECT CASE (LABELA(1:2))

     CASE ('-a','-A')
        ANAL=.TRUE.
        DIAG=.TRUE.
        FOUT='dummy'
     CASE ('-d','-D')
        DIAG=.TRUE.
     CASE ('-i','-I')
        FINP=TRIM(LABELA(3:))
     CASE ('-o','-O')
        FOUT=TRIM(LABELA(3:))
     CASE ('-c','-C')
        K=INDEX(LABELA,':')
        IF(K.NE.0)THEN
           READ(LABELA(3:(K-1)),'(F8.0)')CLAT
           READ(LABELA((K+1):) ,'(F8.0)')CLON
        ELSE
           K=LEN_TRIM(LABELA)
           READ(LABELA(3:K),'(F8.0)')CLAT
        END IF
     CASE ('-g','-G')
        K=LEN_TRIM(LABELA)
        READ(LABELA(3:K),'(F8.0)')GSIZE  
     CASE ('-n','-N')
        K=INDEX(LABELA,':')
        IF(K.NE.0)THEN
           READ(LABELA(3:(K-1)),'(I8)')NYP  
           READ(LABELA((K+1):) ,'(I8)')NXP  
        ELSE
           K=LEN_TRIM(LABELA)
           READ(LABELA(3:K),'(I8)')NXP    
           NYP=NXP
        END IF
     CASE ('-m','-M')
        READ(LABELA(3:),'(A1)')MAP     
        IF(MAP.EQ.'a')MAP='A'
        IF(MAP.EQ.'p')MAP='P'
        IF(MAP.EQ.'l')MAP='L'
        IF(MAP.EQ.'m')MAP='M'
     CASE ('-p','-P')
        READ(LABELA(3:),'(A1)')PROCESS
        IF(PROCESS.EQ.'n') PROCESS='N'
        IF(PROCESS.EQ.'i') PROCESS='I'
     CASE ('-v','-V')
        READ(LABELA(3:),'(A1)')VERTICAL
        IF(MAP.EQ.'i')VERTICAL='I'
        IF(MAP.EQ.'s')VERTICAL='S'
     CASE ('-s','-S')
        READ(LABELA(3:),'(A4)')KVAR 
        NVSFC=NVSFC+1
        VCHAR0(NVSFC)=KVAR
     CASE ('-u','-U')
        READ(LABELA(3:),'(A4)')KVAR 
        NVLVL=NVLVL+1
        VCHAR1(NVLVL)=KVAR
     END SELECT
     NARG=NARG-1
  END DO

! check arguments
  IF(FINP.EQ.''.OR.FOUT.EQ.'')THEN
     WRITE(*,'(2A)')'Both input and output files need to be defined!'
     STOP
  END IF

! open input and output data files
  INQUIRE(FILE=FINP,EXIST=FTEST)
  IF(FTEST)THEN
     IF(DIAG) WRITE(*,'(2A)')'Started processing: ',TRIM(FINP)  
  ELSE
     WRITE(*,'(2A)')'File not found:',FINP  
     STOP
  END IF

! if the output file exists, delete first to avoid not overwriting all
! existing data if the new file contents are less than the existing file
  INQUIRE(FILE=FOUT,EXIST=FTEST)
  IF(FTEST)THEN
     OPEN(20,FILE=FOUT)
     CLOSE(20,STATUS='DELETE')
     IF(DIAG) WRITE(*,'(2A)')'Deleted file      : ',TRIM(FOUT)  
  END IF

! open file to decode the standard label (50) plus the
! fixed portion (108) of the extended header
  OPEN(10,FILE=FINP,RECL=158,ACCESS='DIRECT',FORM='UNFORMATTED')

! decode the standard portion of the index record
  READ(10,REC=1)LABEL,HEADER(1:108)
  READ(LABEL,'(5I2,2X,A2,A4)')IYR,IMO,IDA,IHR,IFH,CGRID,KVAR
  IF(DIAG) WRITE(*,'(A,4I5)')'Opened file       : ',IYR,IMO,IDA,IHR

! decode extended portion of the header
! GRIDS(1)=POLE_LAT  GRIDS(2)=POLE_LON  GRIDS(3)=REF_LAT   GRIDS(4)=REF_LON
! GRIDS(5)=SIZE      GRIDS(6)=ORIENT    GRIDS(7)=TANG_LAT  GRIDS(12)=DUMMY
! GRIDS(8)=SYNC_XP   GRIDS(9)=SYNC_YP   GRIDS(10)=SYNC_LAT GRIDS(11)=SYNC_LON

  READ(HEADER(1:108),'(A4,I3,I2,12F7.0,3I3,I2,I4)') &
      MODEL,ICX,MN,GRIDS,NX,NY,NZ,KSYS,LENH

! check for grid sizes exceeding 999
  I=ICHAR(CGRID(1:1))
  J=ICHAR(CGRID(2:2))
  IF(I.GE.64.OR.J.GE.64)THEN
     NX=(I-64)*1000+NX
     NY=(J-64)*1000+NY
  END IF

! check if vertical transformation requested, only valid with pressure data
  IF(VERTICAL.EQ.'S'.AND.KSYS.NE.2)THEN
     WRITE(*,*)'Only pressure data sets can be transformed to sigma-pressure!'     
     STOP
  END IF

! determine if this is a lat-lon grid
  IF(GRIDS(5).EQ.0.0)THEN
!    determine if the grid is global
     IF((GRIDS(2)+GRIDS(4)-GRIDS(11).EQ.360.0).OR.   &
        (GRIDS(2)+GRIDS(4)-GRIDS(11).EQ.0.0).AND.    &
         GRIDS(1)-GRIDS(10).EQ.180.0)THEN 
         CLAT1=GRIDS(10)
         CLON1=GRIDS(11)
         DLATLON=GRIDS(4)
     ELSE
         WRITE(*,'(A)')'Regional lat-lon grids are not supported!'
         CLOSE(10)
         STOP
     END IF   

!    latitude-longitude grid spacing at selection point
     IF(GSIZE.LT.0.0)THEN
!       auto set output comparable to input file
        GSIZE=FLOAT(NINT(REARTH*GRIDS(3)/DEGPRD))
     ELSE
!       user set output
        GTEMP=FLOAT(NINT(REARTH*GRIDS(3)/DEGPRD))
        IF(PROCESS.EQ.'N')THEN
           IF(GSIZE.LT.GTEMP/2.0)THEN
              WRITE(*,'(A)')'WARNING: an output grid size less than half the input'
              WRITE(*,'(A)')'grid size is not recommended when using the nearest'
              WRITE(*,'(A)')'neighbor interpolation method, use bilinear instead!'
           END IF
        END IF
     END IF
  ELSE
     WRITE(*,'(A)')'The input grid is NOT lat-lon!'
     WRITE(*,'(A,A10  )')'Model   : ',model    
     WRITE(*,'(A,I10  )')'Forecast: ',icx     
     WRITE(*,'(A,I10  )')'Minutes : ',mn        
     WRITE(*,'(A,F10.4)')'Pole_lat: ',grids(1)
     WRITE(*,'(A,F10.4)')'Pole_lon: ',grids(2)
     WRITE(*,'(A,F10.4)')'Ref_lat : ',grids(3)
     WRITE(*,'(A,F10.4)')'Ref_lon : ',grids(4)
     WRITE(*,'(A,F10.4)')'Size    : ',grids(5)
     WRITE(*,'(A,F10.4)')'Orient  : ',grids(6)
     WRITE(*,'(A,F10.4)')'Tang_lat: ',grids(7)
     WRITE(*,'(A,F10.4)')'Sync_xp : ',grids(8)
     WRITE(*,'(A,F10.4)')'Sync_yp : ',grids(9)
     WRITE(*,'(A,F10.4)')'Sync_lat: ',grids(10)
     WRITE(*,'(A,F10.4)')'Sync_lon: ',grids(11)
     WRITE(*,'(A,F10.4)')'Dummy   : ',grids(12)
     WRITE(*,'(A,3I5  )')'Grid    : ',nx,ny,nz  
     WRITE(*,'(A, I5  )')'Vert_sys: ',ksys      
     WRITE(*,'(A, I5  )')'Len_head: ',lenh      
     WRITE(*,'(A, I5  )')'Num indx: ',LENH/(nx*ny)+1 
     CLOSE(10)
     STOP
  END IF

! close file and reopen with proper length
  CLOSE (10)
  NXY = NX*NY
  LEN = NXY+50
  NDX = LENH/NXY+1 
  OPEN(10,FILE=FINP,RECL=LEN,ACCESS='DIRECT',FORM='UNFORMATTED')

! input file summary    
  IF(DIAG)THEN
     WRITE(*,'(A)')        '  '
     WRITE(*,'(A)')        'Input file structure ...'
     WRITE(*,'(A,2I6,2I8)')'Grid size and lrec: ',NX,NY,NXY,LEN
     WRITE(*,'(A,I6)')     'Header record size: ',LENH
     WRITE(*,'(A,I6)')     'Number index recds: ',NDX  
     WRITE(*,'(A,2F6.1)')  'Lower left corner : ',CLAT1,CLON1 
     WRITE(*,'(A,2F6.1)')  'Upper right corner: ',GRIDS(1),GRIDS(2)
     WRITE(*,'(A, F6.1)')  'Lat-Lon spacing   : ',DLATLON     
  END IF

! decode the header record to determine input variable structure
  IF(NDX.EQ.1)THEN
     READ(10,REC=1)LABEL,HEADER(1:LENH)
  ELSE
!    read extended character string over multiple index records
!    which is only done the first time from record 1 to NDX
     I=1
     DO K=1,NDX
        J=I+NXY-1
        IF(K.EQ.NDX)J=I+(LENH-(NDX-1)*NXY)-1
        READ(10,REC=K)LABEL,HEADER(I:J)
        I=J+1       
     END DO
  END IF

! determine the number of surface and upper level variables
! and the character ID at each level
  I=109
  VARLVL='    '
  SLVL=0.0              !  when all zero indicates pressure output
  IF(NZ.GT.MAXL)THEN
     WRITE(*,*)'Number of levels exceeds compiled maximum: ',NZ,MAXL
     STOP
  END IF

  DO K=1,NZ
     READ(HEADER(I:I+7),'(F6.2,I2)') PLVL(K),NUMVAR(K)
     IF(NUMVAR(K).GT.MAXV)THEN
        WRITE(*,*)'Number of variables exceeds compiled maximum: ',NUMVAR(K),MAXV
        WRITE(*,*)'At level: ',K
        STOP
     END IF

     I=I+8
     DO J=1,NUMVAR(K)
        READ(HEADER(I:I+7),'(A4,I3)') VARLVL(K,J),CHKSUM(K,J)
        I=I+8
     END DO

     IF(VARLVL(K,J).EQ.'DIFF')THEN
        WRITE(*,*)'Extended precision DIFF variable not supported by this version!'
        STOP
     END IF

     IF(VERTICAL.EQ.'S')THEN
        IF(K.EQ.1)THEN
           SLVL(K)=1.0
        ELSEIF(K.EQ.NZ)THEN
!          top sigma level below top pressure level for interpolation
           SLVL(K)=(0.5*(PLVL(K)+PLVL(K-1))-PTOP)/(PSFC-PTOP)
        ELSE
           SLVL(K)=(PLVL(K)-PTOP)/(PSFC-PTOP)
        END IF
     END IF

     IF(DIAG)THEN
        IF(K.EQ.1)THEN
           WRITE(*,'(A,I3,F10.3,32(1X,A4))') &
          'Level,Press,Variab: ',K,PLVL(K),(VARLVL(K,1:NUMVAR(K)))
        ELSE
           WRITE(*,'(A,I3,F10.3,32(1X,A4))') &
          '                    ',K,PLVL(K),(VARLVL(K,1:NUMVAR(K)))
        END IF
     END IF
  END DO

! configure the output grid
  IF(NXP.GE.1000.OR.NYP.GE.1000)THEN
!    header record does not support grids of more than 999, the grid 
!    number is converted to character to represent the 1000s digit
!    e.g. @(64)=<1000, A(65)=1000, B(66)=2000, etc
     I=NXP/1000
     J=NYP/1000
     WRITE(CGRID,'(A2)') CHAR(I+64)//CHAR(J+64)
  ELSE
     WRITE(CGRID,'(A2)') '99' 
  END IF
  ALLOCATE (TLAT(NXP,NYP),TLON(NXP,NYP))
  CALL MKGRID(MAP,GSIZE,CLAT,CLON,GRIDS,PARMAP,TLAT,TLON)

! how many surface (1) variables are matched in the output selection 
! fill list (FCHAR0) of final output variables
  K=0
  DO I=1,NVSFC       ! selection list of output variables
  DO J=1,NUMVAR(1)   ! list of input variables
     IF(VCHAR0(I).EQ.VARLVL(1,J)) THEN 
        K=K+1     
        FCHAR0(K)=VCHAR0(I)       
     END IF
  END DO
  END DO
  NUMVAR(1)=K        

! update values in the VARLVL array for output in the new index record
  DO K=1,NUMVAR(1)
     VARLVL(1,K)=FCHAR0(K)
  END DO

! how many upper (2) variables are matched in the output selection 
! fill list (FCHAR1) of final output variables
  K=0
  DO I=1,NVLVL      ! selection list of output variables
  DO J=1,NUMVAR(2)  ! list of input variables
     IF(VCHAR1(I).EQ.VARLVL(2,J)) THEN
        K=K+1            
        FCHAR1(K)=VCHAR1(I)       
     END IF
  END DO
  END DO
  NUMVAR(2:)=K      ! upper level will repeat at all levels

! update values in the VARLVL array for output in the new index record
  DO J=2,NZ
  DO K=1,NUMVAR(2)
     VARLVL(J,K)=FCHAR1(K)
  END DO
  END DO

  IF(DIAG)THEN
     WRITE(*,'(A)')'  '
     WRITE(*,'(A)')              'Output file structure ...'
     WRITE(*,'(A,I3,32(1X,A4))') 'Numsfc, Vars: ',NUMVAR(1),FCHAR0(1:NUMVAR(1))
     WRITE(*,'(A,I3,32(1X,A4))') 'Numlvl, Vars: ',NUMVAR(2),FCHAR1(1:NUMVAR(2))
  END IF

! determine the output file structure, replace certain input file values
! with their output file equivalents
  LENH=108+8*(NZ+NUMVAR(1)+NUMVAR(2)*(NZ-1))     
  NXYP=NXP*NYP
  ALLOCATE (NPACK(NXYP))
  LEN=NXYP+50
  NDX=LENH/NXYP+1 

! check to insure the output record length is sufficient
  IF(LENH.LE.512)THEN
     WRITE(*,'(A)')'Output grid size too small to save index record data!'
     WRITE(*,'(A,I6)')     'Header record size: ',LENH
     WRITE(*,'(A,2I6,2I8)')'Output grid points: ',NXP,NYP
     NXP = NINT(SQRT(FLOAT(LENH)))+1
     WRITE(*,'(A,2I6,2I8)')'Suggested minimum : ',NXP,NXP
     CLOSE(10)
     STOP
  END IF

  IF(DIAG)THEN
     WRITE(*,'(A,2I6,I8)') 'Output nx,ny,lrec : ',NXP,NYP,LEN
     WRITE(*,'(A,I6)')     'Header record size: ',LENH
     WRITE(*,'(A,I6)')     'Number index recds: ',NDX   
     WRITE(*,'(A,2F6.1)')  'Lower left corner : ',TLAT(1,1),TLON(1,1)
     WRITE(*,'(A,2F6.1)')  'Upper right corner: ',TLAT(NXP,NYP),TLON(NXP,NYP)
     WRITE(*,'(A,2F6.1)')  'Grid spacing (km) : ',GSIZE       
     WRITE(*,'(2A  )')     'Interpolation     : ',PROCESS
     IF(ANAL) STOP
  END IF

  OPEN(20,FILE=FOUT,RECL=LEN,ACCESS='DIRECT',FORM='UNFORMATTED')

! allocate array space for input array but only for variables to be output
  ALLOCATE (RVARB(NX,NY))   
  ALLOCATE (KPACK(NXY))
  ALLOCATE (SMET(NX,NY,NUMVAR(1)))
  ALLOCATE (VMET(NX,NY,NZ-1,NUMVAR(2)))

! output array
  ALLOCATE (SDAT(NXP,NYP,NUMVAR(1)),VDAT(NXP,NYP,NZ-1,NUMVAR(2)))
  SDAT=0.0
  VDAT=0.0

! start data processing section          
  KREC=1  ! input record counter
  NREC=1  ! output record counter
  KRET=0

  IF(DIAG) WRITE(*,'(A)')'  '
  DO WHILE (kret.EQ.0)
     READ(10,REC=KREC,IOSTAT=kret)LABEL,(KPACK(K),K=1,NXY)

     IF(kret.EQ.0)THEN
        READ(LABEL,'(6I2,2X,A4,I4,2E14.7)') IYR,IMO,IDA,IHR,IFH,LEV,KVAR,NEXP,PREC,VAR1

        IF(KVAR.NE.'INDX') THEN
!          not index means it is a data record, find the variable match, if any,
!          in the output array to determine the index number of the output variable
!          and then unpack and save the data in the surface or upper data arrays
           IF(LEV.EQ.0)THEN
              DO J=1,NUMVAR(1)
                 IF(FCHAR0(J).EQ.KVAR)    & 
                    CALL UNPACK(KPACK,SMET(:,:,J),NEXP,VAR1,PREC)
              END DO
           ELSE
              DO J=1,NUMVAR(2)
                 IF(FCHAR1(J).EQ.KVAR)    &
                    CALL UNPACK(KPACK,VMET(:,:,LEV,J),NEXP,VAR1,PREC)
              END DO
           END IF
        END IF
     END IF

     IF(kret.NE.0.OR.KVAR.EQ.'INDX')THEN
!       after the first index record, when the second and subsequent index records
!       are read, the arrays are filled with the previous time's data, output these
!       data before processing this index record, which starts the next time period!

        IF(KREC.GT.1)THEN
!          the second index record (time period #2) and subsequent index records
!          output the data from the previous time period
           IF(DIAG)WRITE(*,'(2A,F10.1)')'Output time: ',LABELS(1:10),SMET(NX/2,NY/2,1)

!          update the label and header fields for the output grid
           WRITE(LABELS(13:14),'(A2)') CGRID

!          write the index record with the updated variable fields in the header
           CHKSUM=0
           IREC=NREC
           CALL WINDEX(GRIDS,IREC,LABELS,HEADER,NUMVAR,VARLVL,CHKSUM,PLVL,SLVL,   & 
                       NXP,NYP,NZ,NXYP,KSYS,LENH,NDX)

           NREC=NREC+NDX     ! each variable and level for output to one record
           LABEL=LABELS      ! set the base variables of the label field

!          data array filled with previous time period when new index record is read
!          bilinear interpolate from lat-lon to conformal grid
           DO J=1,NUMVAR(1)
              CALL REGRID(SMET(:,:,J),SDAT(:,:,J),      &
                   TLAT,TLON,CLAT1,CLON1,DLATLON,PROCESS)
           END DO
           DO K=1,(NZ-1)
           DO J=1,NUMVAR(2)
              CALL REGRID(VMET(:,:,K,J),VDAT(:,:,K,J),  &
                   TLAT,TLON,CLAT1,CLON1,DLATLON,PROCESS)
           END DO
           END DO

!          find all velocity vectors and remap to the conformal grid
!          assume that the U and V vectors always follow each other
           DO J=1,NUMVAR(1)
              IF(FCHAR0(J).EQ.'U10M')CALL U2GRID(PARMAP,SDAT(:,:,J),SDAT(:,:,J+1))
           END DO
           DO K=1,(NZ-1)
           DO J=1,NUMVAR(2) 
              IF(FCHAR1(J).EQ.'UWND')CALL U2GRID(PARMAP,VDAT(:,:,K,J),VDAT(:,:,K,J+1))
           END DO
           END DO

!          optional section to remap from pressure to sigma coordinates
           IF(VERTICAL.EQ.'S')CALL VSIGMA(FCHAR0,FCHAR1,PTOP,PLVL,SLVL,SDAT,VDAT)

!          pack the real data into one-byte integers and output to file
           WRITE(LABEL(11:12),'(A2)')' 0'       ! surface level = 0
           DO J=1,NUMVAR(1)
              CALL REPACK(SDAT(:,:,J),NPACK,PREC,NEXP,VAR1,CHKSUM(1,J))
              WRITE(LABEL(15:),'(A4,I4,2E14.7)')FCHAR0(J),NEXP,PREC,VAR1
              WRITE(20,REC=NREC)LABEL,NPACK
              NREC=NREC+1
           END DO
           SDAT=0.0

           DO K=1,(NZ-1)
              WRITE(LABEL(11:12),'(I2)') K      ! NZ includes all levels + sfc
              DO J=1,NUMVAR(2)
                 CALL REPACK(VDAT(:,:,K,J),NPACK,PREC,NEXP,VAR1,CHKSUM(K+1,J))
                 WRITE(LABEL(15:),'(A4,I4,2E14.7)')FCHAR1(J),NEXP,PREC,VAR1
                 WRITE(20,REC=NREC)LABEL,NPACK
                 NREC=NREC+1
              END DO
           END DO
           VDAT=0.0

!          overwrite the index record with the correct checksum values, note that
!          the following call can be disabled without penalty to skip updating
!          the index record with the correct checksum values
           CALL WINDEX(GRIDS,IREC,LABELS,HEADER,NUMVAR,VARLVL,CHKSUM,PLVL,SLVL,   & 
                       NXP,NYP,NZ,NXYP,KSYS,LENH,NDX)
        END IF   

!       at the beginning of each time period, which is indicated by the INDX
!       record, reread the label field and update the time information from the 
!       header, which are characters,previously input as integers
        IF(kret.EQ.0) READ(10,REC=KREC)LABELS,HEADER(1:9)

     END IF     
     KREC=KREC+1
  END DO      

  CLOSE(10)
  CLOSE(20)

  IF(DIAG)THEN
     CALL DATE_AND_TIME(REAL_CLOCK(1),REAL_CLOCK(2),REAL_CLOCK(3),DATE_TIME)
     READ(REAL_CLOCK(1),'(I4,2I2)')YR,MO,DA   
     READ(REAL_CLOCK(2),'(2I2,F6.3)')HR,MN,SS
     DSEC=HR*3600+MN*60+NINT(SS)-DSEC
!    DSEC=JULSEC(YR,MO,DA,HR,MN,INT(SS))-DSEC    {function for complete time}
     WRITE(*,'(A,I10)')'Elapsed time (s) = ',DSEC
  END IF
END PROGRAM ll2cmap   



!-------------------------------------------------------------------------------
! unpacks the one-dimensional 1-byte integer input array and converts it to the
! two-dimensional real*4 output array

SUBROUTINE UNPACK(kpack,rvarb,nexp,var1,prec)

  IMPLICIT none

  INTEGER(1),  INTENT(IN)  :: kpack(:)    ! packed input array 
  REAL,        INTENT(OUT) :: rvarb(:,:)  ! unpacked output array
  INTEGER,     INTENT(IN)  :: nexp        ! packing exponent
  REAL,        INTENT(IN)  :: var1        ! real array value at 1,1
  REAL,        INTENT(IN)  :: prec        ! precision 

  REAL      :: vold,scale
  INTEGER   :: i,j,k,nx,ny

  NX=SIZE(rvarb,1)
  NY=SIZE(rvarb,2)

! unpacked value at RVARB(1,1) should be identical to the
! value of VAR1 from the header record, therefore the
! difference encoded in KDATA(1) should be zero which 
! applies to the first difference in both I and J

! sum row differences for column one
  VOLD=0.0   
  K=1-NX
  DO J=1,NY
     K=K+NX       
     RVARB(1,J)=(IAND(INT(KPACK(K)),255)-127.0)+VOLD
     VOLD=RVARB(1,J)
  END DO

! sum the column differences for each row
  K=0
  DO J=1,NY
     VOLD=RVARB(1,J)
     K=K+1        
     DO I=2,NX
        K=K+1       
        RVARB(I,J)=(IAND(INT(KPACK(K)),255)-127.0)+VOLD
        VOLD=RVARB(I,J)
     END DO
  END DO

! apply the scaling factor and add the differences
! to the base value
  SCALE=2.0**(7-NEXP)
  RVARB=RVARB/SCALE+VAR1

  DO J=1,NY
  DO I=1,NX
     IF(ABS(RVARB(I,J)).LT.PREC) RVARB(I,J)=0.0
  END DO
  END DO

END SUBROUTINE unpack



!-------------------------------------------------------------------------------
! sets the values for the output grid, including the map transformation
! parameters (parmap), from the options set on the commmand line

SUBROUTINE MKGRID(MAP,GSIZE,CLAT,CLON,GRIDS,PARMAP,TLAT,TLON)  

  IMPLICIT none

  CHARACTER(1), INTENT(INOUT) :: MAP  ! force map projection 
  REAL, INTENT(IN)    :: GSIZE        ! grid spacing in km
  REAL, INTENT(IN)    :: CLAT,CLON    ! center of the output grid
  REAL, INTENT(OUT)   :: GRIDS(12)    ! lat-lon input grid values
  REAL, INTENT(OUT)   :: PARMAP(9)    ! conformal map parameters
  REAL, INTENT(OUT)   :: TLAT(:,:)    ! lat-lon of output grid pts
  REAL, INTENT(OUT)   :: TLON(:,:)

  INTEGER             :: I,J,NXP,NYP

  NXP=SIZE(TLAT,1)
  NYP=SIZE(TLAT,2)

  IF(MAP.EQ.'A')THEN
!    auto set map projection based upon latitude
     IF(ABS(CLAT).GT.60.0)THEN
        MAP='P'
     ELSEIF(ABS(CLAT).LT.30.0)THEN
        MAP='M'
     ELSE
        MAP='L'
     END IF
  END IF
  
! defines a polar sterographic projection
  IF(MAP.EQ.'P')THEN
     GRIDS(1)=SIGN(90.0,CLAT) ! set the pole position at +90 or -90
     GRIDS(2)=CLON            ! pole longtitude (+180 from cut)
     GRIDS(3)=CLAT            ! reference lat/lon (at which grid size specified)
     GRIDS(4)=CLON            ! reference longitude and grid alignment
     GRIDS(7)=GRIDS(1)        ! tangent latitude

! defines a mercator projection
  ELSEIF(MAP.EQ.'M')THEN
     GRIDS(1)=0.0
     GRIDS(2)=CLON            ! pole lat/lon axis through pole
     GRIDS(3)=CLAT            ! reference lat
     GRIDS(4)=CLON            ! reference lon
     GRIDS(7)=0.0             ! tangent latitude

! defines a lambert conformal projection
  ELSE
     GRIDS(1)=CLAT
     GRIDS(2)=CLON            ! pole lat/lon axis through pole
     GRIDS(3)=CLAT            ! reference lat
     GRIDS(4)=CLON            ! reference lon
     GRIDS(7)=CLAT            ! tangent latitude
  END IF

  GRIDS(6)=0.0                ! grid orientation
  GRIDS(5)=GSIZE              ! delta=x grid size in km
  GRIDS(8)=(NXP+1.0)/2.0
  GRIDS(9)=(NYP+1.0)/2.0      ! synch point in x,y coordintes
  GRIDS(10)=CLAT
  GRIDS(11)=CLON              ! synch point in lat/lon coordinates
  GRIDS(12)=0.0               ! variable reserved for future use

! define the tangent latitude and reference longitude
  CALL STLMBR(PARMAP,GRIDS(7),GRIDS(4))

! define the grid by a one-point specification
  CALL STCM1P(PARMAP,GRIDS(8),GRIDS(9),GRIDS(10),GRIDS(11),                &
                     GRIDS(3),GRIDS(4),GRIDS(5),GRIDS(6))

! determine the lat/lon at the output grid locations
  DO I=1,NXP
  DO J=1,NYP
     CALL CXY2LL(PARMAP,FLOAT(I),FLOAT(J),TLAT(I,J),TLON(I,J))
!    cxy2ll returns -180:+180 while input data array from 0:360
     IF(TLON(I,J).LT.0.0)TLON(I,J)=360.0+TLON(I,J)
  END DO
  END DO

END SUBROUTINE mkgrid



!-------------------------------------------------------------------------------
! rotates the winds from earth coordinates on the latitude-longitude grid to
! coordinates relative to the conformal map grid

SUBROUTINE U2GRID(PARMAP,UU,VV)

  IMPLICIT none

  REAL, INTENT(IN)     :: PARMAP(9) ! map conversion constants
  REAL, INTENT(INOUT)  :: UU(:,:)   ! U wind component (S->N)
  REAL, INTENT(INOUT)  :: VV(:,:)   ! V wind component (W->E)

  REAL                 :: UG,VG
  INTEGER              :: I,J,NXP,NYP

  NXP=SIZE(UU,1)
  NYP=SIZE(UU,2)

! convert compass winds to grid-orientation
  DO I=1,NXP
  DO J=1,NYP
     CALL CC2GXY(PARMAP,FLOAT(I),FLOAT(J),UU(I,J),VV(I,J),UG,VG)
     UU(I,J)=UG
     VV(I,J)=VG
  END DO
  END DO

END SUBROUTINE u2grid


!-------------------------------------------------------------------------------
! transfers the data on the latitude-longitude grid of array V1 to the conformal
! map grid of array V2 using bi-linear interpolation or the nearest neighbor

SUBROUTINE REGRID(V1,V2,TLAT,TLON,CLAT1,CLON1,DLATLON,PROCESS)

  IMPLICIT none

  REAL,         INTENT(IN)  :: V1(:,:)     ! old array on lat-lon grid 
  REAL,         INTENT(OUT) :: V2(:,:)     ! new array on conformal grid
  REAL,         INTENT(IN)  :: TLAT(:,:)   ! lat-lon position of conformal points
  REAL,         INTENT(IN)  :: TLON(:,:)
  REAL,         INTENT(IN)  :: CLAT1,CLON1 ! lat-lon lower left corner point
  REAL,         INTENT(IN)  :: DLATLON     ! lat-lon grid spacing
  CHARACTER(1), INTENT(IN)  :: PROCESS     ! iterpolation process

  INTEGER             :: NX1,NY1,NX2,NY2
  INTEGER             :: I,J,ILO,JLO,IHI,JHI
  REAL                :: XP,YP,FXI,FYJ,TOP,BOT

  NX1=SIZE(V1,1)
  NY1=SIZE(V1,2)
  NX2=SIZE(V2,1)
  NY2=SIZE(V2,2)

! interpolate values to new grid
  DO I=1,NX2
  DO J=1,NY2

!    compute adjacent index values on grid 1
     XP=1.0+(TLON(I,J)-CLON1)/DLATLON
     YP=1.0+(TLAT(I,J)-CLAT1)/DLATLON

     IF(PROCESS.EQ.'N')THEN
!       nearest neighbor assigment
        IHI=NINT(XP)
        JHI=NINT(YP)

!       check limits, note global grids wrap at the prime meridian
        IF(IHI.GT.NX1)IHI=1
        IF(JHI.GT.NY1)JHI=NY1-1
        IF(JHI.LT.1  )JHI=1

!       assign output point to nearest input point
        V2(I,J)=V1(IHI,JHI)          

     ELSE
!       compute base index
        ILO=INT(XP)
        JLO=INT(YP)

!       interpolation fractions from base point
        FXI=XP-ILO
        FYJ=YP-JLO

!       compute upper index point
        IHI=ILO+1
        JHI=JLO+1

!       check limits, note global grids wrap at the prime meridian
        IF(IHI.GT.NX1)IHI=1
        IF(JHI.GT.NY1)JHI=NY1-1
        IF(JHI.LT.1  )JHI=1

!       interpolate across at top and bottom
        TOP=(V1(IHI,JHI)-V1(ILO,JHI))*FXI+V1(ILO,JHI)
        BOT=(V1(IHI,JLO)-V1(ILO,JLO))*FXI+V1(ILO,JLO)

!       interpolate between top and bottom
        V2(I,J)=(TOP-BOT)*FYJ+BOT
     END IF

  END DO
  END DO

END SUBROUTINE regrid


!-------------------------------------------------------------------------------
! writes the header contents of the index record, which is always the first 
! record for a time period, which preceeds all the data records for that time

SUBROUTINE WINDEX(GRIDS,NREC,LABEL,HEADER,NUMVAR,VARLVL,CHKSUM,PLVL,SLVL,  &
                  NX,NY,NZ,NXY,KSYS,LENH,NDX)

  IMPLICIT none

  REAL,           INTENT(IN)    :: GRIDS(12)        ! grid definition array
  INTEGER,        INTENT(IN)    :: NREC             ! output record number
  CHARACTER(50),  INTENT(IN)    :: LABEL            ! first 50 bytes
  CHARACTER(5000),INTENT(INOUT) :: HEADER           ! header portion of index
  INTEGER,        INTENT(IN)    :: NUMVAR(:)        ! number of variables
  CHARACTER(4),   INTENT(IN)    :: VARLVL(:,:)      ! variable list
  INTEGER,        INTENT(IN)    :: CHKSUM(:,:)      ! rotating checksum  
  REAL,           INTENT(IN)    :: PLVL(:)          ! pressure vertical coordinate
  REAL,           INTENT(IN)    :: SLVL(:)          ! sigma vertical coordinate
  INTEGER,        INTENT(IN)    :: NX,NY,NZ,NXY     ! grid dimensions
  INTEGER,        INTENT(IN)    :: KSYS             ! vertical coordinate system
  INTEGER,        INTENT(IN)    :: LENH             ! header length
  INTEGER,        INTENT(IN)    :: NDX              ! number of index records

  INTEGER                       :: I,J,K            ! temporary indicies
  INTEGER                       :: KVER             ! vertical coordinate 1=sigma 2=pres
  REAL                          :: VLVL             ! vertical level

!    determine which output coordinate VERTICAL={P|S}
     IF(SLVL(1).EQ.0.0)THEN  ! in sigma coordinates ground (index=1) should be 1.0
        KVER=KSYS            ! output same as input, no transformation to sigma
     ELSE
        KVER=1               ! set to sigma as the output coordinate
     END IF

! update the header with grid structure information (first 9 bytes from main) 
  I=10 
  DO K=1,12
     IF(GRIDS(K).GE.1000.0)THEN
        WRITE(HEADER(I:I+6),'(F7.2)')GRIDS(K)
     ELSEIF(GRIDS(K).GE.100.0)THEN
        WRITE(HEADER(I:I+6),'(F7.3)')GRIDS(K)
     ELSEIF(GRIDS(K).GE.10.0)THEN
        WRITE(HEADER(I:I+6),'(F7.4)')GRIDS(K)
     ELSEIF(GRIDS(K).GE.1.0)THEN
        WRITE(HEADER(I:I+6),'(F7.5)')GRIDS(K)
     ELSEIF(GRIDS(K).GE.0.0)THEN
        WRITE(HEADER(I:I+6),'(F7.6)')GRIDS(K)
     ELSEIF(GRIDS(K).GT.-1.0)THEN
        WRITE(HEADER(I:I+6),'(F7.5)')GRIDS(K)
     ELSEIF(GRIDS(K).GT.-10.0)THEN
        WRITE(HEADER(I:I+6),'(F7.4)')GRIDS(K)
     ELSEIF(GRIDS(K).GT.-100.0)THEN
        WRITE(HEADER(I:I+6),'(F7.3)')GRIDS(K)
     ELSEIF(GRIDS(K).GT.-1000.0)THEN
        WRITE(HEADER(I:I+6),'(F7.2)')GRIDS(K)
     ELSE
        WRITE(HEADER(I:I+6),'(E7.1)')GRIDS(K)
     END IF
     I=I+7
  END DO
  WRITE(HEADER(I:),'(3I3,I2,I4)') NX,NY,NZ,KVER,LENH

! fill the header record with the variable and level information
  I=109
  DO K=1,NZ

!    determine which output coordinate is valid
     IF(KVER.NE.1)THEN
        VLVL=PLVL(K)
     ELSE
!       when kver=1 then transformed to sigma
        VLVL=SLVL(K)
     END IF

!    precision depends upon the height coordinate
     IF(VLVL.GE.10000.0)THEN
        WRITE(HEADER(I:I+7),'(F6.0,I2)') VLVL,NUMVAR(K)
     ELSEIF(VLVL.GE.1000.0)THEN
        WRITE(HEADER(I:I+7),'(F6.1,I2)') VLVL,NUMVAR(K)
     ELSEIF(VLVL.GE.100.0.AND.PLVL(K).LT.1000.0)THEN
        WRITE(HEADER(I:I+7),'(F6.2,I2)') VLVL,NUMVAR(K)
     ELSEIF(VLVL.GE.10.0.AND.PLVL(K).LT.100.0)THEN
        WRITE(HEADER(I:I+7),'(F6.3,I2)') VLVL,NUMVAR(K)
     ELSEIF(VLVL.GE.1.0.AND.PLVL(K).LT.10.0)THEN
        WRITE(HEADER(I:I+7),'(F6.4,I2)') VLVL,NUMVAR(K)
     ELSE
        WRITE(HEADER(I:I+7),'(F6.5,I2)') VLVL,NUMVAR(K)
     END IF
     I=I+8       

     DO J=1,NUMVAR(K)
        WRITE(HEADER(I:I+7),'(A4,I3)') VARLVL(K,J),CHKSUM(K,J)
        I=I+8      
     END DO
  END DO

! write the header information to one or more index records
  IF(NDX.EQ.1)THEN
     WRITE(20,REC=NREC)LABEL,HEADER(1:LENH)
  ELSE
!    write extended character string over multiple index records
     I=1
     DO K=1,NDX
        J=I+NXY-1
        IF(K.EQ.NDX) J=I+(LENH-(NDX-1)*NXY)-1
        WRITE(20,REC=(NREC+K-1))  LABEL,HEADER(I:J)
        I=J+1
     END DO
  END IF

END SUBROUTINE windex

!-------------------------------------------------------------------------------
! packs a two dimensional real*4 array into a one dimensional integer*1 array
! using the differences between succesive elements

SUBROUTINE REPACK(RVAR,KPACK,PREC,NEXP,VAR1,KSUM)

  IMPLICIT none  

  REAL,       INTENT(IN)  :: RVAR(:,:)   ! data array to be packed  
  INTEGER(1), INTENT(OUT) :: KPACK(:)    ! packed int*1 output array
  REAL,       INTENT(OUT) :: PREC        ! precision of packed data array
  INTEGER,    INTENT(OUT) :: NEXP        ! packing scaling exponent
  REAL,       INTENT(OUT) :: VAR1        ! value of real array at position (1,1)
  INTEGER,    INTENT(OUT) :: KSUM        ! rotating checksum

  INTEGER     :: I,J,K,NX,NY,IVAL
  REAL        :: SEXP,RMAX,ROLD,RCOL

  NX=SIZE(rvar,1)
  NY=SIZE(rvar,2)
  VAR1=RVAR(1,1)
  ROLD=VAR1
  RMAX=0.0

! find the maximum difference between adjacent elements
  DO J=1,NY
     DO I=1,NX
!       compute max difference between elements along row
        RMAX=MAX( ABS(RVAR(I,J)-ROLD),RMAX )
        ROLD=RVAR(I,J)
     END DO
!    row element 1 difference always from previous row
     ROLD=RVAR(1,J)
  END DO

! compute the required scaling exponent
  SEXP=0.0
  IF(RMAX.NE.0.0) SEXP=LOG(RMAX)/LOG(2.)
  NEXP=INT(SEXP)  

! positive or whole number scaling round up for lower precision
  IF(SEXP.GE.0.0.OR.MOD(SEXP,1.0).EQ.0.0)NEXP=NEXP+1

! precision range is -127 to 127 or 254
  PREC=(2.0**NEXP)/254.0
  SEXP=2.0**(7-NEXP)

  K=0
  KSUM=0
  RCOL=VAR1   ! set column1 value
! pack the array from low to high
  DO J=1,NY
     ROLD=RCOL
     DO I=1,NX
        K=K+1
        IVAL=INT((RVAR(I,J)-ROLD)*SEXP+127.5) ! packed integer at element
        ROLD=FLOAT(IVAL-127)/SEXP+ROLD        ! previous element unpacked
        IF(I.EQ.1)RCOL=ROLD                   ! save the 1st col element for next row
        KPACK(K)=INT(IVAL,kind=1)             ! convert 4-byte int to 1-byte int

        KSUM=KSUM+IVAL                        ! accumulate rotating checksum
        IF(KSUM.GE.256)KSUM=KSUM-255          ! when carries over the 8th bit then reset
     END DO
  END DO

END SUBROUTINE repack


!-------------------------------------------------------------------------------
! interpolate the input data from absolute pressure coordinates to pressure-    
! sigma coordinates such that sigma(k) = (Pk-Ptop)/(Psfc-Ptop)

SUBROUTINE VSIGMA(FCHAR0,FCHAR1,PTOP,PLVL,SLVL,SM,VM)

  IMPLICIT none

  CHARACTER(4), INTENT(IN)    :: FCHAR0(:)   ! surface variables
  CHARACTER(4), INTENT(IN)    :: FCHAR1(:)   ! upper variables 
  REAL,         INTENT(IN)    :: PTOP        ! nominal model top for sigma 
  REAL,         INTENT(IN)    :: PLVL(:)     ! input level pressure
  REAL,         INTENT(IN)    :: SLVL(:)     ! output sigma level
  REAL,         INTENT(IN)    :: SM(:,:,:)   ! surface meteorology       
  REAL,         INTENT(INOUT) :: VM(:,:,:,:) ! vertical meteorology       

  INTEGER             :: KSPP,KSTT,KSUU,KSVV ! surface variable indicies
  INTEGER             :: KLHH,KLTT,KLUU,KLVV ! upper level variable indicies

  LOGICAL             :: DIAG        ! turn on local diagnostics
  INTEGER             :: I,J,K,L     ! temporary grid loop indicies
  INTEGER             :: KK          ! index on input pressure level above output
  INTEGER             :: NX,NY,NZ    ! number of grid points
  INTEGER             :: NV          ! number of variables in the vertical
  INTEGER             :: NS          ! number of variables at the surface
  REAL                :: FRAC        ! vertical interpolation fraction
  REAL                :: DHDS        ! change in height per sigma
  REAL, ALLOCATABLE   :: VSIG(:,:)   ! vertical data profile for all variables
  REAL, ALLOCATABLE   :: SINP(:)     ! sigma for each input data level 
  REAL, ALLOCATABLE   :: SOUT(:)     ! sigma for each output data level 

  NX=SIZE(VM,1)
  NY=SIZE(VM,2)
  NZ=SIZE(VM,3)
  NV=SIZE(VM,4)
  NS=SIZE(SM,3)

! find the indicies for critical surface variables
  KSPP=0
  KSTT=0
  KSUU=0
  KSVV=0
  DO J=1,NS
     IF(FCHAR0(J).EQ.'PRSS') KSPP=J
     IF(FCHAR0(J).EQ.'T02M') KSTT=J
     IF(FCHAR0(J).EQ.'U10M') KSUU=J
     IF(FCHAR0(J).EQ.'V10M') KSVV=J
  END DO

! find the indicies for critical upper-level variables
  KLHH=0
  KLTT=0
  KLUU=0
  KLVV=0
  DO J=1,NV
     IF(FCHAR1(J).EQ.'HGTS') KLHH=J
     IF(FCHAR1(J).EQ.'TEMP') KLTT=J
     IF(FCHAR1(J).EQ.'UWND') KLUU=J
     IF(FCHAR1(J).EQ.'VWND') KLVV=J
  END DO

  IF(KSPP.EQ.0)THEN
     WRITE(*,*)'Unable to find the surface pressure variable PRSS!'
     WRITE(*,*)'Vertical pressure to sigma coordinate remapping not possible!'
     STOP
  END IF

! allocate the interpolation output array for a grid point
  ALLOCATE (VSIG(NZ,NV),SINP(NZ),SOUT(NZ))

  DIAG=.FALSE.

  DO I=1,NX
  DO J=1,NY

!    note that plvl(1) and slvl(1) are at the surface and dimensioned to NZ+1
!    for data in SM and VM, k=1 is equivalent to plvl(2) and slvl(2) 
     DO K=2,NZ+1
        SINP(K-1)=(PLVL(K)-PTOP)/(SM(I,J,KSPP)-PTOP)
        SOUT(K-1)= SLVL(K)
     END DO

!    turn on diagnostics depending upon which condition is of most interest
!### IF(SM(I,J,KSPP).LT.PLVL(2)) DIAG=.TRUE.  ! ground surface above first level
!### IF(SM(I,J,KSPP).GT.PLVL(2)) DIAG=.TRUE.  ! ground surface below first level
     IF(DIAG)THEN
        WRITE(*,'(A)')' '
        WRITE(*,'(A,2I5)')'Output at I,J: ',I,J
        WRITE(*,'(A,99F6.3)')'Input sigma:  ',SINP(1:NZ)
        WRITE(*,'(A,99F6.3)')'Output sigma: ',SOUT(1:NZ)
        WRITE(*,'(A)')' '
        WRITE(*,'(A15,A22)')'Output Profile','Input profile'         
     END IF

!    loop through each output sigma level with index K
!    and find the input levels needed for interpolation
     KK=1          ! input array
     DO K=1,NZ     ! output array

!       only move up input data profile with KK when output level K is
!       not contained between KK and KK-1
        IF(SOUT(K).LE.SINP(KK))THEN
           DO WHILE (SOUT(K).LE.SINP(KK).AND.KK.LE.NZ)
              KK=KK+1 
           END DO
        END IF

        IF(KK.GT.1.AND.KK.LE.NZ)THEN
!          interpolation fraction where output between two input levels
           FRAC=(SINP(KK-1)-SOUT(K))/(SINP(KK-1)-SINP(KK))
           DO L=1,NV   
              VSIG(K,L)=FRAC*(VM(I,J,KK,L)-VM(I,J,KK-1,L))+VM(I,J,KK-1,L)
           END DO

        ELSEIF(KK.GT.NZ)THEN
!          when the output sigma level is above the highest input level, 
!          extrapolate the output data array from the highest input level
!          except for height and temperature, other values persist
           VSIG(K,:)=VM(I,J,NZ,:)

!          insure height and temperature field change, all others are unchanged
           IF(KLHH.NE.0)THEN
!             compute input data change in height with sigma
              DHDS=(VM(I,J,NZ,KLHH)-VM(I,J,NZ-1,KLHH))/(SINP(NZ-1)-SINP(NZ))
              VSIG(K,KLHH)=(SOUT(K-1)-SOUT(K))*DHDS+VSIG(K-1,KLHH)
           END IF
           IF(KLTT.NE.0)THEN
!             compute input data change in temperature with sigma
              DHDS=(VM(I,J,NZ,KLTT)-VM(I,J,NZ-1,KLTT))/(SINP(NZ-1)-SINP(NZ))
              VSIG(K,KLTT)=(SOUT(K-1)-SOUT(K))*DHDS+VSIG(K-1,KLTT)
           END IF

        ELSE
!          one or more output sigma levels may be below the lowest input level
!          therefore the wind and temperature are interpolated between the 
!          surface and the first input level above the surface 
           VSIG(K,:)=VM(I,J,1,:)   ! load all variables at the lowest level
           FRAC=(1.0-SOUT(K))/(1.0-SINP(1))

!          low-level temperature
           IF(KSTT.NE.0.AND.KLTT.NE.0)   & 
              VSIG(K,KLTT)=FRAC*(VM(I,J,1,KLTT)-SM(I,J,KSTT))+SM(I,J,KSTT)

!          low-level U-wind      
           IF(KSUU.NE.0.AND.KLUU.NE.0)   & 
              VSIG(K,KLUU)=FRAC*(VM(I,J,1,KLUU)-SM(I,J,KSUU))+SM(I,J,KSUU)

!          low-level V-wind      
           IF(KSVV.NE.0.AND.KLVV.NE.0)   & 
              VSIG(K,KLVV)=FRAC*(VM(I,J,1,KLVV)-SM(I,J,KSVV))+SM(I,J,KSVV)
        END IF

!       diagnostic output only once per entry
        IF(DIAG)THEN
           IF(KK.LE.NZ)THEN
               WRITE(*,'(I3,F6.3,F6.1,I4,F6.3,2F6.1)')             &
               K,SOUT(K),VSIG(K,KLTT),KK,SINP(KK),PLVL(KK),VM(I,J,KK,KLTT)
           ELSE
               WRITE(*,'(I3,F6.3,F6.1)') K,SOUT(K),VSIG(K,KLTT)
           END IF
        END IF 
     END DO
     IF(DIAG)THEN
        DIAG=.FALSE.
!###    READ(*,*)
     END IF

!    at each grid point replace input data with interpolated profile
     VM(I,J,:,:)=VSIG     

  END DO
  END DO
  DEALLOCATE (VSIG,SINP,SOUT)

END SUBROUTINE vsigma
