PROGRAM avckf
        
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: SP = selected_real_kind(5)
    INTEGER, PARAMETER :: DP = selected_real_kind(10)
    
    INTEGER  :: iu=100, ierr, n, i, j
    REAL(DP), ALLOCATABLE :: eig(:)
    COMPLEX(DP), ALLOCATABLE :: avck(:,:)
    
    TYPE bse_index
       INTEGER :: nvck
       INTEGER, ALLOCATABLE :: ik(:)
       INTEGER, ALLOCATABLE :: iv(:)
       INTEGER, ALLOCATABLE :: ic(:)
       INTEGER, ALLOCATABLE :: is(:)
    END TYPE
  
    TYPE (bse_index) :: bd
    
    OPEN(UNIT=iu, FILE='AVCKCAR', FORM='UNFORMATTED', STATUS='OLD', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(*, *) 'error while openning AVCKCAR'
       STOP
    END IF
    
    READ(iu) bd%nvck
    
    n = bd%nvck
    ALLOCATE( bd%ik(n), bd%iv(n), bd%ic(n), bd%is(n) )
    ALLOCATE( eig(n) )
    ALLOCATE( avck(n,n) )
    
    READ(iu) bd%ik, bd%iv, bd%ic, bd%is
    READ(iu) eig
    READ(iu) avck
      
    CLOSE(iu)
    
    OPEN(UNIT=iu, FILE='AVCKCARF', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(*, *) 'error while openning AVCKCARF'
       STOP
    END IF
    
    DO i = 1, n
       WRITE(iu,'(I6,A15,F14.8)') i, 'BSE eigenvalue:', eig(i) 
       DO j = 1, n
          WRITE(iu, '(3I6,2F14.6)')bd%ik(j),bd%iv(j),bd%ic(j),avck(j,i)
       ENDDO
    ENDDO
    
    CLOSE(iu)
    
END PROGRAM
