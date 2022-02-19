	PROGRAM EXTERN
	CHARACTER COMMAND*1000, FINAL*1000
	INTEGER INPUT
	CHARACTER(LEN=100) :: STR
	INPUT = 27;
	COMMAND = 'External.exe -override=a=';
	write(STR , *) INPUT
	FINAL = TRIM(COMMAND) // ADJUSTL(STR)
	CALL SYSTEM(FINAL)
	CALL READ_CSV
	END PROGRAM EXTERN
	
		SUBROUTINE READ_CSV
        INTEGER i,j, ncols, nrows, final
        CHARACTER*4 Title(20)
        DOUBLE PRECISION x(100,20)
        OPEN(1,FILE= 'External_res.csv', STATUS='old')
        ncols=2
        nrows=3
        READ (1,*) (Title(j),j=1,ncols)
!       WRITE (*,"(100(1x,a))") "titles:",title(:ncols)
        DO i=1, nrows
           READ (1,*) (x(i,j),j=1,ncols)
!          WRITE (*,"(1x,i4,100(1x,f0.2))") i,x(i,:ncols)
        END DO
		final=x(1,2)
		WRITE (INTEGER , *) final
		END SUBROUTINE READ_CSV
	