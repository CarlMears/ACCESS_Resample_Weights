C     7/20/2000 version was changed on 11/17/2001
C     1. Code was restructured: Inteface to subroutine command was changed to inteface and subroutine was put within interface
C     2. Intefaces blocks were put into routines
C     3. Dimension(:) replace by Dimension(*)
C     4. $nofreeform was removed
C     These changed need to be made to work with new COMPAQ compiler



C	LOADS AND DECOMPRESSES GZIPPED DATA FROM FILE 
C     GIVEN BY FILENAME INTO DATARRY, 
C	AND RETURNS NUMBER OF BYTES DECOMPRESSED IN DATACNT
C	THIS SHOULD BE THE SAME AS RETURNED BY UNCMPSIZ

	SUBROUTINE FASTUNZIP(FILENAME, ISIZE,ABUF,IERR)

	INTERFACE
	SUBROUTINE GZIPINFL [C,ALIAS:'_GZIP_inflate']
	1 (FILENAME, DATARRY,DATACNT,ERR)
	CHARACTER*(*) FILENAME [REFERENCE]
	CHARACTER(1),DIMENSION(*)::DATARRY [REFERENCE]
	INTEGER(4) DATACNT [REFERENCE]
	INTEGER(4) ERR [REFERENCE]
	END	SUBROUTINE
	END INTERFACE

	CHARACTER(1) CHAR
	CHARACTER*(*) FILENAME
	CHARACTER(101) FILENAMEC
	CHARACTER(1) ABUF(*)

      NN=LEN_TRIM(FILENAME)
	FILENAMEC=FILENAME(1:NN) // CHAR(0)
	CALL GZIPINFL(FILENAMEC, ABUF,ISIZE,IERR)
	RETURN
	END




C	COMPRESSES DATA IN DATARRY IN GZIP FORMAT AND
C     SAVES IN FILE GIVEN BY FILENAME, 
C	AND RETURNS NUMBER OF BYTES DECOMPRESSED IN DATACNT
C	THIS SHOULD BE THE SAME AS RETURNED BY UNCMPSIZ

	SUBROUTINE FASTGZIP(FILENAME,ISIZE,ABUF, IERR)

	INTERFACE
	SUBROUTINE GZIPDEFL [C,ALIAS:'_GZIP_deflate'](FILENAME, DATARRY,DATACNT,ERR)
	CHARACTER*(*) FILENAME [REFERENCE]
	CHARACTER(1),DIMENSION(*)::DATARRY [REFERENCE]
	INTEGER(4) DATACNT [REFERENCE]
	INTEGER(4) ERR [REFERENCE]
	END	SUBROUTINE
	END INTERFACE

	CHARACTER(1) CHAR
	CHARACTER*(*) FILENAME
	CHARACTER(101) FILENAMEC
	CHARACTER(1) ABUF(*)

      NN=LEN_TRIM(FILENAME)
	FILENAMEC=FILENAME(1:NN) // CHAR(0)
	CALL GZIPDEFL(FILENAMEC,ABUF,ISIZE, IERR)

	RETURN
	END




C     DETERMINES UNCOMPRESS SIZE

	SUBROUTINE FIND_UNZIP_SIZE(FILENAME, ISIZE)

	INTERFACE
	SUBROUTINE UNCMPSIZ [C,ALIAS:'_uncompressed_size']
	1 (FILENAME, DATASIZE)
	CHARACTER*(*) FILENAME [REFERENCE]
	INTEGER(4) DATASIZE [REFERENCE]
	END	SUBROUTINE
	END INTERFACE

	CHARACTER(1) CHAR
	CHARACTER*(*) FILENAME
	CHARACTER(101) FILENAMEC
      NN=LEN_TRIM(FILENAME)
	FILENAMEC=FILENAME(1:NN) // CHAR(0)
	CALL UNCMPSIZ(FILENAMEC, ISIZE)
	RETURN
	END