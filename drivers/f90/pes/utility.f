
      SUBROUTINE POTINFO

C***********************************************************************
C   VERSION HISTORY
C
C   Version 1.0
C   Author:              R. J. Duchovic 
C   Date:                99-10-15   
C                        99-11-05   
C                        99-11-09   
C                        99-11-10   
C                        99-11-17   
C                        99-11-19   
C                        99-12-09   
C                        99-12-10   
C                        99-12-14   
C                        01-07-18   
C                                  
C   Version 1.0.1
C   Corrections by:      Jose Corchado
C   Date:                May 13, 2003
C   Summary of changes:  This new version of utility.f corrects two
C                        errors in the original version.  
C                        Corrected lines in are marked 'JC05132003' in 
C                        the margin.  
C                        The corrections are:
C                        1) In the part of the subroutine CARTTOR where
C                           the code handles ICARTR = 2, the counter 
C                           "I" was added.  Previously, not all of the 
C                           components of R(I) were getting the correct
C                           values.
C                        2) In the part of the subroutine RTOCART where
C                           the code handles ICARTR = 2, the DO loop
C                           was changed to run from 1 to NATOMS.  
C                           Previously, some of the derivatives with 
C                           resepect to the Cartesian coordinates were 
C                           incorrectly returned as zero.
C
C***********************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER(75) REF(5)
C
      PARAMETER (N3ATOM=75)
      PARAMETER (NATOM=25)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
      PARAMETER (PI = 3.141592653589793D0)
C
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON /PT3CM/  EZERO(ISURF+1)
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)
      COMMON /PT5CM/  ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON/UTILCM/ DGSCARTNU(NATOM,3),DESCARTNU(NATOM,3,ISURF),
     +               DIJCARTNU(NATOM,3,JSURF),CNVRTD,CNVRTE,
     +               CNVRTDE,IREORDER,KSDIAG,KEDIAG,KSOFFD,KEOFFD
C
C LIST CONTENTS OF REF
C
      write(NFLAG(18),96)
 96   format(/)
      do i =1,5
         write(NFLAG(18),97) REF(i)
      end do
 97   format(2x,a75)
C
C LIST AND CHECK ON NUMBER OF SURFACES
C
      WRITE(NFLAG(18),96)
      KMAX = 0
      DO I = 1,ISURF+1
         DO J = 1,ISURF+1
            IF(NASURF(I,J).NE.0.AND.KMAX.LT.MAX(I,J)) KMAX = MAX(I,J)
         ENDDO
      ENDDO
      WRITE(NFLAG(18),101) MSURF,KMAX-1
101   FORMAT(2x,' MAX. AND ACTUAL NO. OF EXCITED SURFACES: ',I3,5x,I3)
      IF(KMAX-1.GT.MSURF) THEN
         WRITE(6,*) ' WRONG INPUT ON NUMBER OF EXCITED SURFACES'
         STOP
      ENDIF
      KSDIAG = 0
      KEDIAG = 0
      DO I = 2,ISURF+1
         IF(NASURF(I,I).NE.0) THEN
            KEDIAG = I-1
            IF(KSDIAG.EQ.0) KSDIAG = I-1
         ENDIF
      ENDDO
      KSOFFD = 0
      KEOFFD = 0
      K = 0
      DO I = 1,ISURF
         DO J = I+1,ISURF+1
            K = K+1
            IF(NASURF(I,J)+NASURF(J,I).NE.0) THEN
               KEOFFD = K
               IF(KSOFFD.EQ.0) KSOFFD = K
            ENDIF
         ENDDO
      ENDDO
C
C LIST AND CHECK ON ORDER OF DERIVATIVES
C
      WRITE(NFLAG(18),103) MDER,NDER
103   FORMAT(2x,' MAX. AND ACTUAL ORDER OF DERIVATIVES:    ',I3,5x,I3)
      IF(NDER.GT.MDER) THEN
         WRITE(6,*) ' WRONG INPUT ON ORDER OF DERIVATIVES'
         STOP
      ENDIF
C
C**********************************************
C                                             *
C   NFLAG(19) is used to determine whether    *
C   or not a detailed description of the      *
C   options is printed.                       *
C                                             *
C**********************************************
C
      IF(NFLAG(19).EQ.1) THEN
         write(NFLAG(18),100)
 100     format(/)
         write(NFLAG(18),120)
 120     format(2x,'Cartesian coordinates are supplied by',/,
     +          2x,'the user in the array CART.',//)
         write(NFLAG(18),125)
 125     format(2x,'Provide cartesian coordinates in the',/,
     +          2x,'following order using the array CART',//,
     +          2x,' CART(1,1)...CART(1,3)   => ATOM 1',/,
     +          2x,' CART(2,1)...CART(2,3)   => ATOM 2',/,
     +          2x,' CART(3,1)...CART(3,3)   => ATOM 3',/,
     +          2x,' CART(N,1)...CART(N,3)   => ATOM N',/,
     +          2x,'CART(25,1)...CART(25,3)  => ATOM 25',/)
C
C**********************************************
C                                             *
C        Cartesian Coordinates                *
C                                             *
C   CART(1,1)  ... CART(1,3)  => ATOM 1       *
C   CART(2,1)  ... CART(2,3)  => ATOM 2       *
C   CART(3,1)  ... CART(3,3)  => ATOM 3       *
C       .               .                     *
C       .               .                     *
C       .               .                     *
C   CART(N,1)  ... CART(N,3)  => ATOM N       *
C       .               .                     *
C       .               .                     *
C       .               .                     *
C   CART(25,1) ... CART(25,3) => ATOM 25      *
C                                             *
C**********************************************
C
         write(NFLAG(18),130)
 130     format(2x,'If the user wishes to relabel the atoms,',/,
     +          2x,'set the variable IREORDER equal to 1',/,
     +          2x,'in the PARAMETER statement.  The user',/,
     +          2x,'must also specify the new labeling',/,
     +          2x,'scheme.  This is done using the array',/,
     +          2x,'NULBL in the following manner:',//,
     +          2x,'NULBL(i) = j',/,
     +          2x,'where:  i => old index',/,
     +          2x,'        j => new index',//)
         write(NFLAG(18),150)
 150     format(2x,'Cartesian coordinates can be provided to',/,
     +          2x,'the potential routine in a variety of units.',/,
     +          2x,'The input units will be converted to Bohr',/,
     +          2x,'based on the following values of the NFLAG',/,
     +          2x,'variable:',//,
     +          2x,'NFLAG(1)  =  1  =>  CARTESIANS IN BOHR (no',/,
     +          2x,'                    conversion required)',/,
     +          2x,'NFLAG(1)  =  2  =>  CARTESIANS IN ANGSTROMS',//)
         write(NFLAG(18),160)
 160     format(2x,'The value of the energy and derivatives',/,
     +          2x,'(if computed) can be reported in a variety',/,
     +          2x,'units.  A units conversion will take place',/,
     +          2x,'as specified by the following values of the',/,
     +          2x,'NFLAG variable:',//,
     +          2x,'NFLAG(2) = 1 =>  ENERGIES REPORTED IN HARTEEE',/,
     +          2x,'NFLAG(2) = 2 =>  ENERGIES REPORTED IN mHARTREE',/,
     +          2x,'NFLAG(2) = 3 =>  ENERGIES REPORTED IN eV',/,
     +          2x,'NFLAG(2) = 4 =>  ENERGIES REPORTED IN kcal/mol',/,
     +          2x,'NFLAG(2) = 5 =>  ENERGIES REPORTED IN cm**-1',//)
         write(NFLAG(18),165)
 165     format(2x,'A units conversion will take place',/,
     +       2x,'as specified by the following values of the',/,
     +       2x,'NFLAG variable:',//,
     +       2x,'NFLAG(1)=1 & NFLAG(2)=1 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           HARTEEE/BOHR',/,
     +       2x,'NFLAG(1)=1 & NFLAG(2)=2 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           mHARTREE/BOHR',/,
     +       2x,'NFLAG(1)=1 & NFLAG(2)=3 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           eV/BOHR',/,
     +       2x,'NFLAG(1)=1 & NFLAG(2)=4 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           kcal/mol/BOHR',/,
     +       2x,'NFLAG(1)=1 & NFLAG(2)=5 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           cm**-1/BOHR',//)
         write(NFLAG(18),170)
 170     format(2x,'A units conversion will take place',/,
     +       2x,'as specified by the following values of the',/,
     +       2x,'NFLAG variable:',//,
     +       2x,'NFLAG(1)=2 & NFLAG(2)=1 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           HARTEEE/ANGSTROM',/,
     +       2x,'NFLAG(1)=2 & NFLAG(2)=2 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           mHARTREE/ANGSTROM',/,
     +       2x,'NFLAG(1)=2 & NFLAG(2)=3 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           eV/ANGSTROM',/,
     +       2x,'NFLAG(1)=2 & NFLAG(2)=4 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           kcal/mol/ANGSTROM',/,
     +       2x,'NFLAG(1)=2 & NFLAG(2)=5 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           cm**-1/ANGSTROM',//)
      ENDIF
      RETURN
      END
 
 
 
 
C
      SUBROUTINE ANCVRT
C
C***********************************
C                                  *
C   WRITTEN BY: R. J. DUCHOVIC     *
C               A. F. WAGNER       *
C                                  *
C   DATE:               99-10-15   *
C                       99-11-09   *
C                       99-11-10   *
C                       99-11-15   *
C                       99-11-16   *
C                       99-11-17   *
C                       99-11-19   *
C                       99-12-09   *
C                       00-01-23   *
C                       01-07-18   *
C                                  *
C***********************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER(75) REF(5)
      CHARACTER(3) PERIODIC_1(7,32)
C
      PARAMETER (N3ATOM=75)
      PARAMETER (NATOM=25)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      CHARACTER(2) NAME1(NATOM)
      CHARACTER(2) NAME2(NATOM)
      CHARACTER(1) IBLANK
      CHARACTER(20) DISTANCE
      CHARACTER(20) UNITS
C
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON /PT3CM/  EZERO(ISURF+1)
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)
      COMMON /PT5CM/  ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/UTILCM/ DGSCARTNU(NATOM,3),DESCARTNU(NATOM,3,ISURF),
     +               DIJCARTNU(NATOM,3,JSURF),CNVRTD,CNVRTE,
     +               CNVRTDE,IREORDER,KSDIAG,KEDIAG,KSOFFD,KEOFFD
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      DIMENSION IANUM(7,32)
      DIMENSION ISAVE(NATOM),JSAVE(NATOM)
C
C***********************************************************************
C                                                                      *
C  PI:                          3.141592653589793                      *
C                                                                      *
C  FUNDAMENTAL  CONSTANTS (1998 CODATA values)                         *
C                                                                      *
C  SPPED OF LIGHT:              2.99792458 x 10**8   ms**-1            *
C  PERMEABILITY  OF VACUUM:     4.0D0*PI*1.0 x 10**-07  NA**-2         *
C  PERMITIVITY  OF VACUUM:      1.0/(CMU_0*CLIGHT**2) Fm**-1           *
C  ELEMENTARY  CHARGE:          1.602176462 x 10**-19  C               *
C  PLANCK'S CONSTANT            6.62606876 x 10**-34  Js               *
C  ELECTRON MASS:               9.10938188 x 10**-31  kg               * 
C  ANGSTROM:                    1.0 X 10**-10 m                        *
C  AVOGADRO'S CONSTANT:         6.02214199 x 10**23 mol**-1            *
C  KILOCALORIE:                 4.184  x 10**10 ergs                   *
C                                                                      *
C  BOHR     = 2.0D0*CPLANCK**2/(2.0D0*PI*CM_E*CLIGHT**2*CMU_0*CE**2)   *
C           = 5.291772083 x 10**-11 m                                  *
C  HTOMILLH = 1000 milliHartree                                        *
C  HTOEV    = CE/(4.0D0*PI*CEPSILON_0*BOHR)                            *
C           = 27.2113834 eV                                            *
C  HTOKCAL  = (CE**2*CAVOGADRO)/(4.0D0*PI*CEPSILON_0*BOHR*CKCAL)       *
C           = 627.509470 kcal/mole                                     *
C  HTOWAVE  = (CE**2)/(4.0D0*PI*CEPSILON_0*BOHR*CPLANCK*CLIGHT)        *
C           = 219474.631 Hartree                                       *
C  HTOKJ    = HTOKCAL*4.184                                            *
C           = 2625.49962                                               *
C                                                                      *
C***********************************************************************
C
      PARAMETER(        PI = 3.141592653589793D0)
      PARAMETER(    CLIGHT = 2.99792458D08)
      PARAMETER(     CMU_0 = 4.0D0*PI*1.0D-07)
      PARAMETER(CEPSILON_0 = 1.0D0/(CMU_0*CLIGHT**2))
      PARAMETER(        CE = 1.602176462D-19)
      PARAMETER(   CPLANCK = 6.62606876D-34)
      PARAMETER(      CM_E = 9.10938188D-31)
      PARAMETER(      CANG = 1.0D-10)
      PARAMETER( CAVOGADRO = 6.02214199D23)
      PARAMETER(     CKCAL = 4.184D10)
      PARAMETER(  HTOMILLH = 1000.D0)
      PARAMETER(     HTOEV = 27.2113834D0)
      PARAMETER(   HTOKCAL = 627.509470D0)
      PARAMETER(   HTOWAVE = 219474.631D0)
      PARAMETER(     HTOKJ = 2625.49962D0)
      PARAMETER(    BOHR_A = .5291772083D0)
C
      DO I=1,7
         DO J=1,32
            IANUM(I,J)=0
            PERIODIC_1(I,J)=' '
         END DO
      END DO
C
      DISTANCE = 'BOHR                '
      UNITS    = 'HARTREE             '
C
      IANUM(1,1)  =  1
      IANUM(1,32) =  2
      IANUM(2,1)  =  3
      IANUM(2,2)  =  4
      IANUM(2,27) =  5
      IANUM(2,28) =  6
      IANUM(2,29) =  7
      IANUM(2,30) =  8
      IANUM(2,31) =  9
      IANUM(2,32) = 10
      IANUM(3,1)  = 11
      IANUM(3,2)  = 12
      IANUM(3,27) = 13
      IANUM(3,28) = 14
      IANUM(3,29) = 15
      IANUM(3,30) = 16
      IANUM(3,31) = 17
      IANUM(3,32) = 18
      IANUM(4,1)  = 19
      IANUM(4,2)  = 20
      IANUM(4,17) = 21
      IANUM(4,18) = 22
      IANUM(4,19) = 23
      IANUM(4,20) = 24
      IANUM(4,21) = 25
      IANUM(4,22) = 26
      IANUM(4,23) = 27
      IANUM(4,24) = 28
      IANUM(4,25) = 29
      IANUM(4,26) = 30
      IANUM(4,27) = 31
      IANUM(4,28) = 32
      IANUM(4,29) = 33
      IANUM(4,30) = 34
      IANUM(4,31) = 35
      IANUM(4,32) = 36
      IANUM(5,1)  = 37
      IANUM(5,2)  = 38
      IANUM(5,17) = 39
      IANUM(5,18) = 40
      IANUM(5,19) = 41
      IANUM(5,20) = 42
      IANUM(5,21) = 43
      IANUM(5,22) = 44
      IANUM(5,23) = 45
      IANUM(5,24) = 46
      IANUM(5,25) = 47
      IANUM(5,26) = 48
      IANUM(5,27) = 49
      IANUM(5,28) = 50
      IANUM(5,29) = 51
      IANUM(5,30) = 52
      IANUM(5,31) = 53
      IANUM(5,32) = 54
      IANUM(6,1)  = 55
      IANUM(6,2)  = 56
      IANUM(6,3)  = 57
      IANUM(6,4)  = 58
      IANUM(6,5)  = 59
      IANUM(6,6)  = 60
      IANUM(6,7)  = 61
      IANUM(6,8)  = 62
      IANUM(6,9)  = 63
      IANUM(6,10) = 64
      IANUM(6,11) = 65
      IANUM(6,12) = 66
      IANUM(6,13) = 67
      IANUM(6,14) = 68
      IANUM(6,15) = 69
      IANUM(6,16) = 70
      IANUM(6,17) = 71
      IANUM(6,18) = 72
      IANUM(6,19) = 73
      IANUM(6,20) = 74
      IANUM(6,21) = 75
      IANUM(6,22) = 76
      IANUM(6,23) = 77
      IANUM(6,24) = 78
      IANUM(6,25) = 79
      IANUM(6,26) = 80
      IANUM(6,27) = 81
      IANUM(6,28) = 82
      IANUM(6,29) = 83
      IANUM(6,30) = 84
      IANUM(6,31) = 85
      IANUM(6,32) = 86
      IANUM(7,1)  = 87
      IANUM(7,2)  = 88
      IANUM(7,3)  = 89
      IANUM(7,4)  = 90
      IANUM(7,5)  = 91
      IANUM(7,6)  = 92
      IANUM(7,7)  = 93
      IANUM(7,8)  = 94
      IANUM(7,9)  = 95
      IANUM(7,10) = 96
      IANUM(7,11) = 97
      IANUM(7,12) = 98
      IANUM(7,13) = 99
      IANUM(7,14) = 100
      IANUM(7,15) = 101
      IANUM(7,16) = 102
      IANUM(7,17) = 103
      IANUM(7,18) = 104
      IANUM(7,19) = 105
      IANUM(7,20) = 106
      IANUM(7,21) = 107
      IANUM(7,22) = 108
      IANUM(7,23) = 109
      IANUM(7,24) = 110
      IANUM(7,25) = 111
      IANUM(7,26) = 112
      IANUM(7,27) = 113
      IANUM(7,28) = 114
      IANUM(7,29) = 115
      IANUM(7,30) = 116
      IANUM(7,31) = 117
      IANUM(7,32) = 120
C
      PERIODIC_1(1,1)   = 'H  '
      PERIODIC_1(1,32)  = 'He '
      PERIODIC_1(2,1)   = 'Li '
      PERIODIC_1(2,2)   = 'Be '
      PERIODIC_1(2,27)  = 'B  '
      PERIODIC_1(2,28)  = 'C  '
      PERIODIC_1(2,29)  = 'N  '
      PERIODIC_1(2,30)  = 'O  '
      PERIODIC_1(2,31)  = 'F  '
      PERIODIC_1(2,32)  = 'Ne '
      PERIODIC_1(3,1)   = 'Na '
      PERIODIC_1(3,2)   = 'Mg '
      PERIODIC_1(3,27)  = 'Al '
      PERIODIC_1(3,28)  = 'Si '
      PERIODIC_1(3,29)  = 'P  '
      PERIODIC_1(3,30)  = 'S  '
      PERIODIC_1(3,31)  = 'Cl '
      PERIODIC_1(3,32)  = 'Ar '
      PERIODIC_1(4,1)   = 'K  '
      PERIODIC_1(4,2)   = 'Ca '
      PERIODIC_1(4,17)  = 'Sc '
      PERIODIC_1(4,18)  = 'Ti '
      PERIODIC_1(4,19)  = 'V  '
      PERIODIC_1(4,20)  = 'Cr '
      PERIODIC_1(4,21)  = 'Mn '
      PERIODIC_1(4,22)  = 'Fe '
      PERIODIC_1(4,23)  = 'Co '
      PERIODIC_1(4,24)  = 'Ni '
      PERIODIC_1(4,25)  = 'Cu '
      PERIODIC_1(4,26)  = 'Zn '
      PERIODIC_1(4,27)  = 'Ga '
      PERIODIC_1(4,28)  = 'Ge '
      PERIODIC_1(4,29)  = 'As '
      PERIODIC_1(4,30)  = 'Se '
      PERIODIC_1(4,31)  = 'Br '
      PERIODIC_1(4,32)  = 'Kr '
      PERIODIC_1(5,1)   = 'Rb '
      PERIODIC_1(5,2)   = 'Sr '
      PERIODIC_1(5,17)  = 'Y  '
      PERIODIC_1(5,18)  = 'Zr '
      PERIODIC_1(5,19)  = 'Nb '
      PERIODIC_1(5,20)  = 'Mo '
      PERIODIC_1(5,21)  = 'Tc '
      PERIODIC_1(5,22)  = 'Ru '
      PERIODIC_1(5,23)  = 'Rh '
      PERIODIC_1(5,24)  = 'Pd '
      PERIODIC_1(5,25)  = 'Ag '
      PERIODIC_1(5,26)  = 'Cd '
      PERIODIC_1(5,27)  = 'In '
      PERIODIC_1(5,28)  = 'Sn '
      PERIODIC_1(5,29)  = 'Sb '
      PERIODIC_1(5,30)  = 'Te '
      PERIODIC_1(5,31)  = 'I  '
      PERIODIC_1(5,32)  = 'Xe '
      PERIODIC_1(5,32)  = 'Xe '
C
C
      DO I=1,NATOMS
         ISAVE(I)=0
         JSAVE(I)=0
         NAME1(I)='  '
         NAME2(I)='  '
      END DO
      IBLANK=' '
C
      DO IND=1,NATOMS
         DO I=1,7
            DO J=1,32
               IF(INDEXES(IND).EQ.IANUM(I,J)) THEN
                  ISAVE(IND)=I
                  JSAVE(IND)=J
               END IF
            END DO
         END DO
      END DO
      write(NFLAG(18),100)
 100  format(/,2x,'The potential routine in POTLIB assumes a',/,
     +         2x,'default labeling of the atoms.  The user has',/,
     +         2x,'made the following selection of atomic labels:',/)
      write(NFLAG(18),125)
 125  format(15x,'DEFAULT',25x,'USER-SELECTED',/)
 
      DO IND=1,NATOMS
         IND2=NULBL(IND)
         IF(IND2.EQ.0) IND2=IND
         WRITE(NFLAG(18),1000) IND,
     +                 PERIODIC_1(ISAVE(IND),JSAVE(IND)),
     +                 PERIODIC_1(ISAVE(IND2),JSAVE(IND2))
 1000    FORMAT(2x,'Atom ',I2,' is ',6x,a2,34x,a2)
      END DO
      WRITE(NFLAG(18),1500)
 1500 FORMAT(/)
C
C*******************************
C                              *
C   Determine Reactant #1      *
C                              *
C*******************************
C
      INC1=0
      DO IND=1,IRCTNT-1
         INC1=INC1+1
         NAME1(INC1)=PERIODIC_1(ISAVE(IND),JSAVE(IND))(:2)
         WRITE(NFLAG(18),2000) IND,
     +                    PERIODIC_1(ISAVE(IND),JSAVE(IND))
 2000    FORMAT(2x,'Atom ',I2,2x,a2,' is a member of Reactant #1')
      END DO
C
C*******************************
C                              *
C   Determine Reactant #2      *
C                              *
C*******************************
C
      INC2=0
      DO IND=IRCTNT,NATOMS
         INC2=INC2+1
         NAME2(INC2)=PERIODIC_1(ISAVE(IND),JSAVE(IND))(:2)
         WRITE(NFLAG(18),2005) IND,
     +                    PERIODIC_1(ISAVE(IND),JSAVE(IND))
 2005    FORMAT(2x,'Atom ',I2,2x,a2,' is a member of Reactant #2')
      END DO
      write(NFLAG(18),140)
 140  format(2x,/,2x,'The Potential Energy is Zero for the following',/,
     +       2x,'geometric arrangement:  The reactants are at',/,
     +       2x,'their equilibrium geometries, widely separated.',/)
      write(NFLAG(18),145)
 145  format(2x,/,2x,'Reactant #1')
      write(NFLAG(18),*) IBLANK,(NAME1(I),I=1,INC1)
      write(NFLAG(18),147)
 147  format(2x,/,2x,'Reactant #2')
      write(NFLAG(18),*) IBLANK,(NAME2(I),I=1,INC2)
      write(NFLAG(18),148)
 148  format(2x,/,2x,'The default units are:',/,
     +       5x,'BOHR FOR DISTANCE',/,
     +       5x,'HARTREE FOR ENERGY',/)
C
      IF(NFLAG(1).EQ.2) DISTANCE = 'ANGSTROMS           '
      IF(NFLAG(2).EQ.2) THEN
         UNITS = 'MILLIHARTREE        '
      ELSEIF(NFLAG(2).EQ.3) THEN
         UNITS = 'EV                  '
      ELSEIF(NFLAG(2).EQ.4) THEN
         UNITS = 'KCAL PER MOLE       '
      ELSEIF(NFLAG(2).EQ.5) THEN
         UNITS = 'WAVENUMBERS         '
      ELSEIF(NFLAG(2).EQ.6) THEN
         UNITS = 'KILOJOULES PER MOLE '
      ENDIF
      write(NFLAG(18),149) DISTANCE,UNITS
 149  format(2x,/,2x,'The user has chosen units :',/,
     +       5x,a20,1x,'FOR DISTANCE',/,
     +       5x,a20,1x,'FOR ENERGY',/)
C
      write(NFLAG(18),150)
 
 150  format(2x,/,2x,'To reset the zero-of-energy, assign a',/,
     +            2x,'non-zero value to the variable ANUZERO.',/,
     +            2x,'The value of ANUZERO will be SUBTRACTED from',/,
     +            2x,'the energy calculated by this routine.',/)
C
C*********************************************
C                                            *
C   SET UP CONSTANTS FOR UNITS CONVERSIONS   *
C                                            *
C*********************************************
C
C*****************************************************************
C                                                                *
C   NFLAG(1) = 1    => UNITS IN BOHR, NO CONVERSION REQUIRED     *
C   NFLAG(1) = 2    => UNITS IN ANGSTROMS, CONVERSION REQUIRED   *
C                                                                *
C   NFLAG(2) = 1 => UNITS IN HARTREE, NO CONVERSION REQUIRED     *
C            = 2 => UNITS IN MILLIHARTREE, CONVERSION REQUIRED   *
C            = 3 => UNITS IN EV, CONVERSION REQUIRED             *
C            = 4 => UNITS IN KCLAL, CONVERSION REQUIRED          *
C            = 5 => UNITS IN WAVENUMBER, CONVERSION REQUIRED     *
C            = 6 => UNITS IN KILOJOULES/MOLE, CONVERSION REQUIRED*
C                                                                *
C*****************************************************************
C
      CNVRTD = 1.D0
      CNVRTE = 1.D0
      CNVRTDE = 1.D0
C
C**************************
C                         *
C   DISTANCE CONVERSION   *
C                         *
C**************************
C
      IF(NFLAG(1).EQ.2) CNVRTD = BOHR_A
C
C**************************
C                         *
C   ENERGY CONVERSION     *
C                         *
C**************************
C
      IF(NFLAG(2).EQ.2) THEN
         CNVRTE = CNVRTE*HTOMILLH
      ELSEIF(NFLAG(2).EQ.3) THEN
         CNVRTE = CNVRTE*HTOEV
      ELSEIF(NFLAG(2).EQ.4) THEN
         CNVRTE = CNVRTE*HTOKCAL
      ELSEIF(NFLAG(2).EQ.5) THEN
         CNVRTE = CNVRTE*HTOWAVE
      ELSEIF(NFLAG(2).EQ.6) THEN
         CNVRTE = CNVRTE*HTOKJ
      ENDIF
C
C******************************
C                             *
C   DERIVATIVE CONVERSION     *
C                             *
C******************************
C
      CNVRTDE = CNVRTE/CNVRTD
C
C*********************************
C                                *
C   ASSIGN VALUE TO "IREORDER"   *
C                                *
C*********************************
C
      ISUM = 0
      DO INU=1,25
         ISUM=ISUM + NULBL(INU)
      END DO
      IREORDER = 0
      IF(ISUM.NE.0) IREORDER = 1
C
      RETURN
      END
C
      SUBROUTINE CARTOU
C
C********************************************************
C                                                       *
C   WRITTEN BY: R. J. DUCHOVIC,A. F. WAGNER             *
C                                                       *
C   DATE:               00-01-22                        *
C                       01-07-18                        *
C                                                       *
C   THIS SUBROUTINE IS DESIGNED TO CONVERT              *
C   POSITIONS CART IN USER ORDER AND UNITS INTO         *
C   POSITIONS CARTNU IN DEVELOPER ORDER AND UNITS       *
C   (DEVELOPER UNITS ARE FIXED AT BOHR).                *
C   CONVERSION FACTOR DONE IN ANCVRT.F                  *
C                                                       *
C********************************************************
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      CHARACTER(75) REF(5)
C
      PARAMETER (N3ATOM=75)
      PARAMETER (NATOM=25)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/UTILCM/ DGSCARTNU(NATOM,3),DESCARTNU(NATOM,3,ISURF),
     +               DIJCARTNU(NATOM,3,JSURF),CNVRTD,CNVRTE,
     +               CNVRTDE,IREORDER,KSDIAG,KEDIAG,KSOFFD,KEOFFD
c      COMMON /UTILCM/ DGSCARTNU(NATOM,3),DESCARTNU(NATOM,3,ISURF),
c     +                DIJCARTNU(NATOM,3,JSURF),CNVRTD,CNVRTE,CNVRTDE,
c     +                IREORDER,KSDIAG,KEDIAG,KSOFFD,KEOFFD
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
C**********************************************
C                                             *
C   The following code will convert units     *
C   and relabel the                           *
C   atoms to match the the numbering          *
C   convention chosen by the user             *
C                                             *
C   If the user wishes to relabel the         *
C   atoms set the variable IREORDER equal     *
C   to 1 in the PARAMETER statement.  The     *
C   user must also specify the new labeling   *
C   scheme.  This is done using the array     *
C   NULBL in the following manner:            *
C                                             *
C   NULBL(i) = j                              *
C                                             *
C      where:  i => old index                 *
C              j => new index                 *
C                                             *
C**********************************************
C
      IF (IREORDER.EQ.1) THEN
          DO I=1,NATOMS
             DO J=1,3
                CARTNU(NULBL(I),J)=CART(I,J)/CNVRTD
             END DO
          END DO
      ELSE
          DO I=1,NATOMS
             DO J=1,3
                CARTNU(I,J)=CART(I,J)/CNVRTD
             END DO
          END DO
      END IF
C
      RETURN
      END
 
C
      SUBROUTINE CARTTOR
C
C*****************************************************************
C                                                                *
C   WRITTEN BY: R. J. DUCHOVIC, A. F. WAGNER                     *
C                                                                *
C   DATE:               00-01-22                                 *
C                       00-08-03                                 *
C                       00-08-11                                 *
C                       01-07-18                                 * 
C                                                                *
C   THIS SUBROUTINE IS DESIGNED TO CONVERT FROM                  *
C   POSITIONS CARTNU IN DEVELOPER UNITS AND ORDER                *
C   TO DIFFERENT SETS OF INTERNAL COORDINATES R                  *
C   AS CONTROLLED BY ICARTR                                      *
C                                                                *
C   CURRENT ICARTR OPTIONS ARE:                                  *
C                                                                *
C   ICARTR = 0                                                   *
C     NOTHING HAPPENS (R VECTOR IS NOT FILLED)                   *
C                                                                *
C   ICARTR = 1                                                   *
C     ARRAY CARTNU WILL BE STORED IN THE VECTOR R AS             *
C      R(1) ....... R(3)  CARTESIANS (X,Y,Z) ATOM 1              *
C      R(4) ....... R(6)  CARTESIANS (X,Y,Z) ATOM 2              *
C      R(7) ....... R(9)  CARTESIANS (X,Y,Z) ATOM 3              *
C       .            .                                           *
C      R(3*N-2) ... R(3*N)  CARTESIANS (X,Y,Z) ATOM N            *
C       .            .                                           *
C      R(73) ...... R(75)  CARTESIANS (X,Y,Z) ATOM 25            *
C                                                                *
C   ICARTR = 2                                                   *
C     ARRAY CARTNU WILL BE CONVERTED INTO DISTANCE AS            *
C      R(1) .................... DISTANCE OF 1 TO 2              *
C      R(2) .................... DISTANCE OF 1 TO 3              *
C       .                        .                               *
C      R(NATOMS) ............... DISTANCE OF 1 TO NATOMS         *
C      R(NATOMS+1) ............. DISTANCE OF 2 TO 3              *
C      R(NATOMS+2) ............. DISTANCE OF 2 TO 4              *
C       .                        .                               *
C      R(2*NATOMS-1) ........... DISTANCE OF 2 TO NATOMS         *
C      R(2*NATOMS)   ........... DISTANCE OF 3 TO 4              *
C       .                        .                               *
C      R(NATOMS*(NATOMS-1)/2 ... DISTANCE OF NATOMS-1 TO NATOMS  *
C                                                                *
C   ICARTR = 3                                                   *
C     !!!FOR NATOMS = 3 ONLY !!!                                 *
C      R(1) .................... DISTANCE OF 1 TO 2              *
C      R(2) .................... DISTANCE OF 2 TO 3              *
C      R(3) .................... DISTANCE OF 1 TO 3              *
C                                                                *
C   ICARTR = 4                                                   *
C     Transformation for PES requiring Distances and Angles      *
C      HF Dimer PESs                                             *
C      R(1) .................... Center-of-Mass DISTANCE         *
C      R(2) .................... Bond DISTANCE First Diatom      *
C      R(3) .................... Bond DISTANCE Second Diatom     *
C      R(4) .................... Theta one angle                 *
C      R(5) .................... Theta two angle                 *
C      R(6) .................... Dihedral Tau angle              *
C                                                                *
C   NATOMS  =>  NUMBER OF ATOMS                                  *
C                                                                *
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      CHARACTER(75) REF(5)
C
      PARAMETER (N3ATOM=75)
      PARAMETER (NATOM=25)
      PARAMETER (ISURF=5)
C
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
C************************************************************
C                                                           *
C   SELECT CARTNU -> R OPTION WITH ICARTR                   *
C                                                           *
C************************************************************
C
      IF(ICARTR.EQ.1) THEN
C
C************************************************************
C                                                           *
C   R = VECTOR MAP OF CARTNU ARRAY                          *
C                                                           *
C************************************************************
C
         DO I=1,NATOMS
            IND=3*I-2
            R(IND)   = CARTNU(I,1)
            R(IND+1) = CARTNU(I,2)
            R(IND+2) = CARTNU(I,3)
         END DO
      ELSEIF(ICARTR.EQ.2) THEN
C
C************************************************************
C                                                           *
C   R = INTERPARTICLE DISTANCES IN CANONICAL ORDER          *
C                                                           *
C   R(1) 1-2 DISTANCE                                       *
C   R(2) 1-3 DISTANCE                                       *
C   ....                                                    *
C   R(2) 2-3 DISTANCE                                       *
C   ....                                                    *
C                                                           *
C************************************************************
C
         I = 1                                                          
         DO K=1,NATOMS-1
            DO L = K+1,NATOMS                                          
               R(I) = SQRT( (CARTNU(K,1)-CARTNU(L,1))**2 +
     +                      (CARTNU(K,2)-CARTNU(L,2))**2 +
     +                      (CARTNU(K,3)-CARTNU(L,3))**2 )
               I = I + 1                                              
            END DO
         ENDDO
      ELSEIF(ICARTR.EQ.3) THEN
C
C************************************************************
C                                                           *
C   !!! FOR NATOMS = 3 ONLY!!!                              *
C   R = INTERPARTICLE DISTANCES IN THE ORDER                *
C       R(1)=1-2 DIST., R(2)=2-3 DIST., R(3)=1-3 DIST.      *
C                                                           *
C************************************************************
C
         R(1) = SQRT( (CARTNU(1,1)-CARTNU(2,1))**2 +
     +                (CARTNU(1,2)-CARTNU(2,2))**2 +
     +                (CARTNU(1,3)-CARTNU(2,3))**2 )
         R(2) = SQRT( (CARTNU(2,1)-CARTNU(3,1))**2 +
     +                (CARTNU(2,2)-CARTNU(3,2))**2 +
     +                (CARTNU(2,3)-CARTNU(3,3))**2 )
         R(3) = SQRT( (CARTNU(1,1)-CARTNU(3,1))**2 +
     +                (CARTNU(1,2)-CARTNU(3,2))**2 +
     +                (CARTNU(1,3)-CARTNU(3,3))**2 )
C
      ELSEIF(ICARTR.EQ.4) THEN
C
C*************************************************************
C                                                            *
C   ICARTR = 4                                               *
C     Transformation for PES requiring Distances and Angles  *
C      HF Dimer PESs                                         *
C                                                            *
C*************************************************************
C
C      COMMON//V,XB(30),DXB(30)
C      COMMON/SFCOM/IFLAG1,ICOR,ISFCT
C
C*********************************************
C                                            *
C   THE CARTESIAN COORDINATES WERE READ IN   *
C   Calculate INTERNAL COORDINATES           * 
C                                            *
C   FIND CENTER-OF-MASS OF EACH MONOMER      *
C                                            *
C   MASS OF HYDROGEN AND OF FLUORINE         *
C                                            *
C*********************************************
C
      FLM=18.99840D0
      HYM=1.007825D0
C
C************************************************************************
C                                                                       *
C   CENTER OF MASS OF FIRST MONOMER, DO X, Y, AND Z COORDINATE IN TURN  *
C                                                                       *
C************************************************************************
C
C      XCM1=(HYM*XB(1)+FLM*XB(4))/(FLM+HYM)
C      YCM1=(HYM*XB(2)+FLM*XB(5))/(FLM+HYM)
C      ZCM1=(HYM*XB(3)+FLM*XB(6))/(FLM+HYM)
C
      XCM1=(HYM*CARTNU(1,1)+FLM*CARTNU(2,1))/(FLM+HYM)
      YCM1=(HYM*CARTNU(1,2)+FLM*CARTNU(2,2))/(FLM+HYM)
      ZCM1=(HYM*CARTNU(1,3)+FLM*CARTNU(2,3))/(FLM+HYM)
C
C****************************************
C                                       *
C   CENTER OF MASS OF SECOND MONOMER    *
C                                       *
C****************************************
C
C      XCM2=(HYM*XB(7)+FLM*XB(10))/(FLM+HYM)
C      YCM2=(HYM*XB(8)+FLM*XB(11))/(FLM+HYM)
C      ZCM2=(HYM*XB(9)+FLM*XB(12))/(FLM+HYM)
C
      XCM2=(HYM*CARTNU(3,1)+FLM*CARTNU(4,1))/(FLM+HYM)
      YCM2=(HYM*CARTNU(3,2)+FLM*CARTNU(4,2))/(FLM+HYM)
      ZCM2=(HYM*CARTNU(3,3)+FLM*CARTNU(4,3))/(FLM+HYM)
C
C*******************************************
C                                          *
C   CENTER-OF-MASS SEPARATION OF MONOMERS  *
C                                          *
C*******************************************
C
      XCM3=XCM2-XCM1
      YCM3=YCM2-YCM1
      ZCM3=ZCM2-ZCM1
C
C******************************************************
C                                                     *
C   DEFINE THE VECTOR FROM CENTER OF MASS OF FIRST    *
C   MONOMER (XCM1,YCM1,ZCM1) TO ATOM H1               *
C                                                     *
C******************************************************
C                        
C      XRM1=XB(1)-XCM1
C      YRM1=XB(2)-YCM1
C      ZRM1=XB(3)-ZCM1
C
      XRM1=CARTNU(1,1)-XCM1
      YRM1=CARTNU(1,2)-YCM1
      ZRM1=CARTNU(1,3)-ZCM1
C
C*************************
C                        *
C   DEFINE COS(THETA1)   *
C                        *
C*************************
C
      THETA1=(XRM1*XCM3+YRM1*YCM3+ZRM1*ZCM3)
      THETA1=THETA1/(SQRT(XRM1**2+YRM1**2+ZRM1**2))
      THETA1=THETA1/(SQRT(XCM3**2+YCM3**2+ZCM3**2))
      IF(THETA1.GT.1.0D0)THETA1=1.0D0
      IF(THETA1.LT.-1.0D0)THETA1=-1.0D0
      THETA1=ACOS(THETA1)
C
C******************************************************
C                                                     *
C   DEFINE THE VECTOR FROM CENTER OF MASS OF SECOND   *
C   MONOMER (XCM2,YCM2,ZCM2) TO ATOM H2               *
C                                                     *
C******************************************************
C
C      XRM2=XB(7)-XCM2
C      YRM2=XB(8)-YCM2
C      ZRM2=XB(9)-ZCM2
C
      XRM2=CARTNU(3,1)-XCM2
      YRM2=CARTNU(3,2)-YCM2
      ZRM2=CARTNU(3,3)-ZCM2
C
C*************************
C                        *
C   DEFINE COS(THETA2)   *
C                        *
C*************************
C
      THETA2=(XRM2*(-XCM3)+YRM2*(-YCM3)+ZRM2*(-ZCM3))
      THETA2=THETA2/(SQRT(XRM2**2+YRM2**2+ZRM2**2))
      THETA2=THETA2/(SQRT(XCM3**2+YCM3**2+ZCM3**2))
      IF(THETA2.GT.1.0D0)THETA2=1.0D0
      IF(THETA2.LT.-1.0D0)THETA2=-1.0D0
      THETA2=ACOS(THETA2)
      PI=ACOS(-1.0D0)
      THETA2=PI-THETA2
C
C*****************
C                *
C   DEFINE PHI   *
C                *
C*****************
C
      Q1=SQRT(XRM1**2+YRM1**2+ZRM1**2)
      Q2=SQRT(XRM2**2+YRM2**2+ZRM2**2)
      CMM=(XCM3**2+YCM3**2+ZCM3**2)
      CMM=SQRT(CMM)
C      HHD=(XB(1)-XB(7))**2+(XB(2)-XB(8))**2+(XB(3)-XB(9))**2
      HHD=(CARTNU(1,1)-CARTNU(3,1))**2 +
     +    (CARTNU(1,2)-CARTNU(3,2))**2 +
     +    (CARTNU(1,3)-CARTNU(3,3))**2
      HHD=SQRT(HHD)
      Q=CMM-Q1*COS(THETA1)+Q2*COS(THETA2)
      Q3=SQRT(ABS(HHD**2-Q**2))
      Q1=Q1*SIN(THETA1)
      Q2=Q2*SIN(THETA2)
      CPHI=(Q1**2+Q2**2-Q3**2)/(2.*Q1*Q2)
      IF(CPHI.LT.-1.0D0)CPHI=-1.0D0
      IF(CPHI.GT.1.0D0)CPHI=1.0D0
      PHI=ACOS(CPHI)
C
C
C*************************************************
C                                                *
C   ASSIGN INTERNAL COORDINATES TO THE ARRAY R   *
C                                                *
C*************************************************
C
C      RCM=SQRT(XCM3**2+YCM3**2+ZCM3**2)
C      R1=(SQRT(XRM1**2+YRM1**2+ZRM1**2))*(FLM+HYM)/FLM
C      R2=(SQRT(XRM2**2+YRM2**2+ZRM2**2))*(FLM+HYM)/FLM
C     WRITE(6,2001)RCM,R1,R2,THETA1,THETA2,PHI
C 2001 FORMAT(6F12.8)
C
C************************************************************
C                                                           *
C   Assign internal coordinates to elements of the R array  *
C                                                           *
C************************************************************
C
      R(1)=SQRT(XCM3**2+YCM3**2+ZCM3**2)
      R(2)=(SQRT(XRM1**2+YRM1**2+ZRM1**2))*(FLM+HYM)/FLM
      R(3)=(SQRT(XRM2**2+YRM2**2+ZRM2**2))*(FLM+HYM)/FLM
      R(4)=THETA1
      R(5)=THETA2
      R(6)=PHI
C
C      RETURN
C      END
C
      ELSEIF(ICARTR.NE.0) THEN
C
         WRITE(NFLAG(18),1000) ICARTR
 1000    FORMAT(2X,'WRONG ICARTR FOR CARTNU; ICARTR =',I5//)
         STOP
C
      ENDIF
      RETURN
      END
 
 
C
      SUBROUTINE EUNITZERO
C
C************************************************************
C                                                           *
C   WRITTEN BY: R. J. DUCHOVIC, A. F. WAGNER                *
C                                                           *
C   DATE:               00-01-23                            *
C                       01-07-18                            *
C                                                           *
C   THIS PROGRAM CONVERTS OUTPUT ENERGIES FROM              *
C   DEVELOPER UNITS (NAMELY, ATOMIC UNITS) TO USER UNITS.   *
C   CONVERSION CALCULATED IN ANCVRT.F.                      *
C                                                           *
C************************************************************
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C      CHARACTER(75) REF(5)
C
      PARAMETER (N3ATOM=75)
      PARAMETER (NATOM=25)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
 
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON /PT3CM/  EZERO(ISURF+1)
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)
      COMMON /PT5CM/  ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)
C
      COMMON/UTILCM/ DGSCARTNU(NATOM,3),DESCARTNU(NATOM,3,ISURF),
     +               DIJCARTNU(NATOM,3,JSURF),CNVRTD,CNVRTE,
     +               CNVRTDE,IREORDER,KSDIAG,KEDIAG,KSOFFD,KEOFFD
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      PENGYGS = ENGYGS * CNVRTE - ANUZERO
      IF(KSDIAG.NE.0) THEN
         DO I=KSDIAG,KEDIAG
            PENGYES(I) = ENGYES(I) * CNVRTE - ANUZERO
         END DO
      ENDIF
      IF(KSOFFD.NE.0) THEN
         DO J=KSOFFD,KEOFFD
            PENGYIJ(J) = ENGYIJ(J) * CNVRTE
         END DO
      ENDIF
C
      RETURN
      END
 
C
      SUBROUTINE RTOCART
C
C************************************************************
C                                                           *
C   WRITTEN BY: R. J. DUCHOVIC, A. F. WAGNER                *
C                                                           *
C   DATE:               00-01-23                            *
C                       00-08-04                            *
C                       00-08-11                            *
C                       00-08-14                            *
C                       01-07-18                            *
C                                                           *
C   THIS PROGRAM IS DESIGNED TO CONVERT                     *
C   DEGSDR,DEESDR,DEIJDR IN INTERNAL COORDINATES R INTO     *
C   DGSCARTNU,DESCARTNU,DIJCARTNU IN DEVELOPER UNITS AND    *
C   ORDER.                                                  *
C   THIS CODE IS ONLY CALLED IF DERIVATIVES ARE REQUIRED.   *
C   THIS CODE IS THE INVERSE OF CARTTOR.                    *
C                                                           *
C************************************************************
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      CHARACTER(75) REF(5)
C
      PARAMETER (N3ATOM=75)
      PARAMETER (NATOM=25)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON /PT3CM/  EZERO(ISURF+1)
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)
      COMMON /PT5CM/  ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)
C
      COMMON /UTILCM/ DGSCARTNU(NATOM,3),DESCARTNU(NATOM,3,ISURF),
     +                DIJCARTNU(NATOM,3,JSURF),CNVRTD,CNVRTE,CNVRTDE,
     +                IREORDER,KSDIAG,KEDIAG,KSOFFD,KEOFFD
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      DIMENSION YGS(N3ATOM),YES(N3ATOM,ISURF),YIJ(N3ATOM,JSURF)
C
C************************************************************
C                                                           *
C   SELECT R -> CARTNU OPTION WITH ICARTR                   *
C                                                           *
C************************************************************
C
      IF(ICARTR.EQ.1) THEN
C
C********************************************************
C                                                       *
C   CARTNU = ARRAY MAP OF VECTOR R                      *
C   Do in order Ground State, Excited State, Coupling   *
C                                                       *
C********************************************************
C
         DO I = 1, NATOMS
            IND=3*I-2
            DGSCARTNU(I,1) = DEGSDR(IND)
            DGSCARTNU(I,2) = DEGSDR(IND+1)
            DGSCARTNU(I,3) = DEGSDR(IND+2)
            IF(KSDIAG.NE.0) THEN
               DO J = KSDIAG,KEDIAG
                  DESCARTNU(I,1,J) = DEESDR(IND,J)
                  DESCARTNU(I,2,J) = DEESDR(IND+1,J)
                  DESCARTNU(I,3,J) = DEESDR(IND+2,J)
               END DO
            ENDIF
            IF(KEOFFD.NE.0) THEN
               DO K = KSOFFD,KEOFFD
                  DIJCARTNU(I,1,K) = DEIJDR(IND,K)
                  DIJCARTNU(I,2,K) = DEIJDR(IND+1,K)
                  DIJCARTNU(I,3,K) = DEIJDR(IND+2,K)
               END DO
            ENDIF
         END DO
C
      ELSEIF(ICARTR.EQ.2) THEN
C
C**************************************************************
C                                                             *
C   DERIVATIVE WRT CARTNU DERIVED FROM DV/DR VIA CHAIN RULE   *
C     WHERE R= INTERPARTICLE DISTANCES IN CANONICAL ORDER:    *
C   R(1)=1-2 DIST.,R(2)=1-3 DIST.,...,R(NATOMS) 2-3 DIST.,... *
C   Do in order Ground State, Excited State, Coupling         *
C                                                             *
C**************************************************************
C JC05132003
         DO I = 1, NATOMS 
            DGSCARTNU(I,1) = 0.D0
            DGSCARTNU(I,2) = 0.D0
            DGSCARTNU(I,3) = 0.D0
            IF(KSDIAG.NE.0) THEN
               DO J1=KSDIAG,KEDIAG
                  DESCARTNU(I,1,J1) = 0.D0
                  DESCARTNU(I,2,J1) = 0.D0
                  DESCARTNU(I,3,J1) = 0.D0
               ENDDO
            ENDIF
            IF(KSOFFD.NE.0) THEN
               DO J2=KSOFFD,KEOFFD
                  DIJCARTNU(I,1,J2) = 0.D0
                  DIJCARTNU(I,2,J2) = 0.D0
                  DIJCARTNU(I,3,J2) = 0.D0
               ENDDO
            ENDIF
            DO J = 1,NATOMS
               IF(J.LT.I) THEN
                  M1 = NATOMS*(J-1) - (J*(J-1))/2 + I-J
               ELSEIF(J.GT.I) THEN
                  M1 = NATOMS*(I-1) - (I*(I-1))/2 + J-I
               ELSE
                  GO TO 20
               ENDIF
               Y = DEGSDR(M1)
               TERMX = (CARTNU(I,1)-CARTNU(J,1))/R(M1)
               TERMY = (CARTNU(I,2)-CARTNU(J,2))/R(M1)
               TERMZ = (CARTNU(I,3)-CARTNU(J,3))/R(M1)
               DGSCARTNU(I,1) = DGSCARTNU(I,1) + TERMX*Y
               DGSCARTNU(I,2) = DGSCARTNU(I,2) + TERMY*Y
               DGSCARTNU(I,3) = DGSCARTNU(I,3) + TERMZ*Y
               IF(KSDIAG.GT.0) THEN
                  Y = DEESDR(M1,J1)
                  DO J1=KSDIAG,KEDIAG
                     DESCARTNU(I,1,J1)=DESCARTNU(I,1,J1) + TERMX*Y
                     DESCARTNU(I,2,J1)=DESCARTNU(I,2,J1) + TERMY*Y
                     DESCARTNU(I,3,J1)=DESCARTNU(I,3,J1) + TERMZ*Y
                  ENDDO
               ELSEIF(KSOFFD.GT.0) THEN
                  DO J2=KSOFFD,KEOFFD
                     Y = DEIJDR(M1,J2)
                     DIJCARTNU(I,1,J2)=DIJCARTNU(I,1,J2) + TERMX*Y
                     DIJCARTNU(I,2,J2)=DIJCARTNU(I,2,J2) + TERMY*Y
                     DIJCARTNU(I,3,J2)=DIJCARTNU(I,3,J2) + TERMZ*Y
                  ENDDO
               ENDIF
20             CONTINUE
            ENDDO
         ENDDO
C
      ELSEIF(ICARTR.EQ.3) THEN
C
C************************************************************
C                                                           *
C   !!!! FOR 3 ATOMS ONLY !!!                               *
C   DERIVATIVE WRT CARTNU DERIVED FROM DV/DR VIA CHAIN RULE *
C     WHERE R(1)=1-2 DIST.,R(2)=2-3 DIST.,R(3)=1-3 DIST.    *
C                                                           *
C   First do prepartory work                                *
C                                                           *
C************************************************************
C
         DO I = 1, NATOMS
            YGS(I) = DEGSDR(I)/R(I)
            IF(KSDIAG.NE.0) THEN
               DO J=KSDIAG,KEDIAG
                  YES(I,J) = DEESDR(I,J)/R(I)
               ENDDO
            ENDIF
            IF(KSOFFD.NE.0) THEN
               DO K=KSOFFD,KEOFFD
                  YIJ(I,K) = DEIJDR(I,K)/R(I)
               ENDDO
            ENDIF
         ENDDO
C
C****************************************************************
C                                                               *
C   Second, do in order Ground State, Excited State, Coupling   *
C                                                               *
C****************************************************************
C
         DO K = 1,3
            TERM12 = CARTNU(1,K)-CARTNU(2,K)
            TERM23 = CARTNU(2,K)-CARTNU(3,K)
            TERM13 = CARTNU(1,K)-CARTNU(3,K)
            DGSCARTNU(1,K) = TERM12*YGS(1) + TERM13*YGS(3)
            DGSCARTNU(2,K) =-TERM12*YGS(1) + TERM23*YGS(2)
            DGSCARTNU(3,K) =-TERM13*YGS(3) - TERM23*YGS(2)
            IF(KSDIAG.NE.0) THEN
               DO J1=KSDIAG,KEDIAG
                 DESCARTNU(1,K,J1) = TERM12*YES(1,J1) + TERM13*YES(3,J1)
                 DESCARTNU(2,K,J1) =-TERM12*YES(1,J1) + TERM23*YES(2,J1)
                 DESCARTNU(3,K,J1) =-TERM13*YES(3,J1) - TERM23*YES(2,J1)
               ENDDO
            ENDIF
            IF(KSOFFD.NE.0) THEN
               DO J2=KSOFFD,KEOFFD
                 DIJCARTNU(1,K,J2) = TERM12*YIJ(1,J2) + TERM13*YIJ(3,J2)
                 DIJCARTNU(2,K,J2) =-TERM12*YIJ(1,J2) + TERM23*YIJ(2,J2)
                 DIJCARTNU(3,K,J2) =-TERM13*YIJ(3,J2) - TERM23*YIJ(2,J2)
               ENDDO
            ENDIF
         ENDDO
C
      ELSEIF(ICARTR.EQ.4) THEN
C
C************************************************************
C                                                           *
C   ICARTR = 4                                              *
C                                                           *
C   Transformation for PES requiring Distances and Angles   *
C                                                           *
C   HF Dimer PESs                                           *
C                                                           *
C************************************************************
C                                                           *
C   CALCULATE THE ATOM-ATOM DISTANCES IN THE HF DIMER       *
C   ARBITRARILY SET ONE HF IN THE X-Z PLANE AND THE OTHER   *
C   ALONG THE X-AXIS, WITH THE FIRST HF'S CENTER OF MASS    *
C   AT THE ORIGIN. FIND THE CARTESIANS OF EACH INDIVIDUAL   * 
C   ATOM AND THEN FIND THE SEPARATIONS OF EACH PAIR OF      *
C   ATOMS AS DESIRED                                        *
C                                                           *     
C   FOR THE CODE BELOW, "A", "B", "C", AND "D"              *
C   REFER TO THE INDIVIDUAL ATOMS OF THE COMPLEX.           *
C                                                           *
C   A IS THE FIRST FLUORINE                                 *
C   B IS THE FIRST HYDRDOGEN                                *
C   C IS THE SECOND FLUORINE                                *
C   D IS THE SECOND HYDROGEN                                *
C                                                           *
C   RCM IS SEPARATION OF THE TWO CENTERS OF MASS            *
C   R1 IS FIRST HF BOND LENGTH                              *
C   R2 IS SECOND HF BOND LENGTH                             *
C   U1 IS COSINE OF THETA1                                  *
C   U2 IS COSINE OF THETA2                                  *
C   U3 IS COSINE OF PHI                                     *
C   SS1 IS SINE OF THETA1                                   *
C   SS2 IS SINE OF THETA2                                   *
C   SS3 IS SINE OF PHI                                      *
C   WH IS MASS OF HYDROGEN IN AMU'S                         *
C   WF IS MASS OF FLUORINE IN AMU'S                         *
C                                                           *
C************************************************************
C
      WH=1.007825D0
      WF=18.99840D0
C
      SUM=WH+WF
      EPS=WF/SUM
      EPSP=WH/SUM
C
      U1=COS(THETA1)
      U2=COS(THETA2)
      U3=COS(PHI)
C
      SS1=SIN(THETA1)
      SS2=SIN(THETA2)
      SS3=SIN(PHI)
C
      YA=0.0D0
      YB=0.0D0
C
      T0=R1*U1
      ZA=-EPSP*T0
      ZB=EPS*T0
C
      T0=R1*SS1
      XA=-EPSP*T0
      XBB=EPS*T0
C
      T0=R2*SS2
      T1=T0*U3
      XC=-EPSP*T1
      XD=EPS*T1
C
      T1=T0*SS3
      YC=-EPSP*T1
      YD=EPS*T1
C
      T0=R2*U2
      ZC=-EPSP*T0+RCM
      ZD=EPS*T0+RCM
C
      RFF=SQRT((XA-XC)**2+YC**2+(ZA-ZC)**2)
C
C      RETURN
C      END
C
      ELSE
C
          WRITE(NFLAG(18),1000) ICARTR
1000      FORMAT(2X,' WRONG ICARTR FOR DERIVATIVE; ICARTR =',I5//)
          STOP
C
      ENDIF
C
      RETURN
      END
 
C
      SUBROUTINE DEDCOU
C
C********************************************************************
C                                                                   *
C   WRITTEN BY: R. J. DUCHOVIC, A. F. WAGNER                        *
C                                                                   *
C   DATE:               00-01-23                                    *
C                       01-07-18                                    * 
C                                                                   * 
C   THIS PROGRAM CONVERTS DERIVATIVES                               *
C   DGSCARTNU,DESCARTNU,DIJCARTNU IN DEVELOPER ORDER AND UNITS INTO *
C   DCSCART,DESCART,DIJCART IN USER ODER AND UNITS                  *
C   (DEVELOPER UNITS ARE ATOMIC UNITS)                              *
C   CONVERSION CALCULATED IN ANCVRT.F                               *
C                                                                   * 
C********************************************************************
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      CHARACTER(75) REF(5)
C
      PARAMETER (N3ATOM=75)
      PARAMETER (NATOM=25)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      COMMON /UTILCM/ DGSCARTNU(NATOM,3),DESCARTNU(NATOM,3,ISURF),
     +                DIJCARTNU(NATOM,3,JSURF),CNVRTD,CNVRTE,CNVRTDE,
     +                IREORDER,KSDIAG,KEDIAG,KSOFFD,KEOFFD
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
C***********************************************************************
C                                                                      *
C   Do in order, derivatives of Ground State, Excited State, Couplings *
C                                                                      *
C***********************************************************************
C
      IF (IREORDER.EQ.1) THEN
         DO I = 1, NATOMS
            DGSCART(I,1) = DGSCARTNU(NULBL(I),1) * CNVRTDE
            DGSCART(I,2) = DGSCARTNU(NULBL(I),2) * CNVRTDE
            DGSCART(I,3) = DGSCARTNU(NULBL(I),3) * CNVRTDE
            IF(KSDIAG.NE.0) THEN
               DO J=KSDIAG,KEDIAG
                  DESCART(I,1,J) = DESCARTNU(NULBL(I),1,J) * CNVRTDE
                  DESCART(I,2,J) = DESCARTNU(NULBL(I),2,J) * CNVRTDE
                  DESCART(I,3,J) = DESCARTNU(NULBL(I),3,J) * CNVRTDE
               END DO
            ENDIF
            IF(KSOFFD.NE.0) THEN
               DO K=KSOFFD,KEOFFD
                  DIJCART(I,1,K) = DIJCARTNU(NULBL(I),1,K) * CNVRTDE
                  DIJCART(I,2,K) = DIJCARTNU(NULBL(I),2,K) * CNVRTDE
                  DIJCART(I,3,K) = DIJCARTNU(NULBL(I),3,K) * CNVRTDE
               END DO
            ENDIF
         END DO
      ELSE
         DO I = 1, NATOMS
            DGSCART(I,1) = DGSCARTNU(I,1) * CNVRTDE
            DGSCART(I,2) = DGSCARTNU(I,2) * CNVRTDE
            DGSCART(I,3) = DGSCARTNU(I,3) * CNVRTDE
            IF(KSDIAG.NE.0) THEN
               DO J=KSDIAG,KEDIAG
                  DESCART(I,1,J) = DESCARTNU(I,1,J) * CNVRTDE
                  DESCART(I,2,J) = DESCARTNU(I,2,J) * CNVRTDE
                  DESCART(I,3,J) = DESCARTNU(I,3,J) * CNVRTDE
               END DO
            ENDIF
            IF(KSOFFD.NE.0) THEN
               DO K=KSOFFD,KEOFFD
                  DIJCART(I,1,K) = DIJCARTNU(I,1,K) * CNVRTDE
                  DIJCART(I,2,K) = DIJCARTNU(I,2,K) * CNVRTDE
                  DIJCART(I,3,K) = DIJCARTNU(I,3,K) * CNVRTDE
               END DO
            ENDIF
         END DO
      ENDIF
C
      RETURN
      END
 
