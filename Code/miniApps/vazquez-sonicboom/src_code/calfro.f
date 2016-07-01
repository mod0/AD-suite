      SUBROUTINE CALFRO
c---------------------------------------------------------------------   
c Defines the set of physical boundary mesh vertices
c Renumerotes them in order to improve vectorization
c 1..nf1  : slipping vertices
c 1..nf2  : upstream vertices with enforced boundary conditions
c 1..nf3  : freestream vertices for which boundary conditions are 
c           treated through upwinding
c 1..nf11 : non-slipping vertices
c N.B: mesh vertices can be at the same time (non)-slipping 
c      and freestream boundary vertices
c---------------------------------------------------------------------   
      INCLUDE 'Param3D.h'
c---------------------------------------------------------------------   
c     Local variables definition
c
      INTEGER ifac, logf  , logic
      INTEGER is  , nsb2
      INTEGER noea(nfcmax), noe11(nfcmax), noe2(nfcmax), noe3(nfcmax)
c
c*** Falcon : logfac=200 => coque : faces glissantes
c             logfac=53 => faces glissantes
c             logfac=11 => conditions de Steger-Warming

c    Aile   : logfac= -2  => coque : faces glissantes           <<<<------  I always use this!
c             logfac=  2  => faces glissantes
c             logfac=  4  => conditions de Steger-Warming
c             logfac=  5  => conditions de Dirichlet

c    Tuyere : logfac=4 => conditions de Steger-Warming
c             logfac=2 => faces glissantes
c             logfac=7 => coque : faces glissantes
c
c     Initialisations
c
      nf1                          = 0
      nf2                          = 0
      nf3                          = 0
      nf11                         = 0
c
      DO 101 ifac=1,nfac
c
         logf                      = logfac(ifac)
         logic                     = 0
c
         IF (ivis .EQ. 0) THEN
c
            IF ((ABS(logf) .EQ. 2) .OR. (logf .EQ. 0) .OR.
     $  (logf.eq.200) .OR. (logf.eq.53) .OR. (logf.eq.7)) THEN
               logic = 1
            ENDIF
c
         ELSE
c
            IF ((ABS(logf) .EQ. 2) .OR. (logf .EQ. 0)) THEN 
c
               IF (coefm1 .EQ. 1) THEN
c
c                 Shock tube problem
c 
                  logic            = 1
c
               ELSE
c
                  IF (coefm1 .EQ. 2) THEN
c
c                    External flow around an M6 wing
c 
                     IF (logf .EQ. -2) logic =-1
                     IF (logf .EQ. 2)  logic = 1 
c
                  ELSE
c
                     logic         =-1
c
                  ENDIF
               ENDIF
            ENDIF   
         ENDIF
c
c*** Cas du Falcon, de l'aile ou de la tuyere
c
         IF ((coefm1.eq.5).OR.(coefm1.eq.10)) THEN
            if ((logf.eq.4) .OR. (logf.eq.11)) then
               logic = 4
            endif
            if ((logf.eq.5)) then
               logic = 5
            endif
         ENDIF
c
         IF (logic .EQ. -1) THEN
            nf11                   = nf11 + 1
            noe11(nf11)            = ifac 
         ENDIF
c
c*** Cas des faces glissantes
c
         IF (logic .EQ. 1) THEN
            nf1                    = nf1 + 1
            noea(nf1)              = ifac 
         ENDIF
c
c        Initializing the free stream uniform physical state
c        for external flows around a wing or a falcon.
c
         IF (((coefm1 .EQ. 2) .OR.
     .     (coefm1 .EQ. 5) .OR. (coefm1.EQ.10))
     .     .AND.
     .     ((logic .EQ. 4) .or. (logic .EQ. 5))) THEN
c
           if (logic.eq.4) then
             nf3 = nf3 + 1                                  ! Ste-Warm
             noe3(nf3)              = ifac
             ub3(1,nf3)             = roin/rhoref
             ub3(2,nf3)             = uxin/vref
             ub3(3,nf3)             = uyin/vref
             ub3(4,nf3)             = uzin/vref
             ub3(5,nf3)             = pin/pref           
           else if (logic.eq.5) then
             nf2 = nf2 + 1                                  ! Diric
             noe2(nf2)              = ifac
             ub3(1,nf3)             = roin/rhoref
             ub3(2,nf3)             = uxin/vref
             ub3(3,nf3)             = uyin/vref
             ub3(4,nf3)             = uzin/vref
             ub3(5,nf3)             = pin/pref           
           end if
c
c
         ENDIF
c
101   CONTINUE
c

      
C [llh]      WRITE(6,*) '              3- Les conditions aux bords :'
C [llh]      WRITE(6,*) '              ***************************** '
C [llh]      WRITE(6,*) ' '
C [llh]      WRITE(6, *) 'Number of non-slipping faces              : ', nf11
C [llh]      WRITE(6, *) 'Number of slipping faces                  : ', nf1
C [llh]      WRITE(6, *) 'Number of Dirichlet contition faces       : ', nf2
C [llh]      WRITE(6, *) 'Number of Steger-Warming condition faces  : ', nf3
c      WRITE(6, *) 'Number of upstream enforced faces         : ', nf2
c      WRITE(6, *) 'Number of freestream faces with upwinding : ', nf3
C [llh]      WRITE(6, *) ' '
C [llh]      WRITE(6, *) ' '
c
      IF (nfac .NE. (nf11 + nf1 + nf2 + nf3)) THEN
         WRITE(6, *) 'Error in CALFRO'
         CALL TILT
      ENDIF
c
      IF (nf11 .NE. 0) THEN
         DO 111 ifac=1,nf11
            noe1(ifac)             = noe11(ifac)
111      CONTINUE
      ENDIF
c
      IF (nf1 .NE. 0) THEN
         DO 112 ifac=1,nf1
            noe1(nf11+ifac)        = noea(ifac) 
112      CONTINUE
      ENDIF
c
      IF (nf2 .NE. 0) THEN
         DO 113 ifac=1,nf2
            noe1(nf11+nf1+ifac)    = noe2(ifac) 
113      CONTINUE
      ENDIF
c
      IF (nf3 .NE. 0) THEN
         DO 114 ifac=1,nf3
            noe1(nf11+nf1+nf2+ifac)= noe3(ifac) 
114      CONTINUE
      ENDIF
c
      RETURN
      END
