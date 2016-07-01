c     -----------------------------------------------------------------
c     Parameter and common definitions for the 3D Navier-Stokes solver 
c     -----------------------------------------------------------------
      IMPLICIT NONE
c     -----------------------------------------------------------------
c     Degres de liberte
      INTEGER NEQUATION
      PARAMETER (NEQUATION=5)
c
c     Maximum number of vertices 
      INTEGER nsmax
      PARAMETER (nsmax   = 12966)
c     Maximum number of tetraedra
      INTEGER ntmax
      PARAMETER (ntmax   = 65421)
c     Maximum number of boundary faces
      INTEGER nfcmax
      PARAMETER (nfcmax  = 7350)
c     Maximum number of edges 
      INTEGER nsgmax
      PARAMETER (nsgmax = 82056)
c     Maximum number of edges attached to a mesh vertex
      INTEGER ndmax
      PARAMETER (ndmax   = 150)
c     Maximum number of tetraedra attached to a mesh vertex
      INTEGER nvmax
      PARAMETER (nvmax   = 49)
c     Allowable values for boundary vertices logical reference
      INTEGER klogmin, klogmax
      PARAMETER (klogmin = 1, klogmax = 30)
c     Maximum numbers of colors for the tetraedra and edges 
      INTEGER nctmax, ncamax 
      PARAMETER (nctmax  = 710, ncamax  = 710)
c     Characteristic length for vectorisation
      INTEGER lvect
      PARAMETER (lvect   = 128)
c     Maximum length of arrays for vector operations
      INTEGER lvmax
      PARAMETER (lvmax   = 128)
c     -----------------------------------------------------------------
c     Maximum number of timers
      INTEGER icpmax
      PARAMETER (icpmax  = 21)
c     Timers definition
      REAL*8 time(icpmax)
c     Mflop rates
      REAL*8 flop(icpmax)
c
      COMMON/rtimrs/time, flop
c
c     -----------------------------------------------------------------
c     Effective numbers of tetraedra and edges colors 
      INTEGER nca, nct
c     Data structures for the tetraedra and edges colors
      INTEGER icolt(0:nctmax), icola(0:ncamax)
      INTEGER mark(nsgmax)
c
      COMMON/icolor/nca, nct, icolt, icola, mark
c
c     -----------------------------------------------------------------
c     Effective numbers of vertices, tetraedra, edges and 
c     boundary faces
      INTEGER ns, nt, nseg, nfac
c     Tetraedra connectivity table
      INTEGER nu(4,ntmax)
c     Boundary faces connectivity table
      INTEGER nsfac(3,nfcmax)
c     Boundary faces logical reference table
      INTEGER logfac(nfcmax),log1fac(nfcmax) 
c     Boundary vertices logical reference table
      INTEGER logfr(nsmax)
c     Auxiliary table for boundary vertices logical reference 
      INTEGER fv(nsmax)
c     Boundary faces index table
      INTEGER noe1(nfcmax)
c     Effective numbers of faces for each allowable value of
c     boundary face logical reference
      INTEGER nblog(klogmin:klogmax)
c     Edges connectivity table
      INTEGER nubo(2,nsgmax)
c     Effective number of edges attached to a mesh vertex
      INTEGER ndeg(nsmax), inew(nsmax)
c     Identification of the set of edges attached to 
c     a mesh vertex
      INTEGER jaret(nsmax,ndmax)
c     Effective number of tetraedra attached to a mesh vertex
      INTEGER nbvoi(nsmax)
c     Identification of the set of tetraedra attached to 
c     a mesh vertex
      INTEGER ivoi(nsmax,nvmax)
c     Identification of the upstream and downstream tetraedra 
c     associated to an edge
      INTEGER jta(2,nsgmax)
c     Effective number of faces with a given boundary logical 
c     reference 
      INTEGER nf1, nf2, nf3, nf11
c   
      COMMON/igeom0/ns, nt, nseg, nfac
      COMMON/igeom1/nu, nsfac, logfac, nubo, logfr, log1fac
      COMMON/igeom2/noe1, nblog, fv
      COMMON/igeom3/nf1, nf2, nf3, nf11
      COMMON/igeom4/ndeg, inew, jaret, nbvoi, ivoi, jta
c
c     Initial and instantaneous coordinates of the mesh vertices
      REAL*8 coor(3,nsmax), coco(3,nsmax) 
c     Initial mesh coordinates (zz)
      real*8 coin(3,nsmax)
c     Control volume and tetraedra volumes
      REAL*8 vols(nsmax), volt(ntmax)
c     Control volume boundary normal components
      REAL*8 vnocl(3,nsgmax) 
c     Boundary faces normal components
      REAL*8 vnfac(3,nfcmax)
c     Mean values of the normal to tetradra faces
      REAL*8 tvno(ntmax,4,3)
c 
      COMMON/rgeom1/coor, coco,coin
      COMMON/rgeom2/vols, volt, vnocl, vnfac 
      COMMON/rgeom3/tvno
c
c     -----------------------------------------------------------------
c     Physical solution
      REAL*8 ua(5,nsmax), un(5,nsmax)
c     Nodal fluxes (gathered convective and diffusive fluxes)
      REAL*8 ce(5,nsmax)
c     Nodal gradients
      REAL*8 dx(5,nsmax), dy(5,nsmax), dz(5,nsmax)
c     Local time steps
      REAL*8 dtl(nsmax)
c     Mesh vertices velocities
      REAL*8 xw1(nsmax), xw2(nsmax), xw3(nsmax)
c     Mean value of mesh edges and boundary faces velocities
      REAL*8 sigma(nsgmax), vnsig(nfcmax)
c     Auxiliary array used during mesh vertices based operations
      REAL*8 ce1(nsmax)
c
      COMMON/rsol1/ua, un, ce, dx, dy, dz
      COMMON/rsol2/dtl
      COMMON/rsol3/xw1, xw2, xw3
      COMMON/rsol4/sigma, vnsig
      COMMON/rsol5/ce1
c
c     -----------------------------------------------------------------
c     Reinitialisation flag 
      INTEGER ncont
c     Local time-step strategy flag
      INTEGER iloc
c     Courant number law parameters
      REAL*8 cfl, xcfl, ycfl, zcfl, cflmax
c     Convective flux solver identifier
      INTEGER iflux
c     Spatial approximation precision order
      INTEGER nordre
c     Predictor-corrector time integration scheme flag
      INTEGER ipred
c     Maximum, initial and effective number of time steps
      INTEGER ktmax, kt0, kt
c     Flag for a viscous computation
      INTEGER ivis
c     Flag for a moving mesh computation
      INTEGER idefor
c     Flag used when deforming the mesh
      INTEGER imcas
c     Flag for explicit/implicit time integration 
c     Maximum and effective number of relaxations in the 
c     implicit phase
c     Relaxations pour l'etat-adjoint
      INTEGER nexp, nbrel, nit, nbrelPi
c     Linear system solver identifier
      INTEGER irlax
c     Flag for the type of storage used for the extra-diagonal 
c     blocks 
      INTEGER istok
c     Residual tolerance for the linear iteration (implicit phase)
      REAL*8 err
c     Maximal, initial and effective  physical time 
      REAL*8 tmax, t0, t
c     Global and maximal allowable time step
      REAL*8 dt, ddtmax
c     Residual tolerance for the non-linear iteration
      REAL*8 resf
c     Flow characteristic quantities
      REAL*8 rhoref, pref, lref, vref
c     Free stream Mach number
      REAL*8 xmach
c     Free stream quantities used to define the free stream 
c     physical states
      REAL*8 roin , uxin , uyin , uzin , pin
      REAL*8 roout, uxout, uyout, uzout, pout
c     Free stream physical states
      REAL*8 ub2(5,nfcmax), ub3(5,nfcmax)
c     Real gaz physical quantities
      REAL*8 gam, gam1
c     Reynolds and Prandtl numbers, viscosisty and temperature 
c     on the body
      REAL*8 rey, pr, xmu, tbrd
c
      INTEGER oneshot
      COMMON/ressimul/oneshot
c
c     Nombre d'iterations de jacobi pour resoudre l'etat et l'etat-adjoint
c     en one-shot
      INTEGER NBRELOS
c
      COMMON/icmflg/ncont, iloc, ivis, idefor, ipred
      COMMON/icmflx/iflux
      COMMON/icmgrd/imcas
      COMMON/icmimp/nexp, nbrel, nit, irlax, istok, nbrelos, nbrelPi
      COMMON/icmord/nordre
c
c     Nombre d'iterations maximales pour l'etat en one-shot
      INTEGER ktmaxos

      COMMON/icmstp/ktmax, kt0, kt, ktmaxos
c
c     Residu a atteindre pour resoudre l'etat-adjoint
      Real*8 errjacPi
c
      COMMON/rcmstp/tmax, t0, t, resf
      COMMON/rcmdtg/dt, ddtmax
      COMMON/rcmimp/err,errjacPi
      COMMON/comcfl/cfl, xcfl, ycfl, zcfl, cflmax
      COMMON/comref/rhoref, pref , lref , vref
      COMMON/frsmch/xmach
      COMMON/frstr1/roin , uxin , uyin , uzin , pin , 
     &              roout, uxout, uyout, uzout, pout 
      COMMON/frstr2/ub2, ub3
      COMMON/comgam/gam, gam1
      COMMON/comvis/rey, pr, xmu, tbrd
c
c     -----------------------------------------------------------------
c     Frequency parameters used when printing and saving the physical 
c     solution
      INTEGER ifre, ifre1, impre, ifrel
c
      COMMON/icmedt/ifre, ifre1, impre, ifrel
c
c     -----------------------------------------------------------------
c     Selected problem identifier
      INTEGER coefm1
c     Variables related to the piston engine geometry
      REAL*8 cadmax, cads
      REAL*8 thet0, ca, squish, rps, stroke, conrod
      REAL*8 wpist, wspa, wspe
      REAL*8 vpist(3), vspa(3), vspe(3)
      REAL*8 aoa, aoe, rfa, rfe, lma, lme
      REAL*8 hpist, hspa, hspe, hplan, href
c
      COMMON/icmpst/coefm1
c
      COMMON/cmpst0/cadmax, cads
      COMMON/cmpst1/thet0, ca, squish, rps, stroke, conrod
      COMMON/cmpst2/wpist, wspa, wspe, vpist, vspa, vspe
      COMMON/cmpst3/aoa, aoe, rfa, rfe, lma, lme
      COMMON/cmpst4/hpist, hspa, hspe, hplan, href      
c
c     -----------------------------------------------------------------
c     Parameters and data structures related to the implicit time 
c     integration scheme
      INTEGER nsegs
      PARAMETER (nsegs   = 50*nsgmax)
c     Storage of the diagonal blocks
      REAL*8 diag(nsmax,5,5), adh1(nsmax)
c     Storage of the extra diagonal blocks
      REAL*8 stmat(nsegs)
c     Stockage des termes extra-diagonaux de la matrice
      Real*8 ZM(5,5,2,nsgmax)
c
      COMMON/comimp/diag, stmat, adh1, zm
c 
c     -----------------------------------------------------------------
c     Auxiliary variables used when computing the convective fluxes
      INTEGER ient, icoef
      REAL*8 epsiim, epsiex
c
      COMMON/icmfx1/ient, icoef
      COMMON/rcmfx1/epsiim, epsiex
c
c     -----------------------------------------------------------------
c     Auxiliary  variables related to the moving mesh feature
c     Maximum and effective number of relaxations in the 
c     moving mesh phase
      INTEGER maxjac, kjac
c     Residual tolerance and effective residual for the linear 
c     iteration in the moving mesh phase phase
      REAL*8 rsjacf, rsjac
c     Elastic center coordinates
      REAL*8 xce, yce, zce
c     Instantaneous angle of attack
      REAL*8 tetaT 
c     Original coordinates of the mesh vertices
      REAL*8 xt0(nsmax), yt0(nsmax), zt0(nsmax)
c     Displacements of the mesh vertices
      REAL*8 deltx(2,nsmax), delty(2,nsmax), deltz(2,nsmax)
c     Predicted displacements of the mesh vertices
      REAL*8 deltxp(nsmax), deltyp(nsmax), deltzp(nsmax)  
c
      COMMON/imomsh/maxjac, kjac
      COMMON/mvmsh0/rsjacf, rsjac
      COMMON/mvmsh1/xce, yce, zce, tetaT
      COMMON/mvmsh2/xt0, yt0, zt0
      COMMON/mvmsh3/deltx, delty, deltz, deltxp, deltyp, deltzp
c
      REAL*8 som1, dro1, som, dro
      COMMON/residus/som1, dro1, som, dro
c
      REAL*8 omega
      COMMON/Newton/omega
c
      INTEGER diric, bound
      COMMON/boundaries/diric, bound
c
C     NODE3D2D :Tab. de correspondance entre maillage et coque sur les noeuds
c     FACE3D2D :  "            "         "      "          "           faces
c
      INTEGER NODE3D2D(NSMAX), FACE3D2D(NFCMAX)
      COMMON/CORRES3D2D/NODE3D2D, FACE3D2D
C
      REAL*8 PDESP(NSMAX), DJDW(5,NSMAX), PIADJ(5,NSMAX)
      COMMON/PRESSION/PDESP
      COMMON/OPTI/DJDW,PIADJ
c
      integer validpsi, Newton, testder, testadj, testgrad,testderW
      common/verif/validpsi,Newton,testder,testadj,testgrad,testderW
c
      INTEGER ITMAX,ITOPT,NCF,NPRES, NCONTOPT, CONTR, defortuy
      COMMON/OPTIMIS/ITMAX,ITOPT,NCF,NPRES, NCONTOPT, CONTR, defortuy
C
      INTEGER NIVEAU,NIVG,MULTINIV,VCYC,FVC,NBVC,NBVCFIN
      COMMON/HIERA/NIVEAU,NIVG,MULTINIV,VCYC,FVC,NBVC,NBVCFIN
C
      REAL*8 ERREUR, ROOPT, ROINIT, ROOS
      COMMON/RESOPTI/ERREUR, ROOPT, ROINIT, ROOS
c
      REAL*8 coeftrainee, teta, cltarget,cdinit, cdtarget, tetacdcl
      REAL*8 coefpres
      real*8 tetaxy,tetaxz,tetazy
      integer*4  igicc
      common/iangles/igicc
      COMMON/calctrainchoc/coeftrainee,cltarget,cdinit,cdtarget
      COMMON/anglesincid/teta,tetacdcl,coefpres,tetaxy,tetaxz,tetazy

c
c     ZZ PARAMETERS:
c     ===============
c
      integer intzz(30)
      common/integerzz/intzz

c     CPU time evaluation (all real*4)

      real*8 pmnew,pmold
      common/mxmi/pmnew,pmold

C [llh Nan]      real*4 ct0,ctini,cttot,ctoto,cteto,ctato,ctjac,ctjad,ctjet
C [llh Nan]      common/cputimes/ct0,ctini,cttot,ctoto,cteto,ctato,ctjac,
C [llh Nan]     .                ctjad,ctjet
C [llh Nan]      real*4 ct1,ct2,ct3,ct4,ct5,ctf

c     -----------------------------------------------------------------
