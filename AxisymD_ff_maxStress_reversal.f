      subroutine vumat(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock)
C
      character*80 cmname
C      
      integer readflag
      real*4  edgedim(350000,3)
      common /integerbuf/readflag 
      common / realbuf/edgedim

        
        IF ( cmname(1:8).eq.'MATFIBRE' ) THEN
C		
            CALL  vumatfibre( 
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
	 
		ELSEIF ( cmname(1:6).eq.'MATCOH' ) THEN
 
            CALL vumatcoh(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
     
		ENDIF
C         
      RETURN
      END       

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C
C  subroutine: vumatfibre                                            C
C  function: Fibre failure activation using maximum stress criteria  C
C            and evolution                                           C 
C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C SDV1  - Strain 11
C SDV2  - Strain 22
C SDV3  - Strain 33
C SDV4  - Strain 12, engineering shear strain gamma
C SDV5  - Temporary strain 12
C SDV6  - Element deformation gradient, F
C SDV7  - Failure index under tension, FIt
C SDV8  - Failure index under compression, FIc
C SDV9  - Damage variable, D
C SDV10 - Element Delete flag; 0-delete
C SDV11 - Stress 11 at fibre damage initiation under tension, sig_f0_T     
C SDV12 - Strain 11 at fibre damage initiation under tension, eps_f0_T
C SDV13 - Strain 11 at fibre damage completion under tension, eps_ff_T
C SDV14 - Characteristic length, L
C SDV15 - Fibre failure falg: 0-elastic, 1-softening
C SDV16 - Old d1Plus
C SDV17 - Old D1Minus
C SDV18 - Stress 11 at fibre damage initiation under compression, sig_f0_C
C SDV19 - Strain 11 at fibre damage initiation under compression, eps_f0_C
C SDV20 - Strain 11 at fibre damage completion under compression, eps_ff_C
C SDV21 - Fibre failure flag under tension: 0-elastic, 1-softening
C SDV22 - Fibre failure flag under compression: 0-elastic, 1-softening

		 subroutine vumatfibre(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock)
C
       parameter(zero = 0.d0, one = 1.d0, two = 2.d0, three = 3.d0, 
     * third = one/three, half = 0.5d0, twothird = two/three,
     * threehalf = 1.5d0, safety = 1.d-20, tolSgn = 1.d-20, S_res=50.d0)              
C     
       character*80 cmname,cpname
       character*256 outdir,fullpath
       integer intnum,locnum,jrcd,lenoutdir
       double precision E11,E22,E33,v12,v21,v13,v31,v23,v32,G12,G13,G23
       double precision C11 ,C12 ,C13 ,C14 ,C21 ,C22,
     *   C23 , C24,C31 ,C32 ,C33 ,C34,
     *   C41 , C42 ,C43 ,C44
	   double precision n,m,delta
	   double precision Xt,Xc,Yc,Yt,Sl
	   double precision Gft,Gfc
       double precision sgn2,sgn3,sgn4,A,B
       double precision FIt, FIc, FIt_old, FIc_old, FIt_max, FIc_max
       double precision sig_f0, eps_f0, eps_f
	   double precision sig_f0_T, eps_f0_T, eps_ff_T, eps_f_T
	   double precision sig_f0_C, eps_f0_C, eps_ff_C, eps_f_C
	   double precision D, d1Plus, d1Minus
       real*8 Stress(ndir+nshr), Strain(ndir+nshr)
C	   
	   integer i,j
C
C INITIALISATION: MATERIAL CARD PARAMETERS
C =========================================
C
      E11 = props(1) !E11 = 61646
      E22 = props(2) !E33 = 13368
      E33 = props(3) !E22 = 61646
      v12 = props(4) !V13 = 0.3070
      v13 = props(5) !V12 = 0.3187
      v23 = props(6) !V32 = 0.0667 
      G12 = props(7) !G13 = 4575
      G13 = props(8) !G12 = 23373
      G23 = props(9) !G32 = 4575  
C	  
	  A = props(10) ! A = 145
	  B = props(11) ! B = 38
C	 
	  Xt  = props(12)  
	  Xc  = props(13)  
      Yt  = props(14)  
      Yc  = props(15)  
      Sl  = props(16)  
C	  
      Gft = props(17)  
	  Gfc = props(18)  
C
      v21 = v12*E22/E11
      v31 = v13*E33/E11
      v32 = v12
C
C AXISYMMETRIC MODEL STIFFNESS MATRIX 
C ====================================
C
      n = E11/E22
	  m = G12/E22
C
	  delta = (one + v13)*(one - v13 - two*n*v21**2)
C
	  C11 = (E22/delta) * n * (one - n*v21**2)
	  C12 = (E22/delta) * n * v21 * (one + v13)
	  C13 = (E22/delta) * n * (v13 + n* v21**2)
	  C14 = zero
	  C21 = C12
	  C22 = (E22/delta) * (one - v13**2)
	  C23 = (E22/delta) * n * v21 * (one + v13)
	  C24 = zero
	  C31 = C13
	  C32 = C23
	  C33 = (E22/delta) * n * (one - n*v21**2)
	  C34 = zero
	  C41 = zero
	  C42 = zero
	  C43 = zero
	  C44 = G12
	  
C	  
C START CONTINNUM DAMAGE MODEL
C ============================
C
C INITIAL ELASTIC STEP FOR ABAQUS TESTS
C -------------------------------------	  
C        
      IF (totalTime.eq.zero) THEN
      
          DO i = 1, nblock 
		  
		      stressNew(i,1)= stressOld(i,1)+ C11*straininc(i,1)
     *           + C12*straininc(i,2) + C13*straininc(i,3)
        
		      stressNew(i,2)= stressOld(i,2)+ C21*straininc(i,1)
     *           + C22*straininc(i,2) + C23*straininc(i,3) 

              stressNew(i,3)=stressOld(i,3) + C31*straininc(i,1)
     *           + C32*straininc(i,2) + C33*straininc(i,3)
	 
              stressNew(i,4) = stressOld(i,4) + two*C44*straininc(i,4)
		  
              D = stateNew(i,9)
		      stateNew(i,15) = zero ! SDV15 - Fibre failure falg: 0-elastic, 1-softening
			  stateNew(i,21) = zero
			  stateNew(i,22) = zero
              stateNew(i,10) = one  ! SDV10 - Element Delete flag; 0-delete
              stateNew(i,14) = charLength(i) 
		  
          ENDDO
C
C IF NOT INITIAL STEP, CAL. STRESSES ACCRODING TO CDM MODEL
C =========================================================
C      
      ELSE   
C         
          DO i = 1, nblock 
C
		      Do j = 1, nstatev
			      stateNew(i,j)=stateOld(i,j)
		      ENDDO
C		  
		      D = stateNew(i,9)
C
C DEGRADATION OF POISSON'S RATIO
C ------------------------------
			  ! v12=(1-D)*v12
			  ! v13=(1-D)*v13
			  
			  ! v21 = v12*E22/E11
			  ! v31 = v13*E33/E11
			  ! v32 = v12
			  
			  ! n = E11/E22
			  ! m = G12/E22
! C
			  ! delta = (one + v13)*(one - v13 - two*n*v21**2)
! C
			  ! C11 = (E22/delta) * n * (one - n*v21**2)
			  ! C12 = (E22/delta) * n * v21 * (one + v13)
			  ! C13 = (E22/delta) * n * (v13 + n* v21**2)
			  ! C14 = zero
			  ! C21 = C12
			  ! C22 = (E22/delta) * (one - v13**2)
			  ! C23 = (E22/delta) * n * v21 * (one + v13)
			  ! C24 = zero
			  ! C31 = C13
			  ! C32 = C23
			  ! C33 = (E22/delta) * n * (one - n*v21**2)
			  ! C34 = zero
			  ! C41 = zero
			  ! C42 = zero
			  ! C43 = zero
			  ! C44 = G12
C -------------------------------	  
C		  
		      IF ( stateNew(i,10) .eq. one) THEN  ! If the element is deleted
C			  
		          stateNew(i,1)= stateOld(i,1)+ straininc(i,1)
                  stateNew(i,2)= stateOld(i,2)+ straininc(i,2)
                  stateNew(i,3)= stateOld(i,3)+ straininc(i,3)
                  stateNew(i,4)= stateOld(i,4)+ two*straininc(i,4)  
C       
                  stressNew(i,1)= C11*stateNew(i,1)+C12*stateNew(i,2)
     *                  +C13*stateNew(i,3)
                  stressNew(i,2)= C21*stateNew(i,1)+C22*stateNew(i,2)
     *                 + C23*stateNew(i,3)
                  stressNew(i,3)= C31*stateNew(i,1)+C32*stateNew(i,2)
     *                 + C33*stateNew(i,3)  
C		  
				  sgn4 = stateNew(i,4)/(abs(stateNew(i,4)+safety))
				  IF(abs(stateNew(i,4)).gt.stateOld(i,5))THEN
					stateNew(i,5)=abs(stateNew(i,4))
					stressNew(i,4)=sgn4*(A*(one -exp(-B*stateNew(i,5))))  
				  ELSE
					stateNew(i,5)=stateOld(i,5)
					stressNew(i,4)=sgn4*(A*(one -exp(-B*stateNew(i,5))) 
     *                  -G12*(stateNew(i,5)-abs(stateNew(i,4))))
				  ENDIF
C		  
		          Stress(1:4) = stressNew(i,1:4) ! Effective stresses in the undamaged material				  
C
C FIBRE DAMAGE STATRS HERE
C ========================
		          stateNew(i,14) = charLength(i)
C
                  IF (stateNew(i,15) .eq. zero ) THEN  ! In the elastic region, calculate failure indices
C
C CALCULATE FAILURE INDICES
C =========================
				      FIt = zero ! Initialise failure indices
					  FIc = zero
C
C CALCULATE FAILURE INDECES AND RECORD THE MAX
C---------------------------------------------
					  IF ( Stress(1) >= zero ) THEN
						  FIt = Stress(1)/Xt     ! Max. stress criteria
					  ELSEIF ( Stress(1) < zero ) THEN
						  FIc = -Stress(1)/Xc    ! Max. stress criteria
					  ENDIF
C				  
					  FIt_old = stateOld(i,7)     ! Load old state variables
					  FIc_old = stateOld(i,8)
C				  
					  FIt_max = max(FIt, FIt_old) ! Compute max failure indices
					  FIc_max = max(FIc, FIc_old)
C				  
					  stateNew(i,7) = FIt_max     ! Record max failure indices for no damage reversal
					  stateNew(i,8) = FIc_max  
C
C CALCULATE DAMAGE VARIABLES
C ==========================
					  d1Plus = zero ! Initialise damage variables
					  d1Minus = zero
C					  
C LONGITUDINAL TENSILE DAMAGE
C ---------------------------
					  IF ( FIt_max >= one .or. FIc_max >= one ) THEN
C					      
						  stateNew(i,15) = one ! Fibre soften
C						  
						  ! sig_f0 = Stress(1)     ! Stress11 at damage initiation
						  ! stateNew(i,11) = sig_f0  ! SDV11
						  ! eps_f0 = stateNew(i,1)   ! Strain11 at damage initiation
						  ! stateNew(i,12) = eps_f0  ! SDV12
C						  
					      IF ( stateNew(i,1) > zero ) THEN
							  stateNew(i,21) = one
							  sig_f0_T = Stress(1)
							  stateNew(i,11) = sig_f0_T
							  eps_f0_T = stateNew(i,1)
							  stateNew(i,12) = eps_f0_T
						      eps_ff_T = 2.d0*Gft/(Xt*charLength(i)) ! Strain11 at final damage
						      stateNew(i,13) = eps_ff_T							  
						  ELSEIF ( stateNew(i,1) < zero ) THEN
							  stateNew(i,22) = one
							  sig_f0_C = Stress(1)
							  stateNew(i,18) = sig_f0_C
							  eps_f0_C = stateNew(i,1)
							  stateNew(i,19) = eps_f0_C
							  eps_ff_C = -2.d0*Gfc/(Xc*charLength(i)) ! eps_ff_C is minus
							  stateNew(i,20) = eps_ff_C							  
						  ENDIF
C						  
					  ENDIF
C					  
				  ELSEIF ( stateNew(i,15) .eq. one ) THEN
C
					  eps_f = stateNew(i,1)
C
					  IF ( stateNew(i,1) .GT. zero ) THEN
C					  
						  IF ( stateNew(i,21) .EQ. zero ) THEN
C
						      FIt = Stress(1)/Xt
C
							  IF ( FIt .GT. one ) THEN  !If tension damage initiates
							      stateNew(i,21) = one								  
								  sig_f0_T = Stress(1)
								  stateNew(i,11) = sig_f0_T
								  eps_f0_T = stateNew(i,1)
								  stateNew(i,12) = eps_f0_T
								  eps_ff_T = 2.d0*Gft/(Xt*charLength(i)) ! Strain11 at final damage
								  stateNew(i,13) = eps_ff_T
							  ENDIF
C							  
						  ELSEIF ( stateNew(i,21) .EQ. one ) THEN
							  eps_f_T = stateNew(i,1)
							  eps_f0_T = stateNew(i,12)
							  eps_ff_T = stateNew(i,13)
							  d1Plus = (eps_ff_T * (eps_f_T - eps_f0_T))/
     *						       (eps_f_T * (eps_ff_T - eps_f0_T))
							  d1Plus = min(one, max(zero,d1Plus))
							  d1Plus = max(stateNew(i,16), d1Plus)
							  stateNew(i,16) = d1Plus
						  ENDIF
C						  
					  ELSEIF ( stateNew(i,1) .LT. zero ) THEN
						  D = d1Minus    ! For initiates Damage variable when enter compresion stage
						  stateNew(i,9) = D
						  IF ( stateNew(i,22) .EQ. zero ) THEN
C						  
						      FIc = abs(Stress(1)/Xc)
C
							  IF ( FIc .GT. one ) THEN
								  stateNew(i,22) = one
								  sig_f0_C = Stress(1)
								  stateNew(i,18) = sig_f0_C
								  eps_f0_C = stateNew(i,1)
								  stateNew(i,19) = eps_f0_C
								  eps_ff_C = -2.d0*Gfc/(Xc*charLength(i)) ! eps_ff_C is minus
								  stateNew(i,20) = eps_ff_C
							  ENDIF
C							  
						  ELSEIF ( stateNew(i,22) .EQ. one ) THEN
							  eps_f_C = stateNew(i,1)
							  eps_f0_C = stateNew(i,19)
							  eps_ff_C = stateNew(i,20)
							  d1Minus = (eps_ff_C * (eps_f_C - eps_f0_C))/
     *						       (eps_f_C * (eps_ff_C - eps_f0_C))
							  d1Minus = min(one, max(zero,d1Minus))
							  d1Minus = max(stateNew(i,17), d1Minus)
							  stateNew(i,17) = d1Minus
						  ENDIF
C
					  ENDIF
C
C CRACK CLOSURE UNDER LOAD REVERSAL
C ---------------------------------					  
					  IF ( Stress(1) > zero ) THEN
						  D = max(d1Plus, d1Minus)
						  D = min(one, max(zero, D))
						  D = max(stateNew(i,9), D)
						  stateNew(i,9) = D						  
					  ELSEIF ( Stress(1) <= zero ) Then
					      D = d1Minus
						  D = min(one, max(zero, D))						  
						  D = max(stateNew(i,9), D)
						  stateNew(i,9) = D						  
					  ENDIF
C 
C STRESS DEGRADATION
C ==================
                      IF ( Stress(1) >= 0 ) THEN
						  IF ( D .LT. 0.9999d0 ) THEN
							  Stress(1) = (1-D)*Stress(1)
						  ELSE
						      stateNew(i,10) = zero
						  ENDIF
					  ELSEIF ( Stress(1) < 0) THEN
						  Stress(1) = (1-D)*Stress(1)
						  IF ( abs(stateNew(i,1)) > abs(stateOld(i,1)) ) THEN ! If monotonic compression, to avoid wrong degradation during compression reversal to zero. 
							  IF ( ABS(Stress(1)) < S_res ) THEN
								  Stress(1) = -S_res
							  ENDIF
						  ENDIF
					  ENDIF
C
C DEGRADE STRESSES IN OTHER DIRECTIONS
C ------------------------------------
					  ! IF ( D .LT. 0.9999d0 ) THEN
						  ! Stress(1) = (one-D)*Stress(1)						  
						  ! sgn2 = Stress(2)/abs(Stress(2)+safety)
						  ! IF ( sgn2. GT. zero ) THEN
							  ! Stress(2) = (one-D)*Stress(2)
						  ! ELSE
							  ! Stress(2) = Stress(2)
						  ! ENDIF
						  ! sgn3 = Stress(3)/abs(Stress(3)+safety)
						  ! IF ( sgn3. GT. zero ) THEN
							  ! Stress(3) = (one-D)*Stress(3)
						  ! ELSE
							  ! Stress(3) = Stress(3)
						  ! ENDIF
						  ! Stress(4) = (one-D)*Stress(4)
! C					  
					  ! ELSE
						  ! stateNew(i,10) = zero
					  ! ENDIF
C -------------------------------------					  
				  ENDIF ! FOR IF OF FAILURE FALG
C
C ELEMENT DELETION BASED ON DEFORMATION
C -------------------------------------				  
				  detF = det_3x3(defGradNew(i,:))
				  stateNew(i,6) = detF
				  IF (((detF > 0.0d0).and.(detF < 0.8d0)).or.(detF >= 1.6d0)) THEN
				  	  stateNew(i,10) = zero
				  ENDIF
C
C TRANSFER STRESS VECTOR BACK TO ABAQUS
C -------------------------------------
				  stressNew(i,1:4) = Stress(1:4)
C			  
			  ENDIF ! FOR IF OF ELEMENT DELETION FLAG
C				  
	      ENDDO
C	  
	  ENDIF
C	
	  RETURN 
C	  
	  END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                  C
C   FUNCTION det_3x3(A)                                            C
C   AIM: TO CALCULATE THE DETERMINANT OF AXISYMMETRIC ELEMENTS     C
C                                                                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function det_3x3(A) result(det)
        ! Computes the determinant of a 3x3 matrix
        ! matrix should be input as a 9 element vector
        implicit none
        ! input variables
        real*8, dimension(5), intent(in)  :: A ! 3x3 input matrix 
        real*8 :: det  ! determinant of A
        det = A(1)*A(2)*A(3) - A(4)*A(5)*A(3)
      end function det_3x3	
!-------------------------------------------------------------------!
!      Subroutine vumatcoh:                                         !
!-------------------------------------------------------------------!       

#include 'Abaqus_Definitions.f'
#include 'Fatigue_Globals.f'
#include 'Mesh_Utils.f'
#include 'Fatigue_ANN.f'
#include 'Fatigue.f'
#include 'CZM.f'

      !> VUMAT Main entry point from Abaqus
      subroutine vumatcoh(
     *     jblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relSpinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
     *     stressNew, stateNew, enerInternNew, enerInelasNew )

        use CZM, only: vumat_cohesive_fatigue
        include 'vaba_param.inc'

        dimension jblock(*), props(nprops),density(*), coordMp(*),
     *     charLength(*), strainInc(*), relSpinInc(*), tempOld(*),
     *     stretchOld(*), defgradOld(*), fieldOld(*), stressOld(*),
     *     stateOld(*), enerInternOld(*),enerInelasOld(*),
     *     tempNew(*), stretchNew(*), defgradNew(*), fieldNew(*),
     *     stressNew(*), stateNew(*), enerInternNew(*), enerInelasNew(*)

        character*80 cmname

        call vumat_cohesive_fatigue ( jblock(1),
     *     ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relSpinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
     *     stressNew, stateNew, enerInternNew, enerInelasNew,
     *     jblock(2), jblock(3), jblock(4), jblock(5))

      end subroutine vumatcoh

      
      !> Vexternaldb entry from Abaqus
      subroutine vexternaldb(lOp, i_Array, niArray, r_Array, nrArray)

        use Abaqus_Definitions, only: j_int_StartAnalysis
        use CZM, only: CZM_initialisation
        include 'vaba_param.inc'

        dimension i_Array(niArray), r_Array(nrArray)

        ! Initialisation at start of the analysis
        if (lOp .eq. j_int_StartAnalysis) then
          call CZM_initialisation()
        endif

      end subroutine vexternaldb			
