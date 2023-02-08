      subroutine vumat(
     *     jblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relSpinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
     *     stressNew, stateNew, enerInternNew, enerInelasNew )  
      
      include 'vaba_param.inc'
      
      dimension jblock(*), props(nprops),density(*), coordMp(*),
     1     charLength(*), strainInc(*),
     2     relSpinInc(*), tempOld(*),
     3     stretchOld(*),
     4     defgradOld(*),
     5     fieldOld(*), stressOld(*),
     6     stateOld(*), enerInternOld(*),
     7     enerInelasOld(*), tempNew(*),
     8     stretchNew(*),
     9     defgradNew(*),
     1     fieldNew(*),
     2     stressNew(*), stateNew(*),
     3     enerInternNew(*), enerInelasNew(*)

      character*80 cmname
       
      integer readflag
      real*4  edgedim(350000,3)
      common /integerbuf/readflag 
      common / realbuf/edgedim

        
        if(cmname(1:6).eq.'MATCOH')then
        
          
          call vumatcoh(jblock(1),
     *     ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relSpinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
     *     stressNew, stateNew, enerInternNew, enerInelasNew,
     *     jblock(5), jblock(2),
     *     jblock(3), jblock(4))
           
        elseif(cmname(1:8).eq.'MATFIBRE')then

           call  vumatfibre( jblock(1),
     *     ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relSpinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
     *     stressNew, stateNew, enerInternNew, enerInelasNew,
     *     jblock(5), jblock(2),
     *     jblock(3), jblock(4))
     
            endif
 
        
      return
      end       

!--------------------------------------------------------------------!
!                                                                    !
!  subroutine: vumatfibre                                            !
!  function: Fibre failure activation using maximum stress criteria  !
!            and evolution                                           ! 
!                                                                    !
!--------------------------------------------------------------------!  
! SDV1  - strain 11
! SDV2  - strain 22
! SDV3  - strain 33
! SDV4  - strain 12, engineering shear strain gamma
! SDV5  - Temporary strain 12
! SDV6  - Shear softening flag: 0- elastic, 1- softening
! SDV7  - Shear failure indice
! SDV8  - Shear strain at shear damage initiation
! SDV9  - Inelastic shear strain at shear damage initiation
! SDV10 - Shear strain at final shear failure
! SDV11 - Shear damage variable, Ds
! SDV12 - Element characteristic length


C..........fibre_damage starts here..............

		 subroutine vumatfibre(
     *     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relSpinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
     *     stressNew, stateNew, enerInternNew, enerInelasNew,
     *     nElement, nMatPoint, nLayer, nSecPoint )
     
      include 'vaba_param.inc'
      
C
      dimension coordMp(nblock,*), charLength(nblock), props(nprops),
     1     density(nblock), strainInc(nblock,ndir+nshr),
     2     relSpinInc(nblock,nshr), tempOld(nblock),
     3     stretchOld(nblock,ndir+nshr), 
     4     defgradOld(nblock,ndir+nshr+nshr),
     5     fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6     stateOld(nblock,nstatev), enerInternOld(nblock),
     7     enerInelasOld(nblock), tempNew(nblock),
     8     stretchNew(nblock,ndir+nshr),
     9     defgradNew(nblock,ndir+nshr+nshr),
     1     fieldNew(nblock,nfieldv),
     2     stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     3     enerInternNew(nblock), enerInelasNew(nblock),
     4     nElement(nblock)

       parameter(zero = 0.d0, one = 1.d0, two = 2.d0, three = 3.d0, 
     * third = one/three, half = 0.5d0, twothird = two/three,
     * threehalf = 1.5d0, safety = 1.d-20, tolSgn = 1.d-20, S_res=50.d0)              
      
       character*80 cmname,cpname
       character*256 outdir,fullpath
       integer intnum,locnum,jrcd,lenoutdir
       double precision E11,E22,E33,v12,v21,v13,v31,v23,v32,G12,G13,G23
       double precision C11 ,C12 ,C13 ,C14 ,C21 ,C22,
     *   C23 , C24, C31 ,C32 ,C33 ,C34,
     *   C41 , C42 ,C43 ,C44
	   double precision n,m,delta
	   double precision Xt,Xc,Yc,Yt,Sl
	   double precision Gft,Gfc
       double precision sgn4,A,B,TAU12,FIs,FIs_old,FIs_max,d12
	   double precision gamma12_f0,gamma12_inel_f0,gamma12_ff, Ds
	   
	   
	   integer i,j



      E11 = props(1) 
      E22 = props(2) 
      E33 = props(3) 
      v12 = props(4) 
      v13 = props(5) 
      v23 = props(6) 
      G12 = props(7) 
      G13 = props(8) 
      G23 = props(9) 

	! Parameters for nonlinear shear constitutive law
	
      A   = props(10) !145 (paper)
      B   = props(11) !38
	 
	! Ply strengths	
	  Xt  = props(12)  ! 1089.2/1906
	  Xc  = props(13)  ! 675.6/1182
      Yt  = props(14)  ! 11.1/44.4
      Yc  = props(15)  ! 25/100
      Sl  = props(16)  ! 45.2
	  
      Gft = props(17)  ! Fracture energy in lonitudinal tension = 92 kJ/m2
	  Gfc = props(18)  ! Fracture energy in lonitudinal compression = 80

      v21 = v12*E22/E11
      v31 = v13*E33/E11
      v32 = v12

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
      

	  IF (stepTime .eq. zero) THEN

		  Do i = 1, nblock

            stressNew(i,1)= stressOld(i,1) + C11*strainInc(i,1)
     *				+ C12*strainInc(i,2) + C13*strainInc(i,3)
		    stressNew(i,2)= stressOld(i,2) + C21*strainInc(i,1)
     *				+ C22*strainInc(i,2) + C23*strainInc(i,3)
            stressNew(i,3)= stressOld(i,3) + C31*strainInc(i,1)
     *				+ C32*strainInc(i,2) + C33*strainInc(i,3)
            stressNew(i,4) = stressOld(i,4) + two*C44*strainInc(i,4)
			
			Ds = stateNew(i,11)
			stateNew(i,6) = zero
			stateNew(i,12) = charLength(i)
				  
		  ENDDO

	  ELSE

		  DO i = 1, nblock
			  
			  ! DO j = 1, nblock
			      ! stateNew(i,j) = stateOld(i,j)
			  ! ENDDO

			  ! WRITE(*,*) '----------- stepTime:', stepTime, '----------------'
			  
			  stateNew(i,1)= stateOld(i,1)+ strainInc(i,1)
			  stateNew(i,2)= stateOld(i,2)+ strainInc(i,2)
			  stateNew(i,3)= stateOld(i,3)+ strainInc(i,3)
			  stateNew(i,4)= stateOld(i,4)+ two*strainInc(i,4)

			  stressNew(i,1)= C11*stateNew(i,1)+C12*stateNew(i,2)
     *                  +C13*stateNew(i,3)

			  stressNew(i,2)= C21*stateNew(i,1)+C22*stateNew(i,2)
     *                 + C23*stateNew(i,3)

			  stressNew(i,3)= C31*stateNew(i,1)+C32*stateNew(i,2)
     *                 + C33*stateNew(i,3) 		
			  
			  sgn4 = stateNew(i,4)/(abs(stateNew(i,4)+safety))

			  IF(abs(stateNew(i,4)).gt.stateOld(i,5))THEN
				stateNew(i,5)=abs(stateNew(i,4))
				stressNew(i,4)=sgn4*(A*(one -exp(-B*stateNew(i,5))))  
			  ELSE
				stateNew(i,5)=stateOld(i,5)
				stressNew(i,4)=sgn4*(A*(one -exp(-B*stateNew(i,5))) 
     *     				-G12*(stateNew(i,5)-abs(stateNew(i,4))))
				
			  ENDIF
			  
			  TAU12 = stressNew(i,4)
C Calculate response from elastic region to failure initiation point
C ==================================
			  IF (stateNew(i,6) .eq. zero) THEN !SDV6: shear softening flag
				  FIs = abs(TAU12)/S12
				  FIs_old = stateOld(i,7)
					
				  FIs_max = max(FIs,FIs_old)
				  stateNew(i,7) = FIs_max 
				  
				  d12 = zero
				  
				  IF (FIs .GT. one) THEN
				  
					  stateNew(i,6) = one
					  
					! Recod the gamma12 at shear damage initiation
					  gamma12_f0 = stateNew(i,4)
					  stateNew(i,8) = gamma12_f0
					
					! Calculate the inelastic strain at shear damage initiation
					  gamma12_inel_f0 = gamma12_f0 - S12/G12
					  stateNew(i,9) = gamma12_inel_f0
					
					! Calculate the gamma12 at final shear failure
					  gamma12_ff = 2.d0*Gs/(S12*charLength(i))
					  stateNew(i,10) = gamma12_ff
					  
				  ENDIF
			  ELSEIF (stateNew(i,6) .eq. one) THEN
				  gamma12_f0 = stateNew(i,8)
				  gamma12_inel_f0 = stateNew(i,9)
				  gamma12_ff = stateNew(i,10)
				  gamma12 = stateNew(i,4)
				  
				  d12 = ((gamma12_ff - gamma12_inel_f0)/(gamma12_ff - gamma12_f0))
     *             *(one - (gamma12_f0 - gamma12_inel_f0)
     *			   /(gamma12 - gamma12_inel_f0))
	
				  Ds = min(one, max(zero, d12))
				  Ds = max(stateNew(i,11),d12)
				  stateNew(i,11) = Ds
				  
				  TAU12 = TAU12*(1-Ds)*(gamma12 - gamma12_inel_f0)
			  ENDIF
			     
				  
				
		   
			  ! IF(abs(stateNew(i,4)).gt.abs(stateOld(i,5)))THEN
				
				  ! stateNew(i,5)=abs(stateNew(i,4)) ! SDV5: experienced max shear strain
				  ! stressNew(i,4)=sgn4*(A*(one -exp(-B*stateNew(i,5))))
				  ! s12_max = stressNew(i,4)
				
			  ! ELSE
			  
				  ! stateNew(i,5)=abs(stateOld(i,5))
				  ! stressNew(i,4)=sgn4*(s12_max
 				! - A*(one -exp(-B*(abs(stateNew(i,5)-abs(stateNew(i,4)))))))
				
			  ! ENDIF
		  ENDDO
	  ENDIF

      RETURN
      END
	  
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
            
