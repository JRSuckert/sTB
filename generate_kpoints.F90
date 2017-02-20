MODULE global
   IMPLICIT NONE
   !Setup
   INTEGER, PARAMETER :: nptk=2000 !The num. k-points desired (roughly)
   !Generate points only inside 1st BZ? If .false. generate in the 
   ! paralelogram determined by b1 and b2
   LOGICAL, PARAMETER :: firstBZ=.TRUE.
   !Lattice definition
   REAL(KIND=8), PARAMETER :: a1(2)=[0.5d0, SQRT(3.d0)/2.d0]
   REAL(KIND=8), PARAMETER :: a2(2)=[0.5d0,-SQRT(3.d0)/2.d0]
   !Globals, constants and non-ajustables
   REAL(KIND=8), PARAMETER :: pi=2.d0*ASIN(1.d0)
   !eff_nptk = effective num. k-points. nptk is only a start point.
   INTEGER :: eff_nptk
   !To contain the kpoints and its corresponding weight.
   REAL(KIND=8), ALLOCATABLE :: kbz(:,:), wkbz(:)
END MODULE global

SUBROUTINE generate_kpoints()
   !This sub generates kpoints in the BZ. If firstBZ is .TRUE. it will
   ! be generated inside of the 1st BZ, otherwise in the parelogram 
   ! determined by b1 and b2.
   !Author: Flaviano José dos Santos, Dec. 2015
   !Institution: Forschungszentrum Jülich 
   
   !INTEND:          in  in  in    out       out  out   in  in
   USE global, ONLY: a1, a2, nptk, eff_nptk, kbz, wkbz, pi, firstBZ
   IMPLICIT NONE
   REAL(KIND=8) :: BZ(4,2), smallest_dist, distance, diff(2), &
                   b1(2), b2(2), ini_smallest_dist
   REAL(KIND=8), ALLOCATABLE :: inikbz(:,:),iniwkbz(:)
   REAL(KIND=8) :: extrakbz(nptk*10,2), extrawkbz(nptk*10)
   INTEGER :: l,j,m, smallest_indice, numextrakbz
   INTEGER :: nptk_perdim !n. of k point per dimension

   !Determining the amount of point per dimension. It imposes to have
   ! a number of kpoints per dimesion divisible by 6.
   nptk_perdim=CEILING(SQRT(DBLE(nptk))/6.d0)
   nptk_perdim=nptk_perdim*6
   eff_nptk=nptk_perdim**2
   ALLOCATE( iniwkbz(eff_nptk), inikbz(eff_nptk,2) )
   
   !Reciprocal lattice determination
   b1(:)=        [a2(2),-a2(1)]
   b1(:)=2.d0*pi*[a2(2),-a2(1)]/DOT_PRODUCT(a1,b1)
   b2(:)=        [a1(2),-a1(1)]
   b2(:)=2.d0*pi*[a1(2),-a1(1)]/DOT_PRODUCT(a2,b2)

   !These vectors indicate the BZ's in which the k-point will be
   ! genarated. This will be used for redistribute them in the 1st BZ
   BZ(1,:)=0.d0
   BZ(2,:)=b1
   BZ(3,:)=b2
   BZ(4,:)=b1+b2

   !Generate k-points in the paralelogram determined by b1 and b2
   iniwkbz=1.d0
   m=0
   DO l=1, nptk_perdim; DO j=1, nptk_perdim !Run over the entire BZ
      m=m+1
      inikbz(m,:)= ( DBLE(l-1)*b1 + DBLE(j-1)*b2 )/ DBLE(nptk_perdim)
   END DO; END DO

   !Translate the k-points to the 1st BZ.
   IF(firstBZ) THEN
      !10*|b1+b2|, bigger than the distance of any genarated kpoint
      ini_smallest_dist=10.d0*SQRT(DOT_PRODUCT(b1+b2,b1+b2))
      numextrakbz=0
      !Run over all the kpoints generated initially.
      DO l=1, eff_nptk
         smallest_dist=ini_smallest_dist
         m=0
         !Checks to each of the 4 BZ's the kpoint belongs by checking
         ! to each BZ it's closer.
         DO j=1, 4
            diff=inikbz(l,:)-BZ(j,:)
            distance=SQRT(DOT_PRODUCT(diff,diff))
            IF(distance.LT.smallest_dist) THEN
               smallest_dist=distance
               smallest_indice=j
            END IF
         END DO
         !Checks if the kpoint is in the border between two or more
         ! BZ's. If yes, create a clone of it to translate later into
         ! the 1st BZ. 
         DO j=1, 4
            diff=inikbz(l,:)-BZ(j,:)
            distance=SQRT(DOT_PRODUCT(diff,diff))
            IF( ( ABS(distance-smallest_dist) .LT. 1.d-12 ) .AND. &
                                          j.NE.smallest_indice ) THEN
               m=m+1
               numextrakbz=numextrakbz+1
               extrakbz(numextrakbz,:)=inikbz(l,:)-BZ(j,:)
            END IF
         END DO
         IF(m.NE.0) THEN
            !The weight of the kpoint in the border is shared with
            ! its clones.
            iniwkbz(l)=1.d0/DBLE(m+1)
            extrawkbz(numextrakbz-m+1:numextrakbz)=1.d0/DBLE(m+1)
         END IF
         !Translate the kpoint to the 1st BZ
         inikbz(l,:)=inikbz(l,:)-BZ(smallest_indice,:)
      END DO
   END IF !IF(firstBZ)

   ALLOCATE( wkbz(eff_nptk+numextrakbz),kbz(eff_nptk+numextrakbz,2) )
   !The final array of kpoints will be the initial one by the clones
   ! of the ones in the border between BZ's.
   kbz (1:eff_nptk,:)=inikbz (:,:)
   wkbz(1:eff_nptk)  =iniwkbz(:)
   IF(numextrakbz.NE.0) THEN
      kbz(eff_nptk+1:eff_nptk+numextrakbz,:)=extrakbz(1:numextrakbz,:)
      wkbz(eff_nptk+1:eff_nptk+numextrakbz) =extrawkbz(1:numextrakbz)
   END IF
   wkbz=wkbz/DBLE(eff_nptk)
   eff_nptk=eff_nptk+numextrakbz
   
   PRINT *, "firstBZ = ", firstBZ
   PRINT "(a,i0,a,i0)", "nptk=",eff_nptk,"  nptk_perdim=",nptk_perdim
END SUBROUTINE generate_kpoints