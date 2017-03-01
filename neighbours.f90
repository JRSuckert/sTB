MODULE sTB
  implicit none
  real*8, allocatable :: nn(:,:), k(:,:), t(:,:,:), oe(:), rOne(:,:)
  double complex, allocatable :: cOne(:,:), H0(:,:), G(:,:), t_k(:,:)
  integer            :: nn_size, nn_stages, layers, k_size, spin, H_size
  real*8 :: a1(2), a2(2)

CONTAINS

  SUBROUTINE initialize(l_nn_stages, l_layers, l_a1, l_a2, l_spin)
    INTEGER, INTENT(IN) :: l_nn_stages, l_layers, l_spin
    REAL*8, INTENT(IN) :: l_a1(2), l_a2(2)

    INTEGER :: i,j

    nn_stages = l_nn_stages
    layers = l_layers
    a1 = l_a1
    a2 = l_a2

    spin = l_spin

    H_size = spin * layers

    ALLOCATE( rOne(H_size, H_size), cOne(H_size,H_size) )
    rOne = 0.D0
    cOne = (0.D0, 0.D0)
    DO i=1, H_size
      rOne(i,i) = 1.D0
      cOne(i,i) = (1.D0, 0.D0)
    END DO

    ALLOCATE( H0(H_size, H_size), G(H_size, H_size), t_k(layers,layers) )
    H0=0.0D0
    G=0.0D0
    t_k=0.0D0

    CALL init_nn()
    CALL init_t()
    CALL init_oe(0.0D0)
    CALL init_k()

  END SUBROUTINE

  SUBROUTINE init_t()
    IMPLICIT NONE
    INTEGER                          :: i,p,q        ! loop variables

    ALLOCATE(t(nn_stages+1, layers, layers))
    DO p = 1, layers
      DO q = 1, layers
        DO i = 1, nn_size+1
          t(i,p,q) = 1.0D0 / 10.0D0**(i-2 + ABS(p-q))
        END DO
        IF(p == q) THEN
          t(1,p,q) = 0.D0
        END IF
      END DO
    END DO
    print *, t
    RETURN
  END SUBROUTINE

  SUBROUTINE init_oe(E)
    implicit none
    real*8, intent(in) :: E

    allocate( oe(layers) )
    oe = E
    return
  END SUBROUTINE

  SUBROUTINE init_nn()
    IMPLICIT NONE
    INTEGER                          ::  k, i, x, y

    nn_size = nn_stages*(nn_stages+1)*2 + 1
    ALLOCATE( nn(2,nn_size) )

    i = 1
    nn(1,i) = 0.0D0
    nn(2,i) = 0.0D0
    i = i+1
    DO k = 1, nn_stages
      DO x = -k, k
        y = k - ABS(x)
        nn(1,i) = x * a1(1) + y * a2(1)
        nn(2,i) = x * a1(2) + y * a2(2)
        i = i+1
        IF(y /= 0) THEN
          y = -1 * y
          nn(1,i) = x * a1(1) + y * a2(1)
          nn(2,i) = x * a1(2) + y * a2(2)
          i = i+1
        END IF
      END DO
    END DO
    print *, nn(1,:)
    print *, nn(2,:)
    RETURN
  END SUBROUTINE

  SUBROUTINE init_k()
    USE global, global_a1 => a1, global_a2 => a2
    INTEGER :: i

    global_a1 = a1
    global_a2 = a2

    CALL generate_kpoints()
    k_size = eff_nptk

    !k_size = 1

    ALLOCATE( k(2,k_size) )
    DO i=1, k_size
     k(:,i) = kbz(i,:)
    END DO


      !k(1,1) = 0.D0
      !k(2,1) = 0.D0
    RETURN
  END SUBROUTINE


  SUBROUTINE k_hopping_sum(kp)
    IMPLICIT NONE
    REAL*8,               INTENT(IN)   :: kp(2)

    INTEGER :: i, j, l, p, q
    DOUBLE COMPLEX :: iOne

    iOne = CMPLX(0.D0,1.D0)

    t_k = 0.D0
    DO p = 1, layers
      t_k(p,p) = t_k(p,p) + oe(p)
      DO q = 1, layers
        IF(ABS(p-q) <= nn_stages) THEN
          l = 1                             ! Running Index for next neighbours nn
          t_k(p,q) = t_k(p,q) + t(1,p,q)    ! R = (0,0) term
          l = l+1
          DO i = 1, nn_stages - ABS(p-q)            ! sum over next neighbour stages
            DO j = 1, 4*i                   ! sum over sites with same distance to (0,0)
              t_k(p, q) = t_k(p,q) + EXP(-iOne * DOT_PRODUCT(kp, nn(:,l))) * t(i+1, p, q)
              !print *, "exp", -iOne * DOT_PRODUCT(kp, nn(:,l))
              l = l+1
            END DO
          END DO
          !print *, t_k(p,q)

        END IF
      END DO
    END DO
    RETURN
  END SUBROUTINE

  SUBROUTINE hamilt_0(kp)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: kp(2)
    INTEGER :: i, s

    H0 = 0.0D0

    CALL k_hopping_sum(kp)
    H0(1:layers, 1:layers) = H0(1:layers, 1:layers) + t_k
    if(spin == 2) then
        H0(layers+1:2*layers, layers+1:2*layers) = H0(layers+1:2*layers, layers+1:2*layers) + t_k
    end if

    !H0 = H0 / k_size
    !print *, H0
    ! DO i=1, layers
    !   WRITE(*, *) H(:, i)
    ! END DO

  END SUBROUTINE

  SUBROUTINE green(kp, z)
    DOUBLE COMPLEX, INTENT(IN) :: z
    REAL*8, INTENT(IN) :: kp(2)
    INTEGER :: i, info
    INTEGER :: ipiv(H_size), lwork, work(H_size)

    lwork = H_size
    G = CMPLX(0.0D0, 0.D0)
    CALL hamilt_0(kp)

    DO i = 1,H_size
      G(i,i) = CMPLX(1.0D0, 0.D0)
      H0(i,i) = H0(i,i) - z
    END DO

    H0 = -H0

    CALL ZGETRF(     H_size, H_size, H0, H_size, ipiv,            info)
    CALL ZGETRS('N', H_size, H_size, H0, H_size, ipiv, G, H_size, info)

    !print *, G
  END SUBROUTINE

  SUBROUTINE ldos(z, EF)
    DOUBLE COMPLEX, INTENT(IN) :: z
    DOUBLE COMPLEX, ALLOCATABLE :: G_sum(:,:)
    DOUBLE COMPLEX :: E
    REAL*8 :: h, e1, e2, EF
    REAL*8, allocatable :: occ(:,:)
    INTEGER :: n, i, j

    ALLOCATE(G_sum(H_size,H_size))

    e1 = -EF
    e2 =  EF
    n = 1E+3
    h = (e2-e1) / n


    ALLOCATE(occ(H_size, n+1))

    open(unit=20, file='ldos.txt', status='replace', action='write')


    DO i = 0, n
      E = CMPLX( e1+i*h, IMAG(z) )
      G_sum = 0.D0
      DO j=1, k_size
        CALL green(k(:,j), E)
        G_sum = G_sum + G
      END DO
      G_sum = G_sum / k_size
      DO j=1, H_size
        occ(j,i+1) = -IMAG(G_sum(j,j)) / 3.1415926D0
      END DO
      write (20,*) e1+i*h, occ(1,i+1), occ(layers+1, i+1)
    END DO

    close(unit=20)


  END SUBROUTINE

END MODULE sTB

PROGRAM neighbours
  USE sTB
  IMPLICIT NONE
  INTEGER ::  n, l
  REAL*8  ::  vec1(2), vec2(2), EF
  DOUBLE COMPLEX :: z = 0.0D0

  n = 1
  l = 1

  vec1(1) = 1.0D0
  vec1(2) = 0.0D0

  vec2(1) = 0.0D0
  vec2(2) = 1.0D0

  CALL initialize(n, l, vec1, vec2, 2)
  z = (0.D0, 1.D-5)
  EF = 15.D0
  CALL ldos(z, EF)


END PROGRAM neighbours
