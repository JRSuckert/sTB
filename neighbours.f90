MODULE sTB
CONTAINS

  SUBROUTINE initialize_t(n, layers, t, oE)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: oE                           !on-site energy
    REAL*8, ALLOCATABLE, INTENT(OUT) :: t(:,:,:)       !output hopping matrix
    INTEGER, INTENT(IN) :: n, layers                 ! # next neighbour / layers

    INTEGER :: i,p,q                                 ! loop variables

    ALLOCATE(t(n+1, layers, layers))



    DO p = 1, layers
      DO q = 1, layers
        DO i = 1, n+1
          t(i,p,q) = 1.0 / 10.0**(i-2 + ABS(p-q))
        END DO
        IF(p == q) THEN
          t(1,p,q) = oE
        END IF
      END DO
    END DO

    RETURN
  END SUBROUTINE

  SUBROUTINE initialize_k(Nk, k, vec1, vec2)
    USE global
    INTEGER,              INTENT(OUT) :: Nk
    REAL*8,    ALLOCATABLE, INTENT(OUT) :: k(:,:)
    REAL*8, INTENT(IN) :: vec1(2), vec2(2)
    INTEGER :: i
    a1 = vec1
    a2 = vec2
    CALL generate_kpoints()

    ! ------------>TODO<--------------- !
    Nk = eff_nptk
    ALLOCATE(k(2,Nk))
    DO i=1, Nk
      k(:,i) = kbz(i,:)
    END DO
    PRINT *, k(:,Nk)
    RETURN

  END SUBROUTINE

  SUBROUTINE next_neighbour(n, m, neighbour , vec1, vec2)
    IMPLICIT NONE

    INTEGER,              INTENT(IN)  ::  n    ! number of next neighbour stages
    INTEGER,              INTENT(OUT) ::  m   ! array size
    REAL*8,                 INTENT(IN)  ::  vec1(2), vec2(2)
    REAL*8,    ALLOCATABLE, INTENT(OUT) ::  neighbour(:,:)
    INTEGER                           ::  k, i, x, y

    m = n*(n+1)*2 + 1
    ALLOCATE(neighbour(2,m))

    i = 1
    neighbour(1,i) = 0.0
    neighbour(2,i) = 0.0
    i = i+1
    DO k = 1, n
      DO x = -k, k
        y = k - ABS(x)
        neighbour(1,i) = x * vec1(1) + y * vec2(1)
        neighbour(2,i) = x * vec1(2) + y * vec2(2)
        i = i+1
        IF(y /= 0) THEN
          y = -1 * y
          neighbour(1,i) = x * vec1(1) + y * vec2(1)
          neighbour(2,i) = x * vec1(2) + y * vec2(2)
          i = i+1

        END IF
      END DO
    END DO

    RETURN
  END SUBROUTINE

  SUBROUTINE k_hopping_sum(n, m, layers, k, neighbour, t, t_k)
    IMPLICIT NONE
    INTEGER,              INTENT(IN)   :: n,m, layers
    REAL*8,               INTENT(IN)   :: neighbour(2,m)
    REAL*8,               INTENT(IN)   :: t(n+1, layers, layers) !hopping parameter, +1 for e.g. t_01(0)
    REAL*8,               INTENT(IN)   :: k(2)
    DOUBLE COMPLEX, ALLOCATABLE, INTENT(OUT)  :: t_k(:,:)

    INTEGER :: i, j, l, p, q
    DOUBLE COMPLEX :: iOne
    ALLOCATE(t_k(layers,layers))

    iOne = (0,1)

    t_k(:,:) = 0.0
    DO p = 1, layers
      DO q = 1, layers
        IF(ABS(p-q) <= n) THEN
          l = 1                             ! Running Index for next neighbours nn
          t_k(p,q) = t_k(p,q) + t(1,p,q)    ! R = (0,0) term
          l = l+1
          DO i = 1, n - ABS(p-q)            ! sum over next neighbour stages
            DO j = 1, 4*i                   ! sum over sites with same distance to (0,0)
              t_k(p, q) = t_k(p,q) + EXP(-iOne * DOT_PRODUCT(k, neighbour(:,l))) * t(i+1, p, q)
              l = l+1
            END DO
          END DO
        END IF
      END DO
    END DO
    RETURN
  END SUBROUTINE

  SUBROUTINE build_H0(n, layers, vec1, vec2, H)
    IMPLICIT NONE
    INTEGER,                     INTENT(IN)  :: n, layers
    REAL*8,                      INTENT(IN)  :: vec1(2), vec2(2)
    DOUBLE COMPLEX, ALLOCATABLE, INTENT(OUT) :: H(:,:)
    REAL*8 :: Energy = 3.3
    REAL*8,    ALLOCATABLE :: t(:,:,:), nn(:,:), k(:,:)
    DOUBLE COMPLEX, ALLOCATABLE :: t_k(:,:)

    INTEGER :: m, i, Nk

    ALLOCATE(H(layers,layers))
    H(:,:) = 0

    CALL initialize_t(n,layers, t, Energy)
    CALL initialize_k(Nk, k, vec1, vec2)

    CALL next_neighbour(n, m, nn, vec1, vec2)

    PRINT *, Nk
    DO i = 1, Nk
      CALL k_hopping_sum(n,m, layers, k(:,i), nn, t, t_k)
      H = H + t_k
    END DO
    H = H / layers

    DO i=1, layers
      WRITE(*, *) H(:, i)
    END DO

  END SUBROUTINE

  SUBROUTINE calc_green(n, sites, H, z, G)
    DOUBLE COMPLEX, INTENT(IN) :: H(n,n), z
    DOUBLE COMPLEX, ALLOCATABLE, INTENT(OUT) :: G(:,:)
    DOUBLE COMPLEX, ALLOCATABLE :: Hz(:,:)
    INTEGER, INTENT(IN) :: sites
    INTEGER :: i,j, info
    INTEGER :: ipiv(n), lwork, work(40*n)

    ALLOCATE(G(n,n), Hz(n,n))
    lwork = 40*n
    ! G = 0.0
    ! DO i = 1,n
    !   G(i,i) = 1.0
    !   Hz(i,i) = z - H(i,i)
    ! END DO

    !CALL CPBTRF('U', n, sites, Hz, n, info)
    !CALL CPBTRS('U', n, sites, n, Hz, n, G, n, info)
    G = H
    DO i=1,n
      G(i,i) = G(i,i) - z
    END DO

    CALL ZGETRF(n, n, G, n, ipiv, info)

    CALL ZGETRI(n, G, n, ipiv, work, lwork, info)
    PRINT *, " "
    DO i=1, n
        PRINT *,  G(:, i)
    END DO

  END SUBROUTINE

END MODULE sTB

PROGRAM neighbours
  USE sTB
  IMPLICIT NONE
  INTEGER ::  n, layers, i
  REAL*8  ::  vec1(2), vec2(2)
  DOUBLE COMPLEX, ALLOCATABLE :: H0(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: G(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: Ident(:,:)
  DOUBLE COMPLEX :: E = 3.0
  DOUBLE COMPLEX :: One = 1.0
! n = 1 -> 4 neighbours (4)
! n = 2 -> 12 neighbours (8)
! n = 3 -> 24 neighbours (12)
! n = 4 -> 40 neighbours (16)
! (3*4)/2 * 4 =
! # of neighbours = n*(n+1)*2

! find all x,y s.th. |x| + |y| = k
!

  n = 1
  layers = 4

  vec1(1) = 1
  vec1(2) = 0

  vec2(1) = 0
  vec2(2) = 1

  ALLOCATE(Ident(layers,layers))

  CALL build_H0(n, layers, vec1, vec2, H0)
  CALL calc_green(layers, n, H0, E, G)

  CALL CGEMM('N', 'N', layers, layers, layers, One, H0, layers, G, layers, One, Ident, layers)

  PRINT *, " "
  DO i=1, layers
    WRITE(*, *) Ident(:, i)
  END DO


END PROGRAM neighbours
