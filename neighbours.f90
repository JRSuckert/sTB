MODULE sTB
CONTAINS

  SUBROUTINE intialize_t(n, layers, t)
    IMPLICIT NONE

    REAL, ALLOCATABLE, INTENT(OUT) :: t(:,:,:)
    INTEGER, INTENT(IN) :: n, layers

    INTEGER :: i,p,q

    ALLOCATE(t(n+1, layers, layers))



    DO p = 1, layers
      DO q = 1, layers
        DO i = 1, n+1
          t(i,p,q) = 1.0 / 10.0**(i-2 + ABS(p-q))
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE

  SUBROUTINE next_neighbour(n, m, neighbour , vec1, vec2)
    IMPLICIT NONE

    INTEGER ::  k, i, x, y
    INTEGER, INTENT(IN) :: n    ! number of next neighbour stages
    INTEGER, INTENT(OUT) :: m   ! array size
    REAL, INTENT(IN)  ::  vec1(2), vec2(2)
    REAL, ALLOCATABLE, INTENT(OUT) ::  neighbour(:,:)

    m = n*(n+1)*2
    ALLOCATE(neighbour(2,m))

    i = 1
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

    ! DO i = 1, m
    !   PRINT *, neighbour(1,i), neighbour(2,i)
    ! END DO
    RETURN
  END SUBROUTINE

  SUBROUTINE k_hopping_sum(n, m, layers, k, neighbour, t, t_k)
    INTEGER, INTENT(IN)   :: n,m, layers
    REAL,    INTENT(IN)   :: neighbour(2,m)
    REAL,    INTENT(IN)   :: t(n+1, layers, layers) !hopping parameter, +1 for e.g. t_01(0)
    REAL,    INTENT(IN)   :: k(2)
    COMPLEX, ALLOCATABLE, INTENT(OUT)  :: t_k(:,:)

    INTEGER :: i, j, l, p, q
    COMPLEX :: iOne
    ALLOCATE(t_k(layers,layers))

    iOne = (0,1)

    DO p = 1, layers
      DO q = 1, layers
        l = 1
        t_k(p,q) = 0

        IF(p == q) THEN
          DO i = 1, n
            DO j = 1, 4*i
              t_k(p, q) = t_k(p,q) + EXP(-iOne * DOT_PRODUCT(k, neighbour(:,l))) * t(i+1, p, q)
              ! PRINT *, i, j, "((", k, ")) ((", neighbour(:,l), "))", DOT_PRODUCT(k,neighbour(:,l)), t(i+1, p,q), t_k(p,q)
              l = l+1
            END DO
          END DO
        ELSE IF(ABS(p-q) <= n) THEN
          t_k(p,q) = t_k(p,q) + t(1,p,q)
          DO i = 1, n - ABS(p-q)
            DO j = 1, 4*i
              t_k(p, q) = t_k(p,q) + EXP(-iOne * DOT_PRODUCT(k, neighbour(:,l))) * t(i+1, p, q)
              l = l+1
            END DO
          END DO
        END IF
      END DO
    END DO
    RETURN
  END SUBROUTINE
END MODULE sTB

PROGRAM neighbours
  USE sTB
  IMPLICIT NONE
  INTEGER ::  n, m, layers, i,j
  REAL  ::  vec1(2), vec2(2), k(2)
  REAL, ALLOCATABLE ::  neighbour(:,:), t(:,:,:)
  COMPLEX, ALLOCATABLE :: t_k(:,:)

! n = 1 -> 4 neighbours (4)
! n = 2 -> 12 neighbours (8)
! n = 3 -> 24 neighbours (12)
! n = 4 -> 40 neighbours (16)
! (3*4)/2 * 4 =
! # of neighbours = n*(n+1)*2

! find all x,y s.th. |x| + |y| = k
!

  n = 2
  layers = 3

  vec1(1) = 1
  vec1(2) = 0

  vec2(1) = 0
  vec2(2) = 1

  CALL intialize_t(n, layers, t)

  CALL next_neighbour(n, m, neighbour, vec1,vec2)

  k(1) = SQRT(2.0)
  k(2) = SQRT(2.0)

  CALL k_hopping_sum(n,m, layers, k, neighbour, t, t_k)

  DO i = 1, layers
    DO j = 1, layers
      print *, i, j, t_k(i,j)
    END DO
  END DO
  DEALLOCATE(neighbour)

END PROGRAM neighbours
