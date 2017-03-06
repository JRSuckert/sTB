MODULE sTB
  implicit none
  real*8, allocatable :: nn(:,:), k(:,:), t(:,:,:), oe(:), rOne(:,:), wi(:), xi(:)
  double complex, allocatable :: cOne(:,:), H0(:,:), G(:,:), t_k(:,:)
  integer            :: nn_size, nn_stages, layers, k_size, spin, H_size, xi_size
  real*8 :: a1(3), a2(3), a3(3)
  real*8 :: pi = 3.1415926535897932d0

CONTAINS

  SUBROUTINE initialize(l_nn_stages, l_layers, l_a1, l_a2, l_a3, l_spin, p)
    INTEGER, INTENT(IN) :: l_nn_stages, l_layers, l_spin
    REAL*8, INTENT(IN) :: l_a1(3), l_a2(3), l_a3(3)
    real*8, intent(in) :: p(3)
    INTEGER :: i,j

    nn_stages = l_nn_stages
    layers = l_layers
    a1 = l_a1
    a2 = l_a2
    a3 = l_a3
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

    CALL next_neighbour_init(a1, a2, a3, p)
    CALL init_t()
    CALL init_oe(0.0D0)
    !CALL init_k()
    call generate_chebychev_points(128)
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

  subroutine next_neighbour_init(a1,a2,a3,p)
    implicit none
    real*8, intent(in) :: a1(3)
    real*8, intent(in) :: a2(3)
    real*8, intent(in) :: a3(3)
    real*8, intent(in) :: p(3)
    integer :: i
    integer :: j
    integer :: l
    integer :: m
    integer :: n
    integer :: nnt
    real*8, allocatable :: tmp(:,:,:)          !> @var Array contains all neighbours in the ellipsoid 2*nn_stages*(a1,a2,a3)
    real*8, allocatable :: dist(:)             !> @var Contains all distances in the ellipsoid 2*nn_stages*(a1,a2,a3)
    real*8, allocatable :: on_plane(:,:,:)     !> @var Array of type (3,2,on_cnt), all coordinates and directional cosines of neighbours in plane described by vector
    real*8, allocatable :: off_plane(:,:,:)    !> @var Array of type (3,2,off_cnt), containing half the intra-plane elements
    real*8, allocatable :: on_plane_dist(:)    !> @var Contains the distance of each point to the Origin
    real*8, allocatable :: off_plane_dist(:,:) !> @var Contains the distance of each point to the Origin and to the plane
    integer :: on_plane_aux(nn_stages+2)       !> @var on_plane_aux(i) contains the first element of the i-1 -th neighbour stage
    integer :: off_plane_aux(nn_stages+2)
    integer :: on_cnt
    integer :: off_cnt
    real*8 :: dist_tmp
    real*8 :: pos_tmp(3)
    real*8 :: cos_tmp(3)
    real*8 :: cnt(2, nn_stages+2)              !> @var Contains distance and stage informations for all elements in tmp


    nnt = 2 * nn_stages

    ALLOCATE(tmp(3, 2, (2*nnt+1)**3))
    ALLOCATE(dist((2*nnt+1)**3))

    !Generate all neighbours up to 2 * nn_stages the number of wanted neighbours
    i = 1
    do l = -nnt, nnt
      do m = -nnt, nnt
        do n = -nnt, nnt
          pos_tmp = l*a1 + m*a2 + n*a3
          cos_tmp = pos_tmp / sqrt(pos_tmp(1)**2 + pos_tmp(2)**2 + pos_tmp(3)**2)
          dist_tmp = sqrt(pos_tmp(1)**2 + pos_tmp(2)**2 + pos_tmp(3)**2)
          j = i-1
          do while( j >= 1 .and. dist(j) > dist_tmp)
            tmp(:,:, j+1) = tmp(:,:,j)
            dist(j+1) = dist(j)
            j = j-1
          end do
          tmp(:, 1, j+1) = pos_tmp
          tmp(:, 2, j+1) = cos_tmp
          dist(j+1) = dist_tmp
          i = i+1
        end do
      end do
    end do

    ! Count member of each neighbour stage
    ! cnt(1,i) contains the last member of stage i
    ! cnt(2,i) contains the actual distance of the i-th neighbours
    i = 1
    cnt(1,1) = 1.d0
    do j=1, nn_stages+1
      cnt(2,j) = dist(i)
      do while(cnt(2,j) == dist(i))
        i = i+1
      end do
      cnt(1,j+1) = i
    end do

    nn_size = cnt(1, nn_stages+2) - 1
    print *, nn_size

    ! Sort into on and off plane of interest
    allocate(on_plane(3,2,nn_size), off_plane(3,2,nn_size))
    allocate(on_plane_dist(nn_size), off_plane_dist(2, nn_size))

    on_cnt = 0
    off_cnt = 0
    on_plane_aux = 1
    off_plane_aux = 1

    do i=1, nn_size
      dist_tmp = p(1)*tmp(1,1,i) + p(2)*tmp(2,1,i) + p(3)*tmp(3,1,i)
      if(abs(dist_tmp) < 1d-9) then
        on_cnt = on_cnt + 1
        on_plane(:,:, on_cnt) = tmp(:,:, i)
        on_plane_dist(on_cnt) = dist(i)
        do j=1, nn_stages+1
          if( dist(i) <= cnt(2,j) ) then
            on_plane_aux(j+1) = on_plane_aux(j+1) + 1
          end if
        end do
      else if(dist_tmp < 0) then
        off_cnt = off_cnt + 1
        off_plane(:,:, off_cnt) = tmp(:,:, i)
        off_plane_dist(1,off_cnt) = dist(i)
        off_plane_dist(2,off_cnt) = abs(p(1) * tmp(1,1,i) + p(2) * tmp(2,1,i) + p(3) * tmp(3,1,i))
        do j=1, nn_stages + 1
          if( dist(i) <= cnt(2,j) ) then
            off_plane_aux(j+1) = off_plane_aux(j+1) + 1
          end if
        end do
      end if
    end do

    ! Output for debugging
    print *, "On-plane elements"
    print *, on_cnt, on_plane_aux
    do i=1, nn_stages+1
      print *, i-1, "neighbour"
      do j=on_plane_aux(i), on_plane_aux(i+1)-1
        print *, j, on_plane(:,1,j), on_plane(:,2,j)
        print *, on_plane_dist(j)
      end do
    end do
    print *, "Off-plane elements"
    print *, off_cnt, off_plane_aux
    do i=1, nn_stages+1
      print *, i-1, "neighbour"
      do j=off_plane_aux(i), off_plane_aux(i+1)-1
        print *, j, off_plane(:,1,j), off_plane(:,2,j)
        print *, off_plane_dist(1,j), off_plane_dist(2,j)
      end do
    end do


    print *, "Hey", nn_size
    call abort


    return
  end subroutine

  ! SUBROUTINE init_k()
  !   USE global, global_a1 => a1, global_a2 => a2
  !   INTEGER :: i
  !
  !   global_a1 = a1
  !   global_a2 = a2
  !
  !   CALL generate_kpoints()
  !   k_size = eff_nptk
  !
  !   !k_size = 1
  !
  !   ALLOCATE( k(2,k_size) )
  !   DO i=1, k_size
  !    k(:,i) = kbz(i,:)
  !   END DO
  !
  !
  !     !k(1,1) = 0.D0
  !     !k(2,1) = 0.D0
  !   RETURN
  ! END SUBROUTINE


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

  SUBROUTINE green(GF, kp, z, n)
    INTEGER, INTENT(IN) :: n
    DOUBLE COMPLEX, INTENT(IN) :: z
    DOUBLE COMPLEX, INTENT(OUT) :: GF(n,n)
    REAL*8, INTENT(IN) :: kp(2)
    INTEGER :: i, info
    INTEGER :: ipiv(n), lwork, work(n)

    lwork = n
    GF = CMPLX(0.0D0, 0.D0)
    CALL hamilt_0(kp)

    DO i = 1, n
      GF(i,i) = CMPLX(1.0D0, 0.D0)
      H0(i,i) = H0(i,i) - z
    END DO

    H0 = -H0

    CALL ZGETRF(     n, n, H0, n, ipiv,        info)
    CALL ZGETRS('N', n, n, H0, n, ipiv, GF, n, info)

    !print *, G
  END SUBROUTINE

  subroutine generate_chebychev_points(count)
    integer, intent(in) :: count
    integer :: i

    xi_size = count
    allocate(xi(xi_size+1))
    allocate(wi(xi_size+1))
    wi = pi / (xi_size+1.d0)
    do i=0, xi_size
      xi(i) = cos( (2.d0 * i + 1) / (2.d0 * xi_size + 2.d0) * pi )
      wi(i) = sin( (2.d0 * i + 1) / (2.d0 * xi_size + 2.d0) * pi )
    end do
      wi = wi * pi / (xi_size+1)
  end subroutine



  SUBROUTINE ldos(z, EF, Ec)
    DOUBLE COMPLEX, INTENT(IN) :: z
    DOUBLE COMPLEX, ALLOCATABLE :: G_sum(:,:)
    DOUBLE COMPLEX :: E
    REAL*8 :: h, e1, e2, EF, Ec
    REAL*8, allocatable :: occ(:,:)
    INTEGER :: n, i, j

    ALLOCATE(G_sum(H_size,H_size))

    e1 = Ec
    e2 = EF
    n = 1E+3
    h = (e2-e1) / n


    ALLOCATE(occ(H_size, n+1))

    open(unit=20, file='ldos.txt', status='replace', action='write')


    DO i = 0, n
      E = CMPLX( e1+i*h, IMAG(z) )
      G_sum = 0.D0
      DO j=1, k_size
        CALL green(G, k(:,j), E, H_size)
        G_sum = G_sum + G
      END DO
      G_sum = G_sum / k_size
      DO j=1, H_size
        occ(j,i+1) = -IMAG(G_sum(j,j)) / pi
      END DO
      write (20,*) e1+i*h, occ(1,i+1), occ(layers+1, i+1)
    END DO

    close(unit=20)

  end subroutine

  subroutine exp_occ_imag(z, EF)
    double complex, intent(in) :: z
    real*8, intent(in) :: EF
    double complex, allocatable :: G_sum(:,:), G_loc(:,:)
    double complex :: E
    real*8 :: w
    real*8, allocatable :: occ(:)
    integer :: i,j

    allocate(G_sum(H_size, H_size), G_loc(H_size, H_size))
    allocate(occ(H_size))

    G_sum = 0.d0
    G_loc = 0.d0

    open(unit=20, file='occ_imag.txt', status='replace', action='write')

    do i = 0, xi_size
      do j = 1, k_size
        !E = cmplx(EF, imag(z) + (1+xi(i))/(1-xi(i)))
        !w = 2.d0 / (1.d0 - xi(i))**2 * wi(i)
        E = cmplx(EF, imag(z) + 0.5d0 * (1.d0+xi(i))**2 / (1.d0-xi(i)))
        w  = 0.5d0 * (4.d0 / (1.d0 - xi(i))**2 - 1.d0) * wi(i)
        call green(G_loc, k(:,j), E, H_size)

        G_sum = G_sum + G_loc * w
      end do
!      write(20, *) imag(E), real(G_sum(1,1)), real(G_sum(2,2))
      write (20,*) imag(E), 0.5d0 + real(G_sum(1,1)) / (pi * k_size), 0.5d0 + real(G_sum(2,2)) / (pi * k_size)
    end do

    close(unit=20)

    G_sum = G_sum / k_size

    do i = 0, H_size
      occ(i) = 0.5d0 + real(G_sum(i,i)) / pi
    end do
    print *, occ
  end subroutine

  SUBROUTINE exp_occ_real(z, EF, E_cutoff)
    DOUBLE COMPLEX, INTENT(IN) :: z
    REAL*8, INTENT(IN) :: EF, E_cutoff
    DOUBLE COMPLEX, ALLOCATABLE :: G_sum(:,:), G_loc(:,:)
    DOUBLE COMPLEX :: E
    REAL*8 :: h, e1, e2, w
    REAL*8, allocatable :: occ(:)
    integer :: n, i, j

    allocate(G_sum(H_size,H_size))

    e1 = E_cutoff
    e2 = EF


    allocate(occ(H_size))
    allocate(G_loc(H_size,H_size))
    G_sum = 0.0D0
    G_loc = 0.0D0

    open(unit=20, file='occ.txt', status='replace', action='write')

    !$OMP PARALLEL DO PRIVATE(i,j,G_loc, E, H0, t_k) SHARED(G_sum) COLLAPSE(2)
    do i = 0, xi_size
      do j=1, k_size
        !E  = CMPLX( -0.5D0 * (e2-e1) * xi(i) + 0.5D0 * (e1 + e2), IMAG(z) )
        !w  = (e2-e1) * 0.5d0 * wi(i)
        E  = cmplx(EF - 0.5d0 * (1.d0+xi(i))**2 / (1.d0-xi(i)), IMAG(z))
        w  = 0.5d0 * (4.d0 / (1.d0 - xi(i))**2 - 1.d0) * wi(i)
        !E  = cmplx(EF - (1.d0 + xi(i)) / (1.d0 - xi(i)), imag(z))
        !w  = 2.d0 / (1 - xi(i))**2 * wi(i)
        CALL green(G_loc, k(:,j), E, H_size)
        !$OMP CRITICAL
        G_sum =  G_sum + G_loc * w
        !$OMP END CRITICAL
      end do
      !write (20,*) REAL(E), -IMAG(G_sum(1,1)) * (e2-e1)/(pi*2.D0*k_size), -IMAG(G_sum(2,2)) * (e2-e1)/(pi*2.D0*k_size)
      write (20,*) REAL(E), -IMAG(G_sum(1,1)) / (k_size * pi), -IMAG(G_sum(2,2)) / (k_size * pi)
      !write (20,*) REAL(E), -IMAG(G_sum(1,1)) / (pi * k_size), -IMAG(G_sum(2,2)) / (pi * k_size)

    end do
    !$OMP END PARALLEL DO

    G_sum = G_sum / k_size

    DO j=1, H_size
      occ(j) = -1.d0 / pi * IMAG(G_sum(j,j))
    END DO
    print *, occ

    close(unit=20)


  END SUBROUTINE


END MODULE sTB

PROGRAM neighbours
  USE sTB
  IMPLICIT NONE
  INTEGER ::  n, l
  REAL*8  ::  vec1(3), vec2(3), vec3(3), EF, Ec
  real*8 :: p(3)
  DOUBLE COMPLEX :: z = 0.0D0

  n = 3
  l = 1

  ! ! BCC 1 -1 0
  ! p    = [ 1.0d0, -1.0d0,  0.0d0] / sqrt(2.d0)
  ! vec1 = [-0.5d0,  0.5d0,  0.5d0]
  ! vec2 = [ 0.5d0, -0.5d0,  0.5d0]
  ! vec3 = [ 0.5d0,  0.5d0, -0.5d0]

  ! FCC 1 0 0
  p    = [ 0.0d0,  0.0d0,  1.0d0] !/ sqrt(3.d0)
  vec1 = [ 0.0d0,  0.5d0,  0.5d0]
  vec2 = [ 0.5d0,  0.0d0,  0.5d0]
  vec3 = [ 0.5d0,  0.5d0,  0.0d0]

  CALL initialize(n, l, vec1, vec2, vec3, 2,p)
  z = (0.D0, 5.D-3)
  EF = 15.D0
  Ec = -15.D0
  !CALL ldos(z, EF, Ec)
  !CALL exp_occ_imag(z,10.d0)

END PROGRAM neighbours
