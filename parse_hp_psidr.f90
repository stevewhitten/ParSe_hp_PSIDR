!------------------------------------------------------------------------------
! parse_hp_psidr.f90
!
! ParSe v2 plus a high-pathogenicity PS-IDR scan.
!
! For each 25-residue window the program assigns:
!   F  folded
!   D  conventional (non-phase-separating) intrinsic disorder
!   P  phase-separating intrinsic disorder
!
! Disordered windows are also labeled X or Y in a second track using
! Lifson-Sander sheet propensity and Radzicka-Wolfenden vapor-to-octanol
! hydrophobicity. Contiguous X segments (>=20 residues, >=90% X) that
! also fall on the P side of the helix / nu_model boundary are reported
! as high-pathogenicity PS-IDRs.
!
! Cite:
!   Wilson et al., Protein Sci. 32, e4756 (2023)           ParSe v2
!   Ibrahim et al., J. Biol. Chem. 299, 102801 (2023)      training sets
!   Remo et al., under review                              high-pathogenicity scan
!
! Build:
!   gfortran -O2 -o parse_hp_psidr parse_hp_psidr.f90
!
! Run one sequence from the command line:
!   ./parse_hp_psidr ACDEFGHIKLMNPQRSTVWYACDEFGHIKLM
!
! Run a FASTA file (one or many records):
!   ./parse_hp_psidr proteome.fasta
!
! Limits: 20 standard amino acids only; length 25 to 10000 per record.
!------------------------------------------------------------------------------
program parse_hp_psidr
  implicit none

  integer, parameter :: maxn = 10000
  integer, parameter :: win  = 25
  real,    parameter :: pct_cut = 0.90
  real,    parameter :: hydr_cut = 0.08280152
  real,    parameter :: hydr_sd2 = 0.01679414 * 2.0
  real,    parameter :: pd_m = -0.244078945
  real,    parameter :: pd_b =  0.7885823
  real,    parameter :: helix_ps = 0.9327272
  real,    parameter :: nu_ps    = 0.5416
  real,    parameter :: helix_id = 1.022552
  real,    parameter :: nu_id    = 0.5582901
  real,    parameter :: xy_m = 20.0788
  real,    parameter :: xy_b = -20.0254

  real, parameter :: ppii(20) = (/ &
       0.37,0.25,0.30,0.42,0.17,0.13,0.20,0.39,0.56,0.24, &
       0.36,0.27,1.00,0.53,0.38,0.24,0.32,0.39,0.25,0.25 /)
  real, parameter :: helix_sc(20) = (/ &
       1.42,0.73,1.01,1.63,1.16,0.50,1.20,1.12,1.24,1.29, &
       1.21,0.71,0.65,1.02,1.06,0.71,0.78,0.99,1.05,0.67 /)
  real, parameter :: hydr_sc(20) = (/ &
       0.0728, 0.3557,-0.0552,-0.0295, 0.4201,-0.0589, 0.0874, 0.3805,-0.0053, 0.3819, &
       0.1613,-0.0390,-0.0492, 0.0126, 0.0394,-0.0282, 0.0239, 0.2947, 0.4114, 0.3113 /)
  real, parameter :: sheet_sc(20) = (/ &
       0.90,1.24,0.47,0.62,1.23,0.56,1.12,1.54,0.74,1.26, &
       1.09,0.62,0.42,1.18,1.02,0.87,1.30,1.53,1.75,1.68 /)
  real, parameter :: hydr2_sc(20) = (/ &
       1.42,  0.00,  0.00, -9.45, -2.85,  2.39,-11.22,  0.11, -9.60,  0.52, &
      -2.80, -9.67,  0.00, -9.31,-18.60, -5.10, -5.15,  0.81, -8.39, -7.74 /)

  character(len=11000) :: arg
  character(len=80) :: hdr
  character(len=1) :: seq(maxn)
  integer :: narg, ios, nseq, nrec
  logical :: isfile

  narg = command_argument_count()
  if (narg < 1) stop 'no input argument, exiting program'
  call get_command_argument(1, arg)
  if (len_trim(arg) == 0) stop 'no input argument, exiting program'
  inquire(file=trim(arg), exist=isfile)
  if (isfile) then
     open(unit=10, file=trim(arg), status='old', action='read', iostat=ios)
     if (ios /= 0) stop 'could not open file'
     nrec = 0
     do
        call read_fasta_record(10, hdr, seq, nseq, ios)
        if (ios /= 0) exit
        nrec = nrec + 1
        call analyze_one(hdr, seq, nseq)
     end do
     close(10)
     if (nrec == 0) stop 'no sequences found in file'
  else
     call load_raw(trim(arg), seq, nseq)
     call analyze_one(' ', seq, nseq)
  end if

contains

  integer function aai(c)
    character(len=1), intent(in) :: c
    select case (c)
    case ('A'); aai = 1
    case ('C'); aai = 2
    case ('D'); aai = 3
    case ('E'); aai = 4
    case ('F'); aai = 5
    case ('G'); aai = 6
    case ('H'); aai = 7
    case ('I'); aai = 8
    case ('K'); aai = 9
    case ('L'); aai = 10
    case ('M'); aai = 11
    case ('N'); aai = 12
    case ('P'); aai = 13
    case ('Q'); aai = 14
    case ('R'); aai = 15
    case ('S'); aai = 16
    case ('T'); aai = 17
    case ('V'); aai = 18
    case ('W'); aai = 19
    case ('Y'); aai = 20
    case default; aai = 0
    end select
  end function aai

  subroutine load_raw(s, seq, nseq)
    character(len=*), intent(in) :: s
    character(len=1), intent(out) :: seq(:)
    integer, intent(out) :: nseq
    integer :: i, ic
    character(len=1) :: c
    nseq = 0
    do i = 1, len(s)
       c = s(i:i)
       if (c == ' ') exit
       ic = iachar(c)
       if (ic >= 97 .and. ic <= 122) c = achar(ic - 32)
       nseq = nseq + 1
       if (nseq > maxn) stop 'input sequence is too long'
       seq(nseq) = c
    end do
  end subroutine load_raw

  subroutine read_fasta_record(u, hdr, seq, nseq, ios)
    integer, intent(in) :: u
    character(len=*), intent(out) :: hdr
    character(len=1), intent(out) :: seq(:)
    integer, intent(out) :: nseq, ios
    character(len=4096) :: line
    integer :: i, ic
    character(len=1) :: c
    logical :: started
    hdr = ' '
    nseq = 0
    started = .false.
    do
       read(u, '(A)', iostat=ios) line
       if (ios /= 0) then
          if (started) ios = 0
          return
       end if
       if (len_trim(line) == 0) cycle
       if (line(1:1) == '>') then
          if (started) then
             backspace(u)
             ios = 0
             return
          end if
          hdr = line(2:min(len(hdr)+1, len_trim(line)))
          started = .true.
          cycle
       end if
       if (.not. started) then
          hdr = 'unnamed'
          started = .true.
       end if
       do i = 1, len_trim(line)
          c = line(i:i)
          if (c == ' ') cycle
          ic = iachar(c)
          if (ic >= 97 .and. ic <= 122) c = achar(ic - 32)
          nseq = nseq + 1
          if (nseq > maxn) stop 'input sequence is too long'
          seq(nseq) = c
       end do
    end do
  end subroutine read_fasta_record

  subroutine window_props(cnt, nres, hydr, helix, nu, sheet, hydr2)
    integer, intent(in) :: cnt(20), nres
    real, intent(out) :: hydr, helix, nu, sheet, hydr2
    integer :: k, q
    real :: fppii, vexp, rh
    hydr = 0.0
    helix = 0.0
    sheet = 0.0
    hydr2 = 0.0
    fppii = 0.0
    q = abs((cnt(3) + cnt(4)) - (cnt(9) + cnt(15)))
    do k = 1, 20
       hydr  = hydr  + real(cnt(k)) * hydr_sc(k)
       helix = helix + real(cnt(k)) * helix_sc(k)
       sheet = sheet + real(cnt(k)) * sheet_sc(k)
       hydr2 = hydr2 + real(cnt(k)) * hydr2_sc(k)
       fppii = fppii + real(cnt(k)) * ppii(k)
    end do
    hydr  = hydr  / real(nres)
    helix = helix / real(nres)
    sheet = sheet / real(nres)
    hydr2 = hydr2 / real(nres)
    fppii = fppii / real(nres)
    if (fppii == 1.0) fppii = 0.98
    vexp = 0.503 - 0.11 * log(1.0 - fppii)
    rh = 2.16 * (real(4 * nres)**vexp) + 0.26 * real(4 * q) - 0.29 * sqrt(real(4 * nres))
    nu = log(rh / 2.16) / log(real(4 * nres))
  end subroutine window_props

  subroutine grow_runs(lab, npep, want, s0, s1, nrun)
    ! Original grow-from-20 / >=90% rule (labels 20/30, 40/50, 60/70, 96/97).
    character(len=1), intent(in) :: lab(:), want
    integer, intent(in) :: npep
    integer, intent(out) :: s0(:), s1(:), nrun
    integer :: i, j, count_p, count_w, region, pstart, pend
    real :: percent_p
    nrun = 0
    i = 1
    do while (i + 19 <= npep)
       count_p = 0
       count_w = 0
       region = 0
       do j = i, i + 19
          count_w = count_w + 1
          if (lab(j) == want) count_p = count_p + 1
       end do
       do
          percent_p = real(count_p) / real(count_w)
          if (percent_p < pct_cut) exit
          region = 1
          pstart = i
          pend = j
          if (j >= npep) exit
          j = j + 1
          count_w = count_w + 1
          if (lab(j) == want) count_p = count_p + 1
       end do
       if (region == 1) then
          nrun = nrun + 1
          s0(nrun) = pstart
          s1(nrun) = pend
          i = j
       else
          i = i + 1
       end if
    end do
  end subroutine grow_runs

  subroutine split_embedded(a0, a1, na, b0, b1, nb)
    ! Original labels 90-95: if a B-run lies strictly inside an A-run,
    ! trim A or split A into two runs.
    integer, intent(inout) :: a0(:), a1(:), na
    integer, intent(in)    :: b0(:), b1(:), nb
    integer :: i, j, nterm, cterm
    logical :: restart
    restart = .true.
    do while (restart)
       restart = .false.
       outer: do i = 1, na
          do j = 1, nb
             if (b0(j) > a0(i) .and. b1(j) < a1(i)) then
                nterm = b0(j) - a0(i)
                cterm = a1(i) - b1(j)
                if (nterm > 20 .and. cterm <= 20) then
                   a1(i) = b0(j) - 1
                else if (nterm <= 20 .and. cterm > 20) then
                   a0(i) = b1(j) + 1
                else if (nterm > 20 .and. cterm > 20) then
                   na = na + 1
                   a1(na) = a1(i)
                   a1(i)  = b0(j) - 1
                   a0(na) = b1(j) + 1
                   restart = .true.
                   exit outer
                end if
             end if
          end do
       end do outer
    end do
  end subroutine split_embedded

  subroutine build_order(p0, np, d0, nd, f0, nf, domain, ndom)
    integer, intent(in) :: p0(:), np, d0(:), nd, f0(:), nf
    integer, intent(out) :: domain(:), ndom
    integer :: i, j, low
    ndom = np + nd + nf
    if (ndom == 0) return
    low = 1000000
    do i = 1, np
       if (p0(i) < low) low = p0(i)
    end do
    do i = 1, nd
       if (d0(i) < low) low = d0(i)
    end do
    do i = 1, nf
       if (f0(i) < low) low = f0(i)
    end do
    domain(1) = low
    do j = 2, ndom
       low = 1000000
       do i = 1, np
          if (p0(i) < low .and. p0(i) > domain(j-1)) low = p0(i)
       end do
       do i = 1, nd
          if (d0(i) < low .and. d0(i) > domain(j-1)) low = d0(i)
       end do
       do i = 1, nf
          if (f0(i) < low .and. f0(i) > domain(j-1)) low = f0(i)
       end do
       domain(j) = low
    end do
  end subroutine build_order

  subroutine trim_adjacent(p0, p1, np, d0, d1, nd, f0, f1, nf, domain, ndom)
    ! Original overlap repair: if domain i ends at or past the start of
    ! domain i+1, give the first domain the extra
    ! int((1-pct_cut)*20/2) residues of the overlap and slide the next start.
    integer, intent(inout) :: p0(:), p1(:), np, d0(:), d1(:), nd, f0(:), f1(:), nf
    integer, intent(inout) :: domain(:), ndom
    integer :: i, j, jj, extra
    extra = int((1.0 - pct_cut) * 20.0 / 2.0)
    if (ndom <= 1) return
    do i = 1, ndom - 1
       do j = 1, np
          if (p0(j) == domain(i) .and. p1(j) >= domain(i+1)) then
             p1(j) = domain(i+1) + extra
             do jj = 1, nd
                if (d0(jj) == domain(i+1)) then
                   d0(jj) = p1(j) + 1
                   domain(i+1) = d0(jj)
                end if
             end do
             do jj = 1, nf
                if (f0(jj) == domain(i+1)) then
                   f0(jj) = p1(j) + 1
                   domain(i+1) = f0(jj)
                end if
             end do
          end if
       end do
       do j = 1, nd
          if (d0(j) == domain(i) .and. d1(j) >= domain(i+1)) then
             d1(j) = domain(i+1) + extra
             do jj = 1, np
                if (p0(jj) == domain(i+1)) then
                   p0(jj) = d1(j) + 1
                   domain(i+1) = p0(jj)
                end if
             end do
             do jj = 1, nf
                if (f0(jj) == domain(i+1)) then
                   f0(jj) = d1(j) + 1
                   domain(i+1) = f0(jj)
                end if
             end do
          end if
       end do
       do j = 1, nf
          if (f0(j) == domain(i) .and. f1(j) >= domain(i+1)) then
             f1(j) = domain(i+1) + extra
             do jj = 1, nd
                if (d0(jj) == domain(i+1)) then
                   d0(jj) = f1(j) + 1
                   domain(i+1) = d0(jj)
                end if
             end do
             do jj = 1, np
                if (p0(jj) == domain(i+1)) then
                   p0(jj) = f1(j) + 1
                   domain(i+1) = p0(jj)
                end if
             end do
          end if
       end do
    end do
  end subroutine trim_adjacent

  subroutine analyze_one(hdr, seq, npep)
    character(len=*), intent(in) :: hdr
    character(len=1), intent(in) :: seq(:)
    integer, intent(in) :: npep
    character(len=1) :: lab(maxn), lab2(maxn)
    integer :: cnt(20), i, j, mid
    integer :: p0(maxn), p1(maxn), np
    integer :: d0(maxn), d1(maxn), nd
    integer :: f0(maxn), f1(maxn), nf
    integer :: x0(maxn), x1(maxn), nx
    integer :: xp0(maxn), xp1(maxn), nxp
    integer :: domain(maxn), ndom
    real :: hydr, helix, nu, sheet, hydr2
    real :: id_dist, ps_dist, m, b, x, y, p_dist_sum

    if (npep < win) then
       write(*,*) 'could not parse sequence'
       write(*,*) 'input sequence is too short'
       return
    end if
    do i = 1, npep
       if (aai(seq(i)) == 0) then
          write(*,*) 'could not parse sequence'
          write(*,*) 'sequence contains noncommon amino acid type'
          return
       end if
       lab(i) = ' '
       lab2(i) = ' '
    end do

    m = -1.0 / pd_m
    b = nu_id - m * helix_id
    x = (b - pd_b) / (pd_m - m)
    y = m * x + b
    id_dist = sqrt((helix_id - x)**2 + (nu_id - y)**2)
    b = nu_ps - m * helix_ps
    x = (b - pd_b) / (pd_m - m)
    y = m * x + b
    ps_dist = sqrt((helix_ps - x)**2 + (nu_ps - y)**2)

    p_dist_sum = 0.0
    do i = 1, npep - win + 1
       mid = i + win / 2
       cnt = 0
       do j = i, i + win - 1
          cnt(aai(seq(j))) = cnt(aai(seq(j))) + 1
       end do
       call window_props(cnt, win, hydr, helix, nu, sheet, hydr2)
       if (hydr >= hydr_cut) then
          lab(mid) = 'F'
          lab2(mid) = 'F'
          cycle
       end if
       if (hydr2 > (xy_m * sheet + xy_b)) then
          lab2(mid) = 'X'
       else
          lab2(mid) = 'Y'
       end if
       m = -1.0 / pd_m
       b = nu - m * helix
       x = (b - pd_b) / (pd_m - m)
       y = m * x + b
       if (((nu - pd_b) / pd_m) <= helix) then
          lab(mid) = 'D'
       else
          lab(mid) = 'P'
          p_dist_sum = p_dist_sum + sqrt((helix - x)**2 + (nu - y)**2) / ps_dist
       end if
    end do

    do j = 1, win / 2
       lab(j) = lab(win / 2 + 1)
       lab2(j) = lab2(win / 2 + 1)
    end do
    do j = npep - win / 2 + 1, npep
       lab(j) = lab(npep - win / 2)
       lab2(j) = lab2(npep - win / 2)
    end do

    call grow_runs(lab,  npep, 'P', p0, p1, np)
    call grow_runs(lab,  npep, 'D', d0, d1, nd)
    call grow_runs(lab,  npep, 'F', f0, f1, nf)
    call grow_runs(lab2, npep, 'X', x0, x1, nx)

    ! Embedded shorter domain inside a longer domain of another class.
    call split_embedded(p0, p1, np, d0, d1, nd)  ! D inside P
    call split_embedded(p0, p1, np, f0, f1, nf)  ! F inside P
    call split_embedded(d0, d1, nd, p0, p1, np)  ! P inside D
    call split_embedded(d0, d1, nd, f0, f1, nf)  ! F inside D
    call split_embedded(f0, f1, nf, d0, d1, nd)  ! D inside F
    call split_embedded(f0, f1, nf, p0, p1, np)  ! P inside F

    call build_order(p0, np, d0, nd, f0, nf, domain, ndom)
    call trim_adjacent(p0, p1, np, d0, d1, nd, f0, f1, nf, domain, ndom)

    nxp = 0
    do i = 1, nx
       cnt = 0
       do j = x0(i), x1(i)
          cnt(aai(seq(j))) = cnt(aai(seq(j))) + 1
       end do
       call window_props(cnt, x1(i) - x0(i) + 1, hydr, helix, nu, sheet, hydr2)
       if (((nu - pd_b) / pd_m) > helix) then
          nxp = nxp + 1
          xp0(nxp) = x0(i)
          xp1(nxp) = x1(i)
       end if
    end do

    write(*,*)
    if (len_trim(hdr) > 0) write(*,'(A,A)') '> ', trim(hdr)
    write(*,'(A,I0)') 'length ', npep
    write(*,'(A,F10.3)') 'PS potential (summed P-window classifier distance) ', p_dist_sum
    write(*,'(A)') 'ParSe labels (F/D/P):'
    write(*,'(10000A1)') (lab(j), j = 1, npep)
    write(*,'(A)') 'pathogenic-site ID labels (X) vs other ID (Y) or folded (F):'
    write(*,'(10000A1)') (lab2(j), j = 1, npep)
    write(*,*)
    write(*,'(A)') 'domains (>=20 residues, >=90% one label)'
    do i = 1, ndom
       do j = 1, np
          if (p0(j) == domain(i)) write(*,'(A,I0,A,I0,A,I0)') &
               'PS-ID (P)     first ', p0(j), '  last ', p1(j), '  length ', p1(j)-p0(j)+1
       end do
       do j = 1, nd
          if (d0(j) == domain(i)) write(*,'(A,I0,A,I0,A,I0)') &
               'nonPS-ID (D)  first ', d0(j), '  last ', d1(j), '  length ', d1(j)-d0(j)+1
       end do
       do j = 1, nf
          if (f0(j) == domain(i)) write(*,'(A,I0,A,I0,A,I0)') &
               'folded (F)    first ', f0(j), '  last ', f1(j), '  length ', f1(j)-f0(j)+1
       end do
    end do
    write(*,*)
    write(*,'(A,I0)') 'high-pathogenicity PS-IDRs = ', nxp
    if (nxp > 0) then
       write(*,'(A)') 'index   first    last   length'
       do i = 1, nxp
          write(*,'(I5,2X,I7,2X,I7,2X,I7)') i, xp0(i), xp1(i), xp1(i)-xp0(i)+1
       end do
    end if

  end subroutine analyze_one

end program parse_hp_psidr
