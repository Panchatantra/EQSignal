module basic

    use, intrinsic :: iso_c_binding
    implicit none
    include "fftw3.fi"
    real(8), parameter :: PI     = 3.141592653589793D0
    real(8), parameter :: TWO_PI = 6.283185307179586D0

    interface lstsq
        module procedure leastsqs, leastsqm
    end interface ! lstsq

    interface polyval
        module procedure polyvals, polyvalv
    end interface ! polyval

contains

function mean(x,n) bind(c)
! 平均值
    integer, intent(in) :: n
    real(8), intent(in) :: x(n)
    real(8) :: mean

    mean = sum(x)/dble(n)
end function mean

function rms(x,n) bind(c)
! 均方根值
    integer, intent(in) :: n
    real(8), intent(in) :: x(n)
    real(8) :: rms

    rms = sqrt(sum(x*x)/dble(n))
end function rms

function peak(x,n) bind(c)
! 峰值
    integer, intent(in) :: n
    real(8), intent(in) :: x(n)
    real(8) :: peak

    peak = maxval(abs(x),dim=1)
end function peak

function peakloc(x,n) bind(c)
! 峰值所在位置
    integer, intent(in) :: n
    real(8), intent(in) :: x(n)
    integer :: peakloc

    peakloc = maxloc(abs(x),dim=1)
end function peakloc

subroutine norm(a,n) bind(c)
! 按峰值归一
    integer, intent(in) :: n
    real(8), intent(inout) :: a(n)

    integer :: i
    real(8) :: pk, ai

    pk = 0.d0

    do i = 1, n, 1
        ai = abs(a(i))
        if ( ai>pk ) pk = ai
    end do

    do i = 1, n, 1
        a(i) = a(i)/pk
    end do

    return

end subroutine norm

function nextpow2(n) bind(c)
! 不小于n的2的整数次幂
    integer, intent(in) :: n
    integer :: nextpow2

    nextpow2 = 1

    do while ( nextpow2<n )
        nextpow2 = nextpow2*2
    end do
    return

end function nextpow2

function nextpow(n,base) bind(c)
! 不小于n的base的整数次幂
    integer, intent(in) :: n, base
    integer :: nextpow

    nextpow = 1

    do while ( nextpow<n )
        nextpow = nextpow*base
    end do
    return

end function nextpow

function dpower(x,n,m) bind(c)
! 对x^n求m阶导数后的值
    real(8), intent(in) :: x
    integer, intent(in) :: n, m
    real(8) :: dpower

    integer :: k, l

    if ( m>n ) then
        dpower = 0.d0
        return
    end if

    k = n
    l = m
    dpower = x**(n-m)

    do while ( l>0 )
        dpower = dpower*k
        l = l - 1
        k = k - 1
    end do
end function dpower

function trapz(y,n,dx,y0) bind(c)
! 等间隔dx离散数据梯形法数值积分，初值为y0
    implicit none
    real*8, intent(in) :: dx, y0
    integer, intent(in) :: n
    real*8, intent(in) :: y(n)

    real*8 :: trapz

    integer :: i
    real*8 :: r

    r = 0.d0

    !$OMP PARALLEL PRIVATE(i), SHARED(y,dx,n)
    !$OMP DO REDUCTION(+:r)

    do i = 2, n, 1
        r = r + dx*(y(i)+y(i-1))
    end do
    
    !$OMP END DO
    !$OMP END PARALLEL

    trapz = y0+r*0.5
    return

end function trapz


subroutine cumtrapz(y,z,n,dx,z0) bind(c)
! 等间隔dx离散数据梯形法累积数值积分，结果保存在数组z中，初值为z0

    real*8, intent(in) :: dx, z0
    integer, intent(in) :: n
    real*8, intent(in) :: y(n)
    real*8, intent(out) :: z(n)

    integer :: i
    real*8 :: r

    z(1) = z0
    r = 0.D0
    do i = 2, n, 1
        r = r + dx*(y(i)+y(i-1))*0.5
        z(i) = z0 + r
    end do
    return

end subroutine cumtrapz

subroutine ariasIntensity(a,Ia,n) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: ariasIntensity
! 计算Arias强度，结果记录在数组Ia中
    integer, intent(in) :: n
    real*8, intent(in) :: a(n)
    real*8, intent(out) :: Ia(n)

    call cumtrapz(a*a,Ia,n,1.d0,0.d0)
    Ia = Ia/Ia(n)
    return

end subroutine ariasIntensity

subroutine fftfreqs(Nfft,fs,freqs) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: fftfreqs
! 离散傅里叶变换的频率轴刻度值，fs为采样频率，结果记录在数组freqs中
    integer, intent(in) :: Nfft
    real(8), intent(in) :: fs
    real(8), intent(out) :: freqs(Nfft)

    real(8) :: fn, df
    integer :: i, l

    fn = 0.5D0*fs
    df = fs/dble(Nfft)
    freqs = 0.0D0

    if ( mod(Nfft, 2) == 0 ) then
        l = Nfft/2-1
        freqs(2:Nfft/2) = [ ( (dble(i)*df), i=1,l ) ]
        freqs(Nfft/2+1:Nfft) = [ ( -fn+(dble(i-1)*df), i=1,l+1 ) ]
    else
        l = (Nfft+1)/2-1
        freqs(2:(Nfft+1)/2) = [ ( (dble(i)*df), i=1,l ) ]
        freqs((Nfft+1)/2+1:Nfft) = [ ( -(dble(l-i+1)*df), i=1,l ) ]
    end if

    return
end subroutine fftfreqs

subroutine acc2vd(a,v,d,n,dt,v0,d0) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: acc2vd
! 加速度积分为速度、位移
    real(8), intent(in) :: dt,v0,d0
    integer(4), intent(in) :: n
    real(8), intent(in) :: a(n)
    real(8), intent(out) :: v(n), d(n)

    integer(4) :: i
    real(8) :: vi, di

    v(1) = v0
    d(1) = d0

    vi = 0.D0
    di = 0.D0
    do i = 2, n, 1
        vi = vi + dt*(a(i)+a(i-1))*0.5
        v(i) = v0 + vi
        di = di + dt*(v(i)+v(i-1))*0.5
        d(i) = d0 + di
    end do
    return
end subroutine acc2vd

subroutine ratacc2vd(a,v,d,n,dt,v0,d0) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: ratacc2vd
! 加速度积分为速度、位移，采用预估的初速度和初位移
    real(8), intent(in) :: dt,v0,d0
    integer(4), intent(in) :: n
    real(8), intent(in) :: a(n)
    real(8), intent(out) :: v(n), d(n)

    real(8) :: vn, dn
    integer :: i

    call acc2vd(a,v,d,n,dt,v0,d0)
    vn = sum(v)/dble(n)
    dn = 0.0d0
    do i = 1, n, 1
        dn = dn + d(i) - vn*dble(i-1)*dt
    end do
    dn = dn/dble(n)

    call acc2vd(a,v,d,n,dt,v0-vn,d0-dn)

end subroutine ratacc2vd

subroutine leastsqs(a,b,m,n) bind(c)
! 最小二乘法解方程组 A*x = b, b为一维数组
    integer, intent(in) :: m, n
    real*8, intent(in) :: a(m,n)
    real*8, intent(inout) :: b(m)

    integer :: i, nrhs, lda, ldb, lwork
    integer :: rank, info
    real*8 :: twork(1)
    integer :: tiwork(1)
    real*8, allocatable :: s(:)
    real*8, allocatable :: work(:)
    integer, allocatable :: iwork(:)

    ! m = size(a, dim=1)
    ! n = size(a, dim=2)
    nrhs = 1

    lda = max(1,m)
    ldb = max(1,m,n)

    allocate(s(min(m,n)))

    lwork = -1

    call dgelsd(m, n, nrhs, a, lda, b, ldb, s, -1.d0, &
        rank, twork, lwork, tiwork, info)

    lwork = int(twork(1))
    allocate(work(lwork))
    allocate(iwork(tiwork(1)))

    call dgelsd(m, n, nrhs, a, lda, b, ldb, s, -1.d0, &
        rank, work, lwork, iwork, info)

    deallocate(s)
    deallocate(work)
    deallocate(iwork)

    return

end subroutine leastsqs

subroutine leastsqm(a,b,m,n,nrhs) bind(c)
! 最小二乘法解方程组 A*x = b, b为二维数组
    real*8, intent(in) :: a(m,n)
    real*8, intent(inout) :: b(m,nrhs)
    integer, intent(in) :: m, n, nrhs

    integer :: i, lda, ldb, lwork
    integer :: rank, info
    real*8 :: twork(1)
    integer :: tiwork(1)
    real*8, allocatable :: s(:)
    real*8, allocatable :: work(:)
    integer,allocatable :: iwork(:)

    ! m = size(a, dim=1)
    ! n = size(a, dim=2)
    ! nrhs = size(b, dim=2)

    lda = max(1,m)
    ldb = max(1,m,n)

    allocate(s(min(m,n)))

    lwork = -1

    call dgelsd(m, n, nrhs, a, lda, b, ldb, s, -1.d0, &
        rank, twork, lwork, tiwork, info)

    lwork = int(twork(1))
    allocate(work(lwork))
    allocate(iwork(tiwork(1)))

    call dgelsd(m, n, nrhs, a, lda, b, ldb, s, -1.d0, &
        rank, work, lwork, iwork, info)

    deallocate(s)
    deallocate(work)
    deallocate(iwork)

    return

end subroutine leastsqm

subroutine error(y,y0,n,aerror,merror) bind(c)
! 计算数组y和y0的相对误差（不包括首尾元素），aerror 为平均误差，merror 为最大误差
    integer, intent(in) :: n
    real(8), intent(in) :: y(n), y0(n)
    real(8), intent(out) :: aerror, merror

    real(8), allocatable :: e(:)

    allocate(e(n-2))

    e = (y(2:n-1) - y0(2:n-1))/y0(2:n-1)

    aerror = sqrt(sum(e*e)/dble(n-2))
    merror = maxval(abs(e), dim=1)

    deallocate(e)

end subroutine error

subroutine errora(y,y0,n,aerror,merror) bind(c)
! 计算数组y和y0的相对误差（包括所有元素），aerror 为平均误差，merror 为最大误差
    integer, intent(in) :: n
    real(8), intent(in) :: y(n), y0(n)
    real(8), intent(out) :: aerror, merror

    real(8), allocatable :: e(:)

    allocate(e(n))

    e = (y - y0)/y0

    aerror = sqrt(sum(e*e)/dble(n))
    merror = maxval(abs(e), dim=1)

    deallocate(e)

end subroutine errora

subroutine errora_endur(y,y0,n,m,aerror,merror) bind(c)
! 计算矩阵y和y0的相对误差（包括所有元素），aerror 为平均误差，merror 为最大误差
    integer, intent(in) :: n, m
    real(8), intent(in) :: y(n,m), y0(n,m)
    real(8), intent(out) :: aerror, merror

    real(8) :: eij

    integer :: i,j

    aerror = 0.D0
    merror = 0.D0
    do i = 1, n, 1
        do j = 1, m, 1
            eij = (y(i,j) - y0(i,j))/y0(i,j)
            aerror = aerror + eij*eij
            if ( merror < abs(eij) ) merror = abs(eij)
        end do
    end do

    aerror = aerror/dble(n*m)

end subroutine errora_endur

subroutine incrlininterp(x,y,n,xi,yi,ni) bind(c)
! 线性插值（x, xi为递增数组）
    integer, intent(in) :: n, ni
    real(8), intent(in) :: x(n), y(n), xi(ni)
    real(8), intent(out) :: yi(ni)

    integer :: i, j
    real(8) :: xc, yc, xp, yp, slope

    xp = x(1)
    yp = y(1)
    j = 1
    do i = 2, n, 1
        xc = x(i)
        yc = y(i)
        slope = (yc - yp)/(xc - xp)
        do while (j<=ni)
            if (xi(j)>xc) exit
            yi(j) = yp + slope*(xi(j)-xp)
            j = j + 1
        end do
        xp = xc
        yp = yc
    end do

end subroutine incrlininterp

subroutine incrfindfirst(a,n,x,i) bind(c)
! 找到数组a从头开始第一个超过x的值的位置（a的包络线为递增趋势）
    integer, intent(in) :: n
    real(8), intent(in) :: a(n), x
    integer, intent(out) :: i

    i = 1
    do while ( a(i)<x .and. i<=n )
        i = i + 1
    end do

    return

end subroutine incrfindfirst

subroutine incrfindlast(a,n,x,i) bind(c)
! 找到数组a从尾开始第一个超过x的值的位置（a的包络线为递增趋势）
    integer, intent(in) :: n
    real(8), intent(in) :: a(n), x
    integer, intent(out) :: i

    i = n
    do while ( a(i)>x .and. i>=1 )
        i = i - 1
    end do

    return

end subroutine incrfindlast

subroutine decrlininterp(x,y,n,xi,yi,ni) bind(c)
! 线性插值（x, xi为递减数组）
    integer, intent(in) :: n, ni
    real(8), intent(in) :: x(n), y(n), xi(ni)
    real(8), intent(out) :: yi(ni)

    integer :: i, j
    real(8) :: xc, yc, xp, yp, slope

    xp = x(n)
    yp = y(n)
    j = 1
    do i = 1, n-1, 1
        xc = x(n-i)
        yc = y(n-i)
        slope = (yc - yp)/(xc - xp)
        do while (j<=ni)
            if (xi(j)<xc) exit
            yi(j) = yp + slope*(xi(j)-xp)
            j = j + 1
        end do
        xp = xc
        yp = yc
    end do

end subroutine decrlininterp

subroutine decrfindfirst(a,n,x,i) bind(c)
! 找到数组a从头开始第一个超过x的值的位置（a的包络线为递减趋势）
    integer, intent(in) :: n
    real(8), intent(in) :: a(n), x
    integer, intent(out) :: i

    i = 1
    do while ( a(i)>x .and. i<=n-1 )
        i = i + 1
    end do

    return

end subroutine decrfindfirst

subroutine decrfindlast(a,n,x,i) bind(c)
! 找到数组a从尾开始第一个超过x的值的位置（a的包络线为递减趋势）
    integer, intent(in) :: n
    real(8), intent(in) :: a(n), x
    integer, intent(out) :: i

    i = n
    do while ( a(i)<x .and. i>=1 )
        i = i - 1
    end do

    return

end subroutine decrfindlast

subroutine targetdc(a,td,n,tp,ntp,ph,pl,dt,v0,d0) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: targetdc
! 基线调整，以规定位移为目标
    integer(4), intent(in) :: n
    real(8), intent(in) :: dt,v0,d0
    real(8), intent(in) :: td(n)
    real(8), intent(out) :: a(n)
    integer(4), intent(in) :: ntp,ph,pl
    integer(4), intent(in) :: tp(ntp)

    integer :: i,j,k,pn
    real(8), allocatable :: M(:,:), b(:), tc(:)
    real(8), allocatable :: v(:), d(:), t(:)

    pn = ph - pl + 1
    if ( pn>ntp ) pn = ntp

    allocate(M(ntp,pn))
    allocate(b(ntp))
    allocate(tc(ntp))
    allocate(t(n))
    allocate(v(n))
    allocate(d(n))

    t = [ (((i-1)*dt),i=1,n) ]
    call acc2vd(a,v,d,n,dt,v0,d0)

    do k = 1, ntp, 1
        j = tp(k)
        tc(k) = t(j)
        b(k)  = td(j)-d(j)
        do i = pl, pl+pn-1, 1
            M(k,i-pl+1) = dpower(tc(k),i,0)
        end do
    end do

    call leastsqs(M,b,ntp,pn)

    do k = 1, n, 1
        do i = pl, pl+pn-1, 1
            a(k) = a(k)+b(i-pl+1)*dpower(t(k),i,2)
        end do
    end do

    deallocate(M)
    deallocate(b)
    deallocate(tc)
    deallocate(t)
    deallocate(v)
    deallocate(d)

    return

end subroutine targetdc

subroutine targetdvc(a,td,tv,n,tp,ntp,ph,pl,dt,v0,d0) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: targetdvc
! 基线调整，以规定速度、位移为目标
    integer(4), intent(in) :: n
    real(8), intent(in) :: dt,v0,d0
    real(8), intent(in) :: td(n),tv(n)
    real(8), intent(out) :: a(n)
    integer(4), intent(in) :: ntp,ph,pl
    integer(4), intent(in) :: tp(ntp)

    integer :: i,j,k,pn
    real(8), allocatable :: M(:,:), b(:), tc(:)
    real(8), allocatable :: v(:), d(:), t(:)

    pn = (ph - pl + 1)*2
    if ( pn>2*ntp ) pn = 2*ntp

    allocate(M(2*ntp,pn))
    allocate(b(2*ntp))
    allocate(tc(ntp))
    allocate(t(n))
    allocate(v(n))
    allocate(d(n))

    t = [ (((i-1)*dt),i=1,n) ]
    call acc2vd(a,v,d,n,dt,v0,d0)

    do k = 1, ntp, 1
        j = tp(k)
        tc(k) = t(j)
        b(k)  = td(j)-d(j)
        b(k+ntp)  = tv(j)-v(j)
        do i = pl, pl+pn-1, 1
            M(k,i-pl+1) = dpower(tc(k),i,0)
            M(k+ntp,i-pl+1) = dpower(tc(k),i,1)
        end do
    end do

    call leastsqs(M,b,2*ntp,pn)

    do k = 1, n, 1
        do i = pl, pl+pn-1, 1
            a(k) = a(k)+b(i-pl+1)*dpower(t(k),i,2)
        end do
    end do

    deallocate(M)
    deallocate(b)
    deallocate(tc)
    deallocate(t)
    deallocate(v)
    deallocate(d)

    return

end subroutine targetdvc

subroutine targetdvac(a,td,tv,ta,n,tp,ntp,ph,pl,dt,v0,d0) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: targetdvac
! 基线调整，以规定加速度、速度、位移为目标
    integer(4), intent(in) :: n
    real(8), intent(in) :: dt,v0,d0
    real(8), intent(in) :: td(n),tv(n),ta(n)
    real(8), intent(out) :: a(n)
    integer(4), intent(in) :: ntp,ph,pl
    integer(4), intent(in) :: tp(ntp)

    integer :: i,j,k,pn
    real(8), allocatable :: M(:,:), b(:), tc(:)
    real(8), allocatable :: v(:), d(:), t(:)

    pn = (ph - pl + 1)*3
    if ( pn>3*ntp ) pn = 3*ntp

    allocate(M(3*ntp,pn))
    allocate(b(3*ntp))
    allocate(tc(ntp))
    allocate(t(n))
    allocate(v(n))
    allocate(d(n))

    t = [ (((i-1)*dt),i=1,n) ]
    call acc2vd(a,v,d,n,dt,v0,d0)

    do k = 1, ntp, 1
        j = tp(k)
        tc(k) = t(j)
        b(k)  = td(j)-d(j)
        b(k+ntp)  = tv(j)-v(j)
        b(k+2*ntp)  = ta(j)-a(j)
        do i = pl, pl+pn-1, 1
            M(k,i-pl+1) = dpower(tc(k),i,0)
            M(k+ntp,i-pl+1) = dpower(tc(k),i,1)
            M(k+2*ntp,i-pl+1) = dpower(tc(k),i,2)
        end do
    end do

    call leastsqs(M,b,3*ntp,pn)

    do k = 1, n, 1
        do i = pl, pl+pn-1, 1
            a(k) = a(k)+b(i-pl+1)*dpower(t(k),i,2)
        end do
    end do

    deallocate(M)
    deallocate(b)
    deallocate(tc)
    deallocate(t)
    deallocate(v)
    deallocate(d)

    return

end subroutine targetdvac

subroutine fft(in,out,n) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: fft
! 快速离散傅里叶变换
    integer, intent(in) :: n
    real*8, intent(inout) :: in(n)
    complex*16, intent(inout) :: out(n)

    type(c_ptr) :: plan

    plan = fftw_plan_dft_r2c_1d(n, in, out, fftw_estimate)
    call fftw_execute_dft_r2c(plan, in, out)
    call fftw_destroy_plan(plan)

    return
end subroutine fft

subroutine ifft(in,out,n) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: ifft
! 快速离散傅里叶拟变换
    integer, intent(in) :: n
    complex*16, intent(inout) :: in(n)
    real*8, intent(inout) :: out(n)

    type(c_ptr) :: plan

    plan = fftw_plan_dft_c2r_1d(n, in, out, fftw_estimate)
    call fftw_execute_dft_c2r(plan, in, out)
    call fftw_destroy_plan(plan)
    out = out/dble(n)

    return
end subroutine ifft

subroutine hann(w,n) bind(c)
! hann窗函数
    integer(C_INT), intent(in) :: n
    real(C_DOUBLE), intent(inout) :: w(n)

    integer(C_INT) i

    do i = 1, n, 1
        w(i) = 2.d0*w(i)*(0.5d0 - 0.5d0 * cos(TWO_PI * dble(i - 1) / dble(n - 1)))
    end do

end subroutine hann

subroutine welch(a,n,m,olr,psd,win) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: welch
! Welch 法估计功率谱
    integer(C_INT), intent(in) :: n, m
    real(C_DOUBLE), intent(in) :: a(n)
    real(C_DOUBLE), intent(out) :: psd(m)
    real(C_DOUBLE), intent(in) :: olr
    integer(C_INT), intent(in) :: win

    integer(C_INT) :: nol, i, j, k, q, nsub
    real(C_DOUBLE) :: re, im
    real(C_DOUBLE), allocatable :: c(:)
    complex(C_DOUBLE_COMPLEX), allocatable :: cf(:)

    type(c_ptr) :: plan

    if ( mod(m,2) == 0 ) then
        q = m/2 + 1
    else
        q = (m+1)/2
    end if

    nol = int(m*olr)
    k = m
    nsub = 1
    do while ( k < n )
        k = k - nol + m
        nsub = nsub + 1
    end do

    allocate(c(m))
    allocate(cf(m))

    psd = 0.d0
    c = a(1:m)
    plan = fftw_plan_dft_r2c_1d(m, c, cf, fftw_estimate)
    if (win>0) then 
        call hann(c,m)
    end if
    call fftw_execute_dft_r2c(plan, c, cf)

    do j = 1, q, 1
        re = dble(cf(j))
        im = aimag(cf(j))
        psd(j) = (re*re + im*im)/dble(m*nsub)
    end do

    k = m
    do i = 2, nsub-1, 1
        k = k - nol + m
        c = a((k-m+1):k)
        if (win>0) then 
            call hann(c,m)
        end if
        call fftw_execute_dft_r2c(plan, c, cf)
        do j = 1, q, 1
            re = dble(cf(j))
            im = aimag(cf(j))
            psd(j) = psd(j) + (re*re + im*im)/dble(m*nsub)
        end do
    end do

    k = k - nol + m
    c(1:(n-k+m)) = a((k-m+1):n)
    c((n-k+m+1):m) = 0.d0
    if (win>0) then 
        call hann(c,m)
    end if
    call fftw_execute_dft_r2c(plan, c, cf)
        
    do j = 1, q, 1
        re = dble(cf(j))
        im = aimag(cf(j))
        psd(j) = psd(j) + (re*re + im*im)/dble(m*nsub)
    end do

    call fftw_destroy_plan(plan)

    deallocate(c)
    deallocate(cf)

end subroutine welch

subroutine fftpadding(in,out,n,nfft) bind(c)
! 数据序列后补零
    integer, intent(in) :: n, nfft
    real*8, intent(in) :: in(n)
    real*8, intent(out) :: out(nfft)

    out(1:n) = in(1:n)
    out(n+1:nfft) = 0.d0

    return
end subroutine fftpadding

subroutine fftresample(a,n,r,ar,nr) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: fftresample
! 用FFT的方法进行数据重采样（降低采样频率）
    integer(C_INT), intent(in) :: n, r, nr
    real(C_DOUBLE), intent(in) :: a(n)
    real(C_DOUBLE), intent(out) :: ar(nr)

    integer(C_INT) :: Nfft, Nfftr
    real(C_DOUBLE), allocatable :: ap(:), apr(:)
    complex(C_DOUBLE_COMPLEX), allocatable :: af(:), afr(:)

    if ( r == 1 ) return

    Nfft = nextpow(r,n)

    allocate(ap(Nfft))
    allocate(af(Nfft))

    call fftpadding(a,ap,n,Nfft)
    call fft(ap,af,Nfft)

    Nfftr = Nfft/r
    allocate(apr(Nfft))
    allocate(afr(Nfft))

    if ( mod(Nfftr,2) == 0 ) then
        afr(1:Nfftr/2) = af(1:Nfftr/2)
        afr(Nfftr/2+1) = 1.0*af(Nfftr/2+1)
    else
        afr(1:(Nfftr-1)/2) = af(1:(Nfftr-1)/2)
        afr((Nfftr-1)/2+1) = 1.0*af((Nfftr-1)/2+1)
    end if

    call ifft(afr,apr,Nfftr)

    ar = (1.d0/dble(r))*apr(1:nr)

    deallocate(ap)
    deallocate(af)
    deallocate(apr)
    deallocate(afr)

    return
        
end subroutine fftresample

subroutine polyvals(p,m,x,y) bind(c)
! 多项式求值（单个值）
    integer, intent(in) :: m
    real*8, intent(in) :: p(m), x        
    real*8, intent(out) :: y

    integer :: i, n

    n = idnz(p)

    y = 0.D0
    do i = n, m, 1
        y = x*y + p(i)
    end do

    return

end subroutine polyvals

subroutine polyvalv(p,m,x,y,n) bind(c)
! 多项式求值（数组）
    integer, intent(in) :: m, n
    real*8, intent(in) :: p(m), x(n)
    real*8, intent(out) :: y(n)
    integer :: i, j, k

    k = idnz(p)

    !$OMP PARALLEL DEFAULT(NONE), PRIVATE(i,j), SHARED(p,m,x,y,n,k)
    !$OMP DO
    do i = 1, n, 1
        y(i) = 0.d0
        do j = k, m, 1
            y(i) = x(i)*y(i) + p(j)
        end do
    end do
    !$OMP END PARALLEL

    return
end subroutine polyvalv

subroutine polyder(p,m,n,r,l) bind(c)
! 多项式求导，p, r为求导前后的多项式系数（阶次从高到低，最后为常数项）
    integer, intent(in) :: m,n
    real*8, intent(in) :: p(m)
    real*8, intent(inout) :: r(l)
    integer, intent(inout) :: l

    integer :: i, j, k

    if ( l<0 ) then
        l = m - n
        if ( l<=0 ) l = 1
        return
    end if

    if ( n>=m ) then
        r(1) = 0.d0
        return
    end if

    do i = 1, l, 1
        r(i) = p(i)
        k = m-i
        do j = 1, n, 1
            r(i) = r(i)*k
            k = k-1
        end do
    end do
    return

end subroutine polyder

subroutine polyint(p,m,n,r,l) bind(c)
! 多项式积分，p, r为积分前后的多项式系数
    integer, intent(in) :: m,n
    real*8, intent(in) :: p(m)
    real*8, intent(inout) :: r(l)
    integer, intent(inout) :: l

    integer :: i, j, k

    if ( l<0 ) then
        l = m + n
        return
    end if

    r(m+1:l) = 0.d0

    do i = 1, m, 1
        r(i) = p(i)
        k = m-i
        do j = 1, n, 1
            r(i) = r(i)/dble(k+1)
            k = k+1
        end do
    end do
    return
end subroutine polyint

subroutine polymul(p,m,q,n,r,l) bind(c)
! 多项式乘法，r = p * q
    integer, intent(in) :: m,n
    real*8, intent(in) :: p(m),q(n)
    real*8, intent(inout) :: r(l)
    integer, intent(inout) :: l

    integer :: i, j

    if ( l<0 ) then
        l = m + n - 1
        return
    end if

    do i = 1, l, 1
        r(i) = 0.d0
        do j = 1, m, 1
            if (i-j > -1 .and. i-j < n) &
                r(i) = r(i) + p(i-j+1)*q(j)
        end do
    end do

    return

end subroutine polymul

subroutine polydiv(p,m,q,n,r,l) bind(c)
! 多项式除法，r = p / q
    integer, intent(in) :: m,n
    real*8, intent(in) :: p(m),q(n)
    integer, intent(inout) :: l
    real*8, intent(inout) :: r(l)

    integer :: i, j

    if ( l<0 ) then
        l = m + n - 1
        return
    end if

    do i = 1, l, 1
        r(i) = 0.d0
        do j = 1, m, 1
            if (i-j > -1 .and. i-j < n) &
                r(i) = r(i) + p(i-j+1)*q(j)
        end do
    end do

    return

end subroutine polydiv

subroutine polyadd(p,m,q,n,r,l) bind(c)
! 多项式加法，r = p + q
    integer, intent(in) :: m,n
    real*8, intent(in) :: p(m),q(n)
    integer, intent(inout) :: l
    real*8, intent(inout) :: r(l)

    integer :: i, j

    if ( l<0 ) then
        l = max(m,n)
        return
    end if

    if ( m>=n ) then
        do i = 1, l, 1
            if ( i<=m-n ) then
                r(i) = p(i)
            else
                r(i) = p(i) + q(i-m+n)
            end if
        end do
    else
        do i = 1, l, 1
            if ( i<=n-m ) then
                r(i) = q(i)
            else
                r(i) = q(i) + p(i-n+m)
            end if
        end do
    end if

    return

end subroutine polyadd

subroutine polyroots(p,m,r,l) bind(c)
! 多项式求根
! Find roots of the (m-1)-order polynomial by solving the eigen
! problem of the coresponding companion matrix.
! The roots are all in complex form.

    integer, intent(in) :: m
    real*8, intent(in) :: p(m)
    integer, intent(inout) :: l
    complex*16, intent(inout) :: r(l)

    real*8, allocatable :: A(:,:), work(:), wr(:), wi(:)
    real*8 :: temp(1)
    integer :: i, j, info, lwork

    if ( l<0 ) then
        l = m - 1
        if ( l<1 ) stop "Error: no roots; order too low (<1)!"
        return
    end if

    if ( m == 2 ) then
        r(1) = -p(2)/p(1)
        return
    end if

    allocate(A(l,l))
    allocate(wr(l))
    allocate(wi(l))

    A = 0.d0
    A(1,1) = -p(2)/p(1)
    do i = 2, l, 1
        A(i,i-1) = 1.d0
        A(1,i) = -p(i+1)/p(1)
    end do

    lwork = -1
    call dgeev ( "N", "N", l, A, l, wr, wi, temp, l, temp, l, temp, lwork, info )
    lwork = int(temp(1))
    allocate(work(lwork))
    call dgeev ( "N", "N", l, A, l, wr, wi, temp, l, temp, l, work, lwork, info )

    do i = 1, l, 1
        r(i) = cmplx(wr(i),wi(i))
    end do

    deallocate(A)
    deallocate(work)
    deallocate(wr)
    deallocate(wi)

    return

end subroutine polyroots

function idnz(p)
! 多项式系数数组中第一个非零值的位置
! Return the index of first non-zero element in the array
    real*8, intent(in) :: p(:)
    integer :: idnz

    idnz = 1
    do while ( abs(p(idnz)) < 1.0d-30 )
        idnz = idnz + 1
    end do
    return
end function idnz

subroutine polyfit(a,t,n,c,oh,ol) bind(c)
! 多项式拟合，a(n)为待拟合数据，t(n)为自变量
! oh, ol为拟合所用多项式的最高与最低阶次, c为拟合多项式的系数
    integer, intent(in) :: n
    real(8), intent(in) :: a(n), t(n)
    integer, intent(in) :: oh,ol
    real(8), intent(out) :: c(oh-ol+1)

    integer :: i, j, k
    real(8), allocatable :: b(:)
    real(8), allocatable :: M(:,:)

    k = oh-ol+1
    allocate(M(n,k))
    allocate(b(n))

    c = 0.D0
    b = a

    !$OMP PARALLEL DEFAULT(NONE), PRIVATE(i,j), SHARED(M,t,n,k,oh)
    !$OMP DO
    do i = 1, n, 1
        do j = 1, k, 1
            M(i,j) = t(i)**(oh-j+1)
        end do
    end do
    !$OMP END PARALLEL

    call leastsqs(M,b,n,k)

    c(1:k) = b(1:k)

    deallocate(M)
    deallocate(b)

    return

end subroutine polyfit

subroutine polydetrend(a,n,oh,ol) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: polydetrend
! 消除多项式趋势，oh, ol为拟合所用多项式的最高与最低阶次
    integer(4), intent(in) :: n
    real(8), intent(inout) :: a(n)
    integer(4), intent(in) :: oh,ol

    integer :: i
    real(8), allocatable :: t(:), c(:), ac(:)

    if ( oh == 0 ) then
        a = a - sum(a)/dble(n)
        return
    end if

    allocate(t(n))
    allocate(ac(n))
    allocate(c(oh+1))

    t = [ (((i-1)*1.D0),i=1,n) ]
    call polyfit(a,t,n,c,oh,ol)
    call polyval(c,oh+1,t,ac,n)

    a = a - ac

    deallocate(t)
    deallocate(c)
    deallocate(ac)

    return

end subroutine polydetrend

subroutine polyblc(a,n,oh,ol,dt,v0,d0) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: polyblc
! 改进的基线校正算法：
!       1. 先用多项式拟合速度时程，将拟合所得多项式求导得到第一次修正时程，
!       修正后重新积分计算速度、位移；
!       2. 然后用多项式拟合位移时程，将拟合所得多项式求导两次得到第二次修正时程；
    integer, intent(in) :: n
    real(8), intent(in) :: dt, v0, d0
    real(8), intent(inout) :: a(n)
    integer, intent(in) :: oh, ol

    integer :: i,l
    real(8), allocatable :: t(:), v(:), d(:), c(:), r(:), ac(:)

    allocate(t(n))
    allocate(v(n))
    allocate(d(n))
    allocate(ac(n))

    allocate(c(oh+2))
    l = oh+1
    allocate(r(l))

    t = [ (((i-1)*dt),i=1,n) ]

    call cumtrapz(a,v,n,dt,v0)
    call polyfit(v,t,n,c,oh+1,ol+1)
    call polyder(c,oh+2,1,r,l)
    call polyval(r,oh+1,t,ac,n)
    a = a - ac

    deallocate(c)
    allocate(c(oh+3))

    call acc2vd(a,v,d,n,dt,v0,d0)
    call polyfit(d,t,n,c,oh+2,ol+2)
    call polyder(c,oh+3,2,r,l)
    call polyval(r,oh+1,t,ac,n)

    a = a - ac

    deallocate(t)
    deallocate(ac)
    deallocate(v)
    deallocate(d)
    deallocate(c)
    deallocate(r)

    return
end subroutine polyblc

end module basic

module eqs

    use, intrinsic :: iso_c_binding
    use basic
    use omp_lib
    implicit none
    ! 自动选择频域与时域求解方法时，SDOF周期与采样间隔之比的临界值
    integer, parameter :: MPR = 20
    ! 非线性分析参数
    integer, parameter :: NPARA = 33
    real(8) :: para(NPARA) = 0.D0

    ! 本模块中一些函数的命名规则：
    ! 1. SDOF响应求解函数：rXXX 以r开头，当仅输出加速度响应时后跟a，XXX可能为freq, nmk, mixed
    ! 2. SDOF反应谱求解函数：spXXX 以r开头，当仅输出加速度谱时后跟a，XXX可能为freq, nmk, mixed

contains

subroutine spamixed(acc,n,dt,zeta,P,nP,SPA,SPI) bind(c)
    integer, intent(in) :: n, nP
    real(C_DOUBLE), intent(in) :: dt, zeta
    real(C_DOUBLE), intent(in) :: acc(n), P(nP)
    real(C_DOUBLE), intent(out) :: SPA(nP)
    integer, intent(out) :: SPI(nP)

    integer :: m

    m = 1
    do while ( P(m)<=MPR*dt )
        m = m + 1
    end do

    call spafreq(acc,n,dt,zeta,P(1:m),m,SPA(1:m),SPI(1:m))
    if ( m<nP ) then
        call spanmk(acc,n,dt,zeta,P(m+1:nP),(nP-m),SPA(m+1:nP),SPI(m+1:nP))
    end if
end subroutine spamixed

subroutine spamixed_endur(acc,n,dt,zeta,P,nP,DI,nD,SPA,SPI) bind(c)

    integer(C_INT), intent(in) :: n, nP, nD
    integer(C_INT), intent(in) :: DI(nD)
    real(C_DOUBLE), intent(in) :: dt, zeta
    real(C_DOUBLE), intent(in) :: acc(n), P(nP)
    real(C_DOUBLE), intent(out) :: SPA(nP,nd)
    integer(C_INT), intent(out) :: SPI(nP,nd)

    integer(C_INT) :: k,i,ind
    real(C_DOUBLE) :: SPVk, SPDk, SPEk
    real(C_DOUBLE), allocatable :: ra(:), rv(:), rd(:), ara(:)

    if ( nd < 2 ) return

    allocate(ra(n))
    allocate(ara(n))

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,i,ra,ara)
    do k = 1, nP, 1
        call ramixed(acc,n,dt,zeta,P(k),ra)
        ara = abs(ra)
        i = 1
        ind = maxloc(ara(1:DI(i)),1)
        SPI(k,i) = ind
        SPA(k,i) = ra(ind)
        do i = 2, nD, 1
            ind = maxloc(ara(DI(i-1)+1:DI(i)),1)+DI(i-1)
            SPI(k,i) = ind
            SPA(k,i) = ra(maxloc(ara(1:DI(i)),1))
        end do
    end do
    !$OMP END PARALLEL DO

    deallocate(ra)
    deallocate(ara)

end subroutine spamixed_endur

subroutine pspamixed(acc,n,dt,zeta,P,nP,SPA,SPI) bind(c)
    integer, intent(in) :: n, nP
    real(C_DOUBLE), intent(in) :: dt, zeta
    real(C_DOUBLE), intent(in) :: acc(n), P(nP)
    real(C_DOUBLE), intent(out) :: SPA(nP)
    integer, intent(out) :: SPI(nP)

    integer :: m

    m = 1

    do while ( P(m)<=MPR*dt )
        m = m + 1
    end do

    call pspafreq(acc,n,dt,zeta,P(1:m),m,SPA(1:m),SPI(1:m))
    if ( m<nP ) then
        call pspanmk(acc,n,dt,zeta,P(m+1:nP),(nP-m),SPA(m+1:nP),SPI(m+1:nP))
    end if
end subroutine pspamixed

subroutine spavdmixed(acc,n,dt,zeta,P,nP,SPA,SPI,SPV,SPD,SPE) bind(c)
    integer, intent(in) :: n, nP
    real(C_DOUBLE), intent(in) :: dt, zeta
    real(C_DOUBLE), intent(in) :: acc(n), P(nP)
    real(C_DOUBLE), intent(out) :: SPA(nP), SPV(nP), SPD(nP), SPE(nP)
    integer, intent(out) :: SPI(nP)

    integer :: m
    m = 1
    do while ( P(m)<=MPR*dt )
        m = m + 1
    end do

    call spavdfreq(acc,n,dt,zeta,P(1:m),m,SPA(1:m),SPI(1:m),&
        SPV(1:m),SPD(1:m),SPE(1:m))
    if ( m<nP ) then
        call spavdnmk(acc,n,dt,zeta,P(m+1:nP),(nP-m),SPA(m+1:nP),&
            SPI(m+1:nP),SPV(m+1:nP),SPD(m+1:nP),SPE(m+1:nP))
    end if
end subroutine spavdmixed

subroutine spafreq(acc,n,dt,zeta,P,nP,SPA,SPI) bind(c)

    integer, intent(in) :: n, nP
    real(C_DOUBLE), intent(in) :: dt, zeta
    real(C_DOUBLE), intent(in) :: acc(n), P(nP)
    real(C_DOUBLE), intent(out) :: SPA(nP)
    integer, intent(out) :: SPI(nP)

    complex(C_DOUBLE_COMPLEX), allocatable :: af(:), raf(:)
    real(C_DOUBLE), allocatable :: a0(:), ra(:), w(:), w0(:)
    integer :: i, j, Nfft
    complex(C_DOUBLE_COMPLEX) :: IONE=cmplx(0.d0,1.d0)
    real(C_DOUBLE) :: w0i2, w0iwj, wj2
    type(C_PTR) :: plan, iplan

    SPA = 1.D0
    SPI = 1

    Nfft = nextpow2(n)

    allocate(a0(Nfft))
    allocate(af(Nfft))
    allocate(raf(Nfft))
    allocate(ra(Nfft))
    allocate(w(Nfft))
    allocate(w0(nP))

    call fftfreqs(Nfft,1.d0/dt,w)
    w = w*TWO_PI
    w0 = TWO_PI/P

    ra = 0.d0
    a0 = 0.d0
    a0(1:n) = acc

    plan  = fftw_plan_dft_r2c_1d(Nfft, a0, af, FFTW_ESTIMATE)
    iplan = fftw_plan_dft_c2r_1d(Nfft, af, ra, FFTW_ESTIMATE)
    call fftw_execute_dft_r2c(plan, a0, af)

    !$OMP PARALLEL DO PRIVATE(i, j, w0i2, w0iwj, wj2, raf, ra), SHARED(nP,Nfft,IONE,iplan,SPI,SPA)

    do i = 1, nP, 1
        w0i2 = w0(i)*w0(i)
        do j = 1, Nfft/2+1, 1
            w0iwj = w0(i)*w(j)
            wj2 = w(j)*w(j)
            raf(j) = af(j)*(w0i2+2.d0*zeta*w0iwj*IONE)/&
                    (w0i2-wj2+2.d0*zeta*w0iwj*IONE)
        end do
        call fftw_execute_dft_c2r(iplan, raf, ra)

        SPI(i) = maxloc(dabs(ra(1:n)),1)
        SPA(i) = ra(SPI(i))/dble(Nfft)
    end do

    !$OMP END PARALLEL DO

    call fftw_destroy_plan(plan)
    call fftw_destroy_plan(iplan)

    deallocate( a0)
    deallocate( af)
    deallocate(raf)
    deallocate( ra)
    deallocate(  w)
    deallocate( w0)

    return
end subroutine spafreq

subroutine spafreq_endur(acc,n,dt,zeta,P,nP,DI,nD,SPA,SPI) bind(c)

    integer, intent(in) :: n, nP, nD
    integer, intent(in) :: DI(nD)
    real(8), intent(in) :: dt, zeta
    real(8), intent(in) :: acc(n), P(nP)
    real(8), intent(out) :: SPA(nP,nd)
    integer, intent(out) :: SPI(nP,nd)

    integer :: k,i,ind
    real(8) :: SPVk, SPDk, SPEk
    real(8), allocatable :: ra(:), ara(:)

    if ( nd < 2 ) return

    allocate(ra(n))
    allocate(ara(n))

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,i,ra,ara)
    do k = 1, nP, 1
        call rafreq(acc,n,dt,zeta,P(k),ra)
        ara = abs(ra)
        i = 1
        ind = maxloc(ara(1:DI(i)),1)
        SPI(k,i) = ind
        SPA(k,i) = ra(ind)
        do i = 2, nD, 1
            ind = maxloc(ara(DI(i-1)+1:DI(i)),1)+DI(i-1)
            SPI(k,i) = ind
            SPA(k,i) = ra(maxloc(ara(1:DI(i)),1))
        end do
    end do
    !$OMP END PARALLEL DO

end subroutine spafreq_endur

subroutine pspafreq(acc,n,dt,zeta,P,nP,SPA,SPI) bind(c)

    integer, intent(in) :: n, nP
    real(C_DOUBLE), intent(in) :: dt, zeta
    real(C_DOUBLE), intent(in) :: acc(n), P(nP)
    real(C_DOUBLE), intent(out) :: SPA(nP)
    integer, intent(out) :: SPI(nP)

    complex(C_DOUBLE_COMPLEX), allocatable :: af(:), raf(:)
    real(C_DOUBLE), allocatable :: a0(:), ra(:), w(:), w0(:)
    integer :: i, j, Nfft
    complex(C_DOUBLE_COMPLEX) :: IONE=cmplx(0.d0,1.d0)
    real(C_DOUBLE) :: w0i2, w0iwj, wj2
    type(C_PTR) :: plan, iplan

    SPA = 1.D0
    SPI = 1

    Nfft = nextpow2(n)

    allocate(a0(Nfft))
    allocate(af(Nfft))
    allocate(raf(Nfft))
    allocate(ra(Nfft))
    allocate(w(Nfft))
    allocate(w0(nP))

    call fftfreqs(Nfft,1.d0/dt,w)
    w = w*TWO_PI
    w0 = TWO_PI/P

    ra = 0.d0
    a0 = 0.d0
    a0(1:n) = acc

    plan  = fftw_plan_dft_r2c_1d(Nfft, a0, af, FFTW_ESTIMATE)
    iplan = fftw_plan_dft_c2r_1d(Nfft, af, ra, FFTW_ESTIMATE)
    call fftw_execute_dft_r2c(plan, a0, af)

    do i = 1, nP, 1
        w0i2 = w0(i)*w0(i)
        do j = 1, Nfft/2+1, 1
            w0iwj = w0(i)*w(j)
            wj2 = w(j)*w(j)
            raf(j) = af(j)*(-1.d0/(w0i2-wj2+2.d0*zeta*w0iwj*IONE))
        end do
        call fftw_execute_dft_c2r(iplan, raf, ra)

        SPI(i) = maxloc(dabs(ra(1:n)),1)
        SPA(i) = ra(SPI(i))/dble(Nfft)*39.47841760435743D0/P(i)/P(i)
    end do

    call fftw_destroy_plan(plan)
    call fftw_destroy_plan(iplan)

    deallocate( a0)
    deallocate( af)
    deallocate(raf)
    deallocate( ra)
    deallocate(  w)
    deallocate( w0)

    return
end subroutine pspafreq

subroutine spavdfreq(acc,n,dt,zeta,P,nP,SPA,SPI,SPV,SPD,SPE) bind(c)

    integer, intent(in) :: n, nP
    real(C_DOUBLE), intent(in) :: dt, zeta
    real(C_DOUBLE), intent(in) :: acc(n), P(nP)
    real(C_DOUBLE), intent(out) :: SPA(nP), SPV(nP), SPD(nP), SPE(nP)
    integer, intent(out) :: SPI(nP)

    complex(C_DOUBLE_COMPLEX), allocatable :: af(:), raf(:), rvf(:), rdf(:)
    real(C_DOUBLE), allocatable :: a0(:), ra(:), w(:), w0(:)
    real(C_DOUBLE), allocatable :: rv(:), rd(:)
    integer :: i, j, Nfft
    complex(C_DOUBLE_COMPLEX) :: IONE=cmplx(0.d0,1.d0)
    real(C_DOUBLE) :: w0i2, w0iwj, wj2, iNfft
    type(C_PTR) :: plan, iplan

    SPA = 1.D0
    SPV = 1.D0
    SPD = 1.D0
    SPI = 1

    Nfft = nextpow2(n)
    iNfft = 1.d0/dble(Nfft)

    allocate(a0(Nfft))
    allocate(af(Nfft))
    allocate(raf(Nfft))
    allocate(rvf(Nfft))
    allocate(rdf(Nfft))
    allocate(ra(Nfft))
    allocate(rv(Nfft))
    allocate(rd(Nfft))
    allocate(w(Nfft))
    allocate(w0(nP))

    call fftfreqs(Nfft,1.d0/dt,w)
    w = w*TWO_PI
    w0 = TWO_PI/P

    ra = 0.d0
    rv = 0.d0
    rd = 0.d0
    a0 = 0.d0
    a0(1:n) = acc

    plan  = fftw_plan_dft_r2c_1d(Nfft, a0, af, FFTW_ESTIMATE)
    iplan = fftw_plan_dft_c2r_1d(Nfft, af, ra, FFTW_ESTIMATE)
    call fftw_execute_dft_r2c(plan, a0, af)

    do i = 1, nP, 1
        w0i2 = w0(i)*w0(i)
        do j = 1, Nfft/2+1, 1
            w0iwj = w0(i)*w(j)
            wj2 = w(j)*w(j)
            raf(j) = af(j)*(w0i2+2.D0*zeta*w0iwj*IONE)/&
                    (w0i2-wj2+2.d0*zeta*w0iwj*IONE)
            rvf(j) = af(j)*(-w(j)*IONE/(w0i2-wj2+2.d0*zeta*w0iwj*IONE))
            rdf(j) = af(j)*(-1.d0/(w0i2-wj2+2.d0*zeta*w0iwj*IONE))
        end do
        call fftw_execute_dft_c2r(iplan, raf, ra)
        call fftw_execute_dft_c2r(iplan, rvf, rv)
        call fftw_execute_dft_c2r(iplan, rdf, rd)
        SPI(i) = maxloc(dabs(ra(1:n)),1)
        SPA(i) = dabs(ra(SPI(i)))*iNfft
        SPV(i) = maxval(dabs(rv(1:n)))*iNfft
        SPD(i) = maxval(dabs(rd(1:n)))*iNfft
        SPE(i) = sqrt(2.d0*trapz(-rv(1:n)*acc*iNfft,n,dt,0.d0))
    end do

    call fftw_destroy_plan(plan)
    call fftw_destroy_plan(iplan)

    deallocate( a0)
    deallocate( af)
    deallocate(raf)
    deallocate(rvf)
    deallocate(rdf)
    deallocate( ra)
    deallocate( rv)
    deallocate( rd)
    deallocate(  w)
    deallocate( w0)

    return
end subroutine spavdfreq

subroutine spanmk(acc,n,dt,zeta,P,nP,SPA,SPI) bind(c)

    integer, intent(in) :: n, nP
    real(C_DOUBLE), intent(in) :: dt, zeta
    real(C_DOUBLE), intent(in) :: acc(n), P(nP)
    real(C_DOUBLE), intent(out) :: SPA(nP)
    integer, intent(out) :: SPI(nP)

    integer :: k
    real(8) :: SPVk, SPDk, SPEk

    !$OMP PARALLEL DO PRIVATE(k)
    do k = 1, nP, 1
        call newmark(acc,n,dt,zeta,P(k),SPA(k),SPI(k),SPVk,SPDk,SPEk)
    end do
    !$OMP END PARALLEL DO

end subroutine spanmk

subroutine spanmk_endur(acc,n,dt,zeta,P,nP,DI,nD,SPA,SPI) bind(c)

    integer, intent(in) :: n, nP, nD
    integer, intent(in) :: DI(nD)
    real(8), intent(in) :: dt, zeta
    real(8), intent(in) :: acc(n), P(nP)
    real(8), intent(out) :: SPA(nP,nd)
    integer, intent(out) :: SPI(nP,nd)

    integer :: k,i,ind
    real(8) :: SPVk, SPDk, SPEk
    real(8), allocatable :: ra(:), ara(:)

    if ( nd < 2 ) return

    allocate(ra(n))
    allocate(ara(n))

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,i,ra,ara)
    do k = 1, nP, 1
        call ranmk(acc,n,dt,zeta,P(k),ra)
        ara = abs(ra)
        i = 1
        ind = maxloc(ara(1:DI(i)),1)
        SPI(k,i) = ind
        SPA(k,i) = ra(ind)
        do i = 2, nD, 1
            ind = maxloc(ara(DI(i-1)+1:DI(i)),1)+DI(i-1)
            SPI(k,i) = ind
            SPA(k,i) = ra(maxloc(ara(1:DI(i)),1))
        end do
    end do
    !$OMP END PARALLEL DO

end subroutine spanmk_endur

subroutine pspanmk(acc,n,dt,zeta,P,nP,SPA,SPI) bind(c)

    integer, intent(in) :: n, nP
    real(C_DOUBLE), intent(in) :: dt, zeta
    real(C_DOUBLE), intent(in) :: acc(n), P(nP)
    real(C_DOUBLE), intent(out) :: SPA(nP)
    integer, intent(out) :: SPI(nP)

    integer :: k
    real(8) :: SPVk, SPDk, SPEk

    do k = 1, nP, 1
        call newmark(acc,n,dt,zeta,P(k),SPDk,SPI(k),SPVk,SPA(k),SPEk)
        SPA(k) = SPA(k)*39.47841760435743D0/P(k)/P(k)
    end do
end subroutine pspanmk

subroutine spavdnmk(acc,n,dt,zeta,P,nP,SPA,SPI,SPV,SPD,SPE) bind(c)

    integer, intent(in) :: n, nP
    real(C_DOUBLE), intent(in) :: dt, zeta
    real(C_DOUBLE), intent(in) :: acc(n), P(nP)
    real(C_DOUBLE), intent(out) :: SPA(nP), SPV(nP), SPD(nP), SPE(nP)
    integer, intent(out) :: SPI(nP)

    integer :: k

    do k = 1, nP, 1
        call newmark(acc,n,dt,zeta,P(k),SPA(k),SPI(k),SPV(k),SPD(k),SPE(k))
    end do
end subroutine spavdnmk

subroutine newmark(acc,n,dt0,zeta,P,SPA,SPI,SPV,SPD,SPE) bind(c)

    integer, intent(in) :: n
    real(8), intent(in) :: dt0, zeta
    real(8), intent(in) :: acc(n), P
    real(8), intent(out) :: SPA, SPV, SPD, SPE
    integer, intent(out) :: SPI

    real(8) :: rc(3),rl(3),ac,al,aa,ein,da
    real(8) :: beta = 1.d0/4.d0, gamma = 0.5d0
    real(8) :: k,c,w,keff,kinv,feff,dt,rvl
    real(8) :: b1, b2, b3, b4, b5, b6, b7, b8
    integer :: i,j,r

    if ( dt0*MPR > P ) then
        r  = ceiling(MPR*dt0/P)
        dt = dt0/r
    else
        r = 1
        dt = dt0
    end if

    b1 = 1.d0/(beta*dt*dt)
    b2 = 1.d0/(beta*dt)
    b3 = 1.d0/(beta*2.0d0)-1.d0
    b4 = gamma/(beta*dt)
    b5 = gamma/beta-1.d0
    b6 = 0.5d0*dt*(gamma/beta-2.0d0)
    b7 = dt*(1.0d0-gamma)
    b8 = dt*gamma

    w = TWO_PI/P
    k = w**2.d0
    c = 2.0d0*zeta*w
    keff = k+b1+b4*c
    kinv = 1.0d0/keff

    !  initiation
    rl = 0.d0
    rc = 0.d0
    SPA = 0.0d0
    SPV = 0.0d0
    SPD = 0.0d0
    SPI = 0
    ein = 0.0d0
    al = 0.d0
    rvl = 0.d0

    !  computation of response
    do i = 1,n,1

        ac = acc(i)
        da = (ac - al)/dble(r)

        do j = 1, r, 1
            ac = al + da*j
            feff = ac+(b1*rl(1)+b2*rl(2)+b3*rl(3))&
            +c*(b4*rl(1)+b5*rl(2)+b6*rl(3))
            rc(1) = feff*kinv
            rc(2) = b4*(rc(1)-rl(1))-b5*rl(2)-b6*rl(3)
            rc(3) = ac-k*rc(1)-c*rc(2)
            rl(1) = rc(1)
            rl(2) = rc(2)
            rl(3) = rc(3)
        end do

        ein = ein+0.5d0*dt0*(rvl*al+rc(2)*ac)
        aa = -rc(3)+ac

        if (abs(aa) > abs(SPA)) then
            SPA=aa
            SPI=i
        endif

        if (abs(rc(2)) > abs(SPV)) then
            SPV=rc(2)
        endif

        if (abs(rc(1)) > abs(SPD)) then
            SPD=rc(1)
        endif

        al = ac
        rvl = rc(2)

    end do
    SPE = sqrt(2.d0*ein)
    return
end subroutine newmark

subroutine rnmk(acc,n,dt0,zeta,P,ra,rv,rd) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: rnmk
    integer, intent(in) :: n
    real(8), intent(in) :: dt0, zeta
    real(8), intent(in) :: acc(n), P
    real(8), intent(out) :: ra(n),rv(n),rd(n)

    real(8) :: rc(3),rl(3),ac,al,aa,da
    real(8) :: beta = 1.d0/4.d0, gamma = 0.5d0
    real(8) :: k,c,w,keff,kinv,feff,dt,rvl
    real(8) :: b1, b2, b3, b4, b5, b6, b7, b8
    integer :: i,j,r

    if ( dt0*MPR > P ) then
        r  = ceiling(MPR*dt0/P)
        dt = dt0/r
    else
        r = 1
        dt = dt0
    end if

    b1 = 1.d0/(beta*dt*dt)
    b2 = 1.d0/(beta*dt)
    b3 = 1.d0/(beta*2.0d0)-1.d0
    b4 = gamma/(beta*dt)
    b5 = gamma/beta-1.d0
    b6 = 0.5d0*dt*(gamma/beta-2.0d0)
    b7 = dt*(1.0d0-gamma)
    b8 = dt*gamma

    w = TWO_PI/P
    k = w**2.d0
    c = 2.0d0*zeta*w
    keff = k+b1+b4*c
    kinv = 1.0d0/keff

    !  initiation
    rl = 0.d0
    rc = 0.d0
    ra = 0.d0
    rv = 0.d0
    rd = 0.d0
    al = 0.d0

    !  computation of response
    do i = 1, n, 1

        ac = acc(i)
        da = (ac - al)/dble(r)

        do j = 1, r, 1
            ac = al + da*j
            feff = ac+(b1*rl(1)+b2*rl(2)+b3*rl(3))&
            +c*(b4*rl(1)+b5*rl(2)+b6*rl(3))
            rc(1) = feff*kinv
            rc(2) = b4*(rc(1)-rl(1))-b5*rl(2)-b6*rl(3)
            rc(3) = ac-k*rc(1)-c*rc(2)
            rl(1) = rc(1)
            rl(2) = rc(2)
            rl(3) = rc(3)
        end do

        ra(i) = -rc(3)+ac
        rv(i) = -rc(2)
        rd(i) = -rc(1)

        al = ac

    end do
    return
end subroutine rnmk

subroutine ranmk(acc,n,dt0,zeta,P,ra) bind(c)

    integer, intent(in) :: n
    real(8), intent(in) :: dt0, zeta
    real(8), intent(in) :: acc(n), P
    real(8), intent(out) :: ra(n)

    real(8), allocatable :: rv(:), rd(:)

    allocate(rv(n))
    allocate(rd(n))

    call rnmk(acc,n,dt0,zeta,P,ra,rv,rd)

    deallocate(rv)
    deallocate(rd)

    return
end subroutine ranmk

subroutine energy(acc,n,dt,zeta,P,ra,rv,rd,Ek,Ez,Es,Eh,Ein,ku) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: energy
    integer, intent(in) :: n
    real(8), intent(in) :: acc(n), ra(n), rv(n), rd(n), ku(n)
    real(8), intent(in) :: dt, zeta, P
    real(8), intent(out) :: Ek(n),Ez(n),Es(n),Eh(n),Ein(n)

    integer :: i
    real(8) :: w, fs(2), Ea, c, ddsub(2)

    w = TWO_PI/P
    c = 2.d0*zeta*w

    Ek(1) = rv(1)*rv(1)*0.5D0
    Ez(1) = 0.5D0*dt*(0.d0+c*rv(1)*rv(1))
    fs = [0.D0, -ra(1)-c*rv(1)]
    Ea = 0.5D0*dt*(0.d0+fs(2)*rv(1))
    Es(1) = fs(2)*fs(2)*0.5D0/(ku(1))
    Eh(1) = Ea - Es(1)
    Ein(1) = 0.5D0*dt*(0.d0-acc(1)*rv(1))
!     BL(1) = Ein(1) - Ek(1) - Ez(1) - Ea

    do i = 2, n, 1
        Ek(i) = rv(i)*rv(i)*0.5D0
        Ez(i) = Ez(i-1) + 0.5D0*dt*c*(rv(i-1)*rv(i-1)+rv(i)*rv(i))

        fs(1) = - ra(i-1) - c*rv(i-1)
        fs(2) = - ra(i) - c*rv(i)
        
        if (fs(1)*fs(2)<0.D0) then
            ddsub(1) = (0.D0-fs(1))/ku(i-1)
            ddsub(2) = fs(2)/ku(i)
            Ea = Ea + 0.5D0*ddsub(1)*(fs(1)+0.D0) + 0.5D0*ddsub(2)*(0.D0+fs(2))
        else
            Ea = Ea + 0.5D0*(rd(i)-rd(i-1))*(fs(1)+fs(2))
        end if
        
        !Ea = Ea + 0.5D0*(rd(i)-rd(i-1))*(fs(1)+fs(2))
        !Ea = Ea + 0.5D0*dt*(fs(1)*rv(i-1)+fs(2)*rv(i))
        Es(i) = fs(2)*fs(2)*0.5D0/(ku(i))
        Eh(i) = Ea - Es(i)

        Ein(i) = Ein(i-1) - 0.5D0*dt*(acc(i-1)*rv(i-1)+acc(i)*rv(i))
!         BL(i) = Ein(i) - Ek(i) - Ez(i) - Ea

    end do
end subroutine energy

subroutine rfreq(acc,n,dt,zeta,P,ra,rv,rd) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: rfreq
    integer, intent(in) :: n
    real(C_DOUBLE), intent(in) :: dt, zeta, P
    real(C_DOUBLE), intent(in) :: acc(n)
    real(C_DOUBLE), intent(out) :: ra(n),rv(n),rd(n)

    real(C_DOUBLE), allocatable :: rac(:), rvc(:), rdc(:)
    complex(C_DOUBLE_COMPLEX), allocatable :: af(:), raf(:), rvf(:), rdf(:)
    real(C_DOUBLE), allocatable :: a0(:), w(:)
    integer :: i, j, Nfft
    complex(C_DOUBLE_COMPLEX) :: IONE=cmplx(0.d0,1.d0)
    real(C_DOUBLE) :: w0i2, w0iwj, wj2, iNfft, w0
    type(C_PTR) :: plan, iplan

    Nfft = nextpow2(n)*2
    iNfft = 1.d0/dble(Nfft)

    allocate(a0(Nfft))
    allocate(af(Nfft))
    allocate(rac(Nfft))
    allocate(rvc(Nfft))
    allocate(rdc(Nfft))
    allocate(raf(Nfft))
    allocate(rvf(Nfft))
    allocate(rdf(Nfft))
    allocate(w(Nfft))

    call fftfreqs(Nfft,1.d0/dt,w)
    w = w*TWO_PI
    w0 = TWO_PI/P

    a0 = 0.d0
    a0(1:n) = acc
    
    plan  = fftw_plan_dft_r2c_1d(Nfft, a0, af, FFTW_ESTIMATE)
    iplan = fftw_plan_dft_c2r_1d(Nfft, af, rac, FFTW_ESTIMATE)
    call fftw_execute_dft_r2c(plan, a0, af)

    w0i2 = w0*w0
    do j = 1, Nfft/2+1, 1
        w0iwj = w0*w(j)
        wj2 = w(j)*w(j)
        raf(j) = af(j)*(w0i2+2.D0*zeta*w0iwj*IONE)/&
                (w0i2-wj2+2.d0*zeta*w0iwj*IONE)*iNfft
        rvf(j) = af(j)*(-w(j)*IONE/(w0i2-wj2+2.d0*zeta*w0iwj*IONE))*iNfft
        rdf(j) = af(j)*(-1.d0/(w0i2-wj2+2.d0*zeta*w0iwj*IONE))*iNfft
    end do
    
    call fftw_execute_dft_c2r(iplan, raf, rac)
    call fftw_execute_dft_c2r(iplan, rvf, rvc)
    call fftw_execute_dft_c2r(iplan, rdf, rdc)

    do i = 1, n, 1
        ra(i) = rac(i)
        rv(i) = rvc(i)
        rd(i) = rdc(i)
    end do

    call fftw_destroy_plan(plan)
    call fftw_destroy_plan(iplan)

    deallocate( a0)
    deallocate( af)
    deallocate(raf)
    deallocate(rvf)
    deallocate(rdf)
    deallocate(rac)
    deallocate(rvc)
    deallocate(rdc)
    deallocate(  w)

    return
end subroutine rfreq

subroutine rafreq(acc,n,dt,zeta,P,ra) bind(c)

    integer, intent(in) :: n
    real(C_DOUBLE), intent(in) :: dt, zeta, P
    real(C_DOUBLE), intent(in) :: acc(n)
    real(C_DOUBLE), intent(out) :: ra(n)

    complex(C_DOUBLE_COMPLEX), allocatable :: af(:), raf(:)
    real(C_DOUBLE), allocatable :: a0(:), w(:), rac(:)
    integer :: i, j, Nfft
    complex(C_DOUBLE_COMPLEX) :: IONE=cmplx(0.d0,1.d0)
    real(C_DOUBLE) :: w0i2, w0iwj, wj2, iNfft, w0
    type(C_PTR) :: plan, iplan

    Nfft = nextpow2(n)*2
    iNfft = 1.d0/dble(Nfft)
    
    allocate(rac(Nfft))
    allocate(a0(Nfft))
    allocate(af(Nfft))
    allocate(raf(Nfft))
    allocate(w(Nfft))
    
    call fftfreqs(Nfft,1.d0/dt,w)
    w = w*TWO_PI
    w0 = TWO_PI/P
    
    rac = 0.d0
    a0 = 0.d0
    a0(1:n) = acc
    
    plan  = fftw_plan_dft_r2c_1d(Nfft, a0, af, FFTW_ESTIMATE)
    iplan = fftw_plan_dft_c2r_1d(Nfft, af, rac, FFTW_ESTIMATE)
    call fftw_execute_dft_r2c(plan, a0, af)
    
    w0i2 = w0*w0
    do j = 1, Nfft/2+1, 1
        w0iwj = w0*w(j)
        wj2 = w(j)*w(j)
        raf(j) = af(j)*(w0i2+2.D0*zeta*w0iwj*IONE)/&
                (w0i2-wj2+2.d0*zeta*w0iwj*IONE)*iNfft
    end do
    call fftw_execute_dft_c2r(iplan, raf, rac)
    ra = rac(1:n)
    
    call fftw_destroy_plan(plan)
    call fftw_destroy_plan(iplan)
    
    deallocate(rac)
    deallocate( a0)
    deallocate( af)
    deallocate(raf)
    deallocate(  w)

    return
    
end subroutine rafreq

subroutine rmixed(acc,n,dt,zeta,P,ra,rv,rd) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: rmixed
    integer, intent(in) :: n
    real(C_DOUBLE), intent(in) :: dt, zeta, P
    real(C_DOUBLE), intent(in) :: acc(n)
    real(C_DOUBLE), intent(out) :: ra(n),rv(n),rd(n)

    if (P<=MPR*dt) then
        call rfreq(acc,n,dt,zeta,P,ra,rv,rd)
    else
        call rnmk(acc,n,dt,zeta,P,ra,rv,rd)
    end if
end subroutine rmixed

subroutine ramixed(acc,n,dt,zeta,P,ra) bind(c)

    integer, intent(in) :: n
    real(C_DOUBLE), intent(in) :: dt, zeta, P
    real(C_DOUBLE), intent(in) :: acc(n)
    real(C_DOUBLE), intent(out) :: ra(n)

    if (P<=MPR*dt) then
        call rafreq(acc,n,dt,zeta,P,ra)
    else
        call ranmk(acc,n,dt,zeta,P,ra)
    end if
end subroutine ramixed

subroutine fitspectra(acc,n,dt,zeta,P,nP,SPAT,a,tol,mit,kpb) bind(c)

    integer, intent(in) :: n, nP, mit, kpb
    real(C_DOUBLE), intent(in) :: dt, zeta, tol
    real(C_DOUBLE), intent(in) :: acc(n), P(nP), SPAT(nP)
    real(C_DOUBLE), intent(out) :: a(n)

    real(C_DOUBLE), allocatable :: SPA(:), R(:), Rf(:), Pf(:), best(:)
    integer(C_INT), allocatable :: SPI(:)

    real(C_DOUBLE) :: aerror, merror, peak0, peak, minerr
    complex(C_DOUBLE_COMPLEX), allocatable :: af(:), afs(:)
    real(C_DOUBLE), allocatable :: a0(:), f(:)
    integer :: i, j, Nfft, NPf, IPf1, IPf2, iter, pk
    complex(C_DOUBLE_COMPLEX) :: IONE=cmplx(0.d0,1.d0)
    real(C_DOUBLE) :: pfi, iNfft, p0, pe0, dp, phi, t, ra
    type(C_PTR) :: plan, iplan

    peak0 = maxval(abs(acc), dim=1)
    Nfft = nextpow2(n)*4
    iNfft = 1.d0/dble(Nfft)

    allocate(best(n))
    allocate(a0(Nfft))
    allocate(af(Nfft))
    allocate(afs(Nfft))
    allocate(f(Nfft))
    allocate(Pf(Nfft/2))
    allocate(Rf(Nfft/2))
    a = acc*1.D0

    call fftfreqs(Nfft,1.d0/dt,f)
    Pf(2:Nfft/2) = 1.d0/f(2:Nfft/2)
    Pf(1) = 2.d0*Pf(2)

    call decrfindfirst(Pf,Nfft/2,P(nP),IPf1)
    call decrfindlast(Pf,Nfft/2,P(1), IPf2)

    NPf = IPf2-IPf1+1

    a0 = 0.d0
    a0(1:n) = acc

    plan  = fftw_plan_dft_r2c_1d(Nfft, a0, af, FFTW_ESTIMATE)
    iplan = fftw_plan_dft_c2r_1d(Nfft, af, a0, FFTW_ESTIMATE)
    call fftw_execute_dft_r2c(plan, a0, af)

    allocate(SPA(nP))
    allocate(R(nP))
    allocate(SPI(nP))

    call spamixed(acc,n,dt,zeta,P,nP,SPA,SPI)
    call error(abs(SPA),SPAT,nP,aerror,merror)
    if ( aerror<=tol .and. merror<=3.0d0*tol ) return
    R = SPAT/abs(SPA)
    call decrlininterp(P,R,nP,Pf(IPf1:IPf2),Rf(IPf1:IPf2),NPf)
    ! write(unit=*, fmt="(A29,2F8.4)") "Initial Error: ",aerror,merror
    minerr = merror
    best = a*1.D0
    iter = 1
    do while ( (aerror>tol .or. merror>3.0d0*tol) .and. iter<=mit )

        j = IPf2
        do i = 2, nP-1, 1
            p0 = p(i)
            t = dt*SPI(i)-dt
            dp = min(p0-P(i-1),P(i+1)-p0)

            do while ( Pf(j)<p0+0.5d0*dp .and. j>IPf1)
                pe0 = Pf(j)
                phi = atan2(aimag(af(j)),real(af(j)))
                ra = rsimple(pe0,p0,phi,t,zeta)
                if ( sign(1.d0,ra)*sign(1.d0,SPA(i)) < 0.d0) then
                    Rf(j) = 1.d0/Rf(j)
                end if
                j = j - 1
            end do
        end do

        af(IPf1:IPf2) = af(IPf1:IPf2)*Rf(IPf1:IPf2)
        afs = af
        call fftw_execute_dft_c2r(iplan, afs, a0)
        a = a0(1:n)*iNfft

        call adjustpeak(a,n,peak0)
        ! call adjustbaseline(a,n,dt)
        ! call adjustpeak(a,n,peak0)

        call spamixed(a,n,dt,zeta,P,nP,SPA,SPI)
        call error(abs(SPA),SPAT,nP,aerror,merror)
        R = SPAT/abs(SPA)
        call decrlininterp(P,R,nP,Pf(IPf1:IPf2),Rf(IPf1:IPf2),NPf)
        ! write(unit=*, fmt="(A11,I4,A14,2F8.4)") "Error After",iter,"  Iterations: ",aerror,merror
        iter = iter + 1
        if (merror < minerr) then
            minerr = merror
            best = a*1.D0
        end if
    end do

    if (kpb>0) a = best*1.D0

    call fftw_destroy_plan(plan)
    call fftw_destroy_plan(iplan)

    deallocate(best)
    deallocate(SPA)
    deallocate(R)
    deallocate(Pf)
    deallocate(Rf)
    deallocate(SPI)
    deallocate(a0)
    deallocate(af)
    deallocate(afs)
    deallocate(f)

end subroutine fitspectra

subroutine adjustpeak(a,n,peak0) bind(c)
    integer, intent(in) :: n
    real(C_DOUBLE), intent(inout) :: a(n)
    real(C_DOUBLE), intent(in) :: peak0

    integer :: pk, i
    real(C_DOUBLE) :: peak
    
    pk = 0
    peak = 0.D0

    do i = 1, n, 1
        if ( abs(a(i))>peak ) then
            peak = a(i)
            pk = i
        end if
        if ( a(i)>peak0 ) then
            a(i) = peak0
        elseif ( a(i)<-peak0 ) then
            a(i) = -peak0
        end if
    end do

    if ( peak<peak0 ) then
        a(pk) = sign(peak0,peak)
    end if
    
end subroutine adjustpeak

subroutine adjustbaseline(a,n,dt) bind(c)
    integer, intent(in) :: n
    real(C_DOUBLE), intent(inout) :: a(n)
    real(C_DOUBLE), intent(in) :: dt

    real(C_DOUBLE), allocatable :: v(:), d(:)
    integer, allocatable :: tp(:)
    integer :: ntp = 4, pl = 2, ph, i

    ph = pl + ntp - 1

    allocate(v(n))
    allocate(d(n))
    allocate(tp(ntp))

    do i = 1, ntp, 1
        tp(i) = int(dble(i)/dble(ntp)*dble(n))
    end do

    call ratacc2vd(a,v,d,n,dt,0.d0,0.d0)
    call targetdc(a,d,n,tp,ntp,ph,pl,dt,0.d0,0.d0)

    deallocate(v)
    deallocate(d)
    deallocate(tp)

    return
end subroutine adjustbaseline

subroutine adjustspectra(acc,n,dt,zeta,P,nP,SPAT,a,tol,mit,kpb) bind(c)

    integer, intent(in) :: n, nP, mit, kpb
    real(C_DOUBLE), intent(in) :: dt, zeta, tol
    real(C_DOUBLE), intent(in) :: acc(n), P(nP), SPAT(nP)
    real(C_DOUBLE), intent(out) :: a(n)

    real(C_DOUBLE), allocatable :: SPA(:), dR(:), ra(:,:,:), best(:)
    real(C_DOUBLE), allocatable :: M(:,:), W(:,:)
    integer(C_INT), allocatable :: SPI(:), SPIp(:)

    real(C_DOUBLE) :: aerror, merror, peak0, minerr
    integer :: i,j,iter

    allocate(best(n))
    allocate(SPA(nP))
    allocate(dR(nP))
    allocate(SPI(nP))
    allocate(SPIp(nP))
    allocate(ra(n,nP,nP))
    allocate(M(nP,nP))
    allocate(W(n,nP))

    peak0 = maxval(abs(acc), dim=1)
    a = acc*1.D0

    iter = 1
    call spamixed(a,n,dt,zeta,P,nP,SPA,SPI)
    SPIp = SPI
    call errora(abs(SPA),SPAT,nP,aerror,merror)

    if (aerror <= tol .and. merror<=3.0d0*tol) return
    !write(unit=*, fmt="(A29,2F8.4)") "Initial Error: ",aerror,merror
    minerr = merror
    best = a*1.D0
    do while ( (aerror>tol .or. merror>3.0d0*tol) .and. iter<=mit )

        dR = SPA*(SPAT/abs(SPA)-1.d0)/SPAT

        !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i)
        do i = 1, nP, 1
            call wfunc( n,dt,SPI(i),P(i),zeta,W(:,i) )
        end do
        !$OMP END PARALLEL DO

        !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i, j)
        do i = 1, nP, 1
            do j = 1, nP, 1
                call ramixed(W(:,j),n,dt,zeta,P(i),ra(:,i,j))
                M(i,j) = ra(SPI(i),i,j)/SPAT(i)
                if ( i /= j ) then
                    M(i,j) = M(i,j)*0.618D0
                end if
            end do
        end do
        !$OMP END PARALLEL DO

        call leastsqs(M,dR,nP,nP)

        do i = 1, nP, 1
            a = a + dR(i)*W(:,i)
        end do
        
        call adjustpeak(a,n,peak0)
        ! call adjustbaseline(a,n,dt)
        ! call adjustpeak(a,n,peak0)

        SPIp = SPI
        call spamixed(a,n,dt,zeta,P,nP,SPA,SPI)
        call errora(abs(SPA),SPAT,nP,aerror,merror)
        if (merror<minerr) then
            minerr = merror
            best = a*1.D0
        end if
        !write(unit=*, fmt="(A11,I4,A14,2F8.4)") "Error After",iter,"  Iterations: ",aerror,merror
        iter = iter + 1

    end do
    if (kpb>0) a = best*1.D0

    deallocate(best)
    deallocate(SPA)
    deallocate(dR)
    deallocate(SPI)
    deallocate(SPIp)
    deallocate(ra)
    deallocate(M)
    deallocate(W)
end subroutine adjustspectra

subroutine flat_double(A,n,m,af)

    integer, intent(in) :: n, m
    real(8), intent(in) :: A(n,m)
    real(8), intent(out) :: af(n*m)

    integer :: i, j

    do i = 1, n, 1
        do j = 1, m, 1
            af(n*(j-1)+i) = A(i,j)
        end do
    end do

end subroutine flat_double

subroutine flat_int(A,n,m,af)

    integer, intent(in) :: n, m
    integer, intent(in) :: A(n,m)
    integer, intent(out) :: af(n*m)

    integer :: i, j

    do i = 1, n, 1
        do j = 1, m, 1
            af(n*(j-1)+i) = A(i,j)
        end do
    end do

end subroutine flat_int

subroutine repeat_double(a,n,m,Ar)
    integer, intent(in) :: n,m
    real(8), intent(in) :: a(n)
    real(8), intent(out) :: Ar(n*m)

    integer :: i

    do i = 1, m, 1
        Ar((i-1)*n+1:(i-1)*n+n) = a
    end do
    
end subroutine repeat_double

subroutine adjustspectra_endur(acc,n,dt,zeta,P,nP,DI,nD,SPAT0,a,tol,mit,kpb) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: adjustspectra_endur
    integer, intent(in) :: n, nP, mit, kpb, nD
    integer, intent(in) :: DI(nD)
    real(8), intent(in) :: dt, zeta, tol
    real(8), intent(in) :: acc(n), P(nP), SPAT0(nP)
    real(8), intent(out) :: a(n)

    real(8), allocatable :: SPA1(:), SPAT1(:), dR1(:), P1(:)
    real(8), allocatable :: SPA(:,:), dR(:,:), ra(:,:,:), best(:)
    real(8), allocatable :: M(:,:), W(:,:), SPAT(:,:)
    integer, allocatable :: SPI(:,:), SPI1(:)

    real(8) :: aerror, merror, peak0, minerr
    integer :: i,j,iter,nT

    nT = nP*nD

    allocate(best(n))
    allocate(SPA(nP,nD))
    allocate(SPAT(nP,nD))
    allocate(SPA1(nT))
    allocate(SPAT1(nT))
    allocate(dR(nP,nD))
    allocate(dR1(nT))
    allocate(SPI(nP,nD))
    allocate(SPI1(nT))
    allocate(P1(nT))
    allocate(ra(n,nT,nT))
    allocate(M(nT,nT))
    allocate(W(n,nT))

    call repeat_double(P,nP,nD,P1)

    peak0 = maxval(abs(acc), dim=1)
    a = acc

    do i = 1, nD, 1
        SPAT(:,i) = SPAT0*dble(DI(i))/dble(n)
        P1((i-1)*nP+1:(i-1)*nP+nP) = P
    end do
    
    call flat_double(SPAT,nP,nD,SPAT1)

    iter = 1

    call spamixed_endur(acc,n,dt,zeta,P,nP,DI,nD,SPA,SPI)
    call errora_endur(abs(SPA),SPAT,nP,nD,aerror,merror)

    if (aerror <= tol .and. merror<=3.0d0*tol) return
    !write(unit=*, fmt="(A29,2F8.4)") "Initial Error: ",aerror,merror
    minerr = merror
    best = a

    do while ( (aerror>tol .or. merror>3.0d0*tol) .and. iter<=mit )

        dR = SPA*(SPAT/abs(SPA)-1.d0)!/SPAT

        call flat_double(dR,nP,nD,dR1)
        call flat_double(SPA,nP,nD,SPA1)
        call flat_int(SPI,nP,nD,SPI1)

        !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i)
        do i = 1, nT, 1
            call wfunc( n,dt,SPI1(i),P1(i),zeta,W(:,i) )
        end do
        !$OMP END PARALLEL DO

        !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i, j)
        do i = 1, nT, 1
            do j = 1, nT, 1
                call ramixed(W(:,j),n,dt,zeta,P1(i),ra(:,i,j))
                M(i,j) = ra(SPI1(i),i,j)!/SPAT1(i)
                if ( i /= j ) then
                    M(i,j) = M(i,j)*0.7D0
                end if
            end do
        end do
        !$OMP END PARALLEL DO

        call leastsqs(M,dR1,nT,nT)

        do i = 1, nT, 1
            a = a + dR1(i)*W(:,i)
        end do

        ! call adjustpeak(a,n,peak0)
        ! call adjustbaseline(a,n,dt)
        ! call adjustpeak(a,n,peak0)

        call spamixed_endur(a,n,dt,zeta,P,nP,DI,nD,SPA,SPI)
        call errora_endur(abs(SPA),SPAT,nP,nD,aerror,merror)

        if (merror<minerr) then
            minerr = merror
            best = a
        end if
        !write(unit=*, fmt="(A11,I4,A14,2F8.4)") "Error After",iter,"  Iterations: ",aerror,merror
        iter = iter + 1

    end do
    if (kpb>0) a = best

    deallocate(best)
    deallocate(SPA)
    deallocate(SPAT)
    deallocate(SPA1)
    deallocate(dR)
    deallocate(dR1)
    deallocate(SPI)
    deallocate(ra)
    deallocate(M)
    deallocate(W)
    
end subroutine adjustspectra_endur

subroutine wfunc(n,dt,itm,P,zeta,wf) bind(c)

    integer, intent(in) :: n, itm
    real(8), intent(in) :: dt, P, zeta
    real(8), intent(out) :: wf(n)

    real(8) :: tm, w, f, tmp1, gamma, deltaT
    real(8), allocatable :: tmp2(:)
    integer :: i

    allocate(tmp2(n))

    tm = dble(itm-1)*dt
    w = TWO_PI/P
    f = 1.d0/P

    tmp1 = sqrt(1.d0-zeta**2)
    gamma = 1.178d0*(f*tmp1)**(-0.93d0)
    deltaT = atan(tmp1/zeta)/(w*tmp1)

    tmp2 = [ ( ( dble(i-1)*dt-tm+deltaT ),i=1,n ) ]

    wf = cos(w*tmp1*tmp2)*exp(-(tmp2/gamma)**2)

    deallocate(tmp2)

    return
end subroutine wfunc

function rsimple(pe0,p0,phi,t,zeta) bind(c)

    real(C_DOUBLE), intent(in) :: pe0,p0,phi,t,zeta
    real(C_DOUBLE) :: rsimple

    real(C_DOUBLE) :: we, w
    real(C_DOUBLE) :: we2, w2, zeta2
    real(C_DOUBLE) :: we3, w3, zeta3
    real(C_DOUBLE) :: we4, w4, zeta4
    real(C_DOUBLE) :: sinphi, cosphi
    real(C_DOUBLE) :: two = 2.d0, three = 3.d0, four = 4.d0
    real(C_DOUBLE) :: eight = 8.d0, ten = 1.d1, twelve = 12.d0

    we = TWO_PI/pe0
    w = TWO_PI/p0
    we2 = we*we
    w2 = w*w
    we3 = we2*we
    w3 = w2*w
    we4 = we3*we
    w4 = w3*w

    zeta2 = zeta*zeta
    zeta3 = zeta2*zeta
    zeta4 = zeta3*zeta

    sinphi = sin(phi)
    cosphi = cos(phi)

    rsimple = cos(we*t+phi) - exp(-w*zeta*t)*((four*w*we3*zeta &
        -four*w*we3*zeta3)*exp(w*zeta*t)*sin(we*t+phi)+((two*we4 &
        -two*w2*we2)*zeta2-two*we4+two*w2*we2)*exp(w*zeta*t)*cos(we*t &
        +phi)+sqrt(four*w2-four*w2*zeta2)*(four*cosphi*w*we2*zeta3 &
        +two*sinphi*we3*zeta2+(cosphi*w3-three*cosphi*w*we2) &
        *zeta-sinphi*we3+sinphi*w2*we)*sin(sqrt(four*w2-four*w2 &
        *zeta2)*t/two)+(eight*cosphi*w2*we2*zeta4+four*sinphi*w &
        *we3*zeta3+(two*cosphi*w4-ten*cosphi*w2*we2)*zeta2-four &
        *sinphi*w*we3*zeta+two*cosphi*w2*we2-two*cosphi*w4) &
        *cos(sqrt(four*w2-four*w2*zeta2)*t/two))/(eight*w2*we2*zeta4 &
        +(two*we4-twelve*w2*we2+two*w4)*zeta2-two*we4+four*w2*we2 &
        -two*w4)
end function rsimple

subroutine rnmknl(acc,n,dt0,zeta,P,ra,rv,rd,ku,model,tol,maxiter)  bind(c)
    !
    !  SUBROUTINE FOR COMPUTATION OF RESPONSE OF SDOF UNDER EARTHQUAKE
    !  acc -----  GROUND ACC
    !  N  -----  LENGTH OF acc
    !  DEL  -----  ORIGINAL TIME INTERVAL
    !  zeta  ----  DAMPING RATIO
    !  PARA  ----  (/K0,FY,K1,XP,FP,KP,FEFF,KEFF/)
    !     K0 -- initial stiffness
    !     FY -- yielding force
    !     K1 -- post-yielding stiffness
    !     XP -- commited displacement of last step
    !     FP -- commited force of last step
    !     KP -- commited stiffness of last step
    !   FEFF -- right-hand-side "force" of the equivalent static balance equation
    !   KEFF -- left-hand-side "stiffness" of the equivalent static balance equation

        integer, intent(in) :: n
        real(8), intent(in) :: dt0, zeta, tol
        real(8), intent(in) :: acc(n), P
        real(8), intent(out) :: ra(n),rv(n),rd(n),ku(n)
        integer, intent(in) :: maxiter, model

        real(8) :: rc(3),rl(3),ac,al,da,x
        real(8) :: k,c,w,feff,keff,beta,gamma,dt
        real(8) :: b1, b2, b3, b4, b5, b6, b7, b8
        integer :: i,j,r

        if ( dt0*MPR > P ) then
            r  = ceiling(MPR*dt0/P)
            dt = dt0/r
        else
            r = 1
            dt = dt0
        end if

        beta = 1.d0/4.d0
        gamma = 0.5d0

        b1 = 1.d0/(beta*dt*dt)
        b2 = 1.d0/(beta*dt)
        b3 = 1.d0/(beta*2.0d0)-1.d0
        b4 = gamma/(beta*dt)
        b5 = gamma/beta-1.d0
        b6 = 0.5d0*dt*(gamma/beta-2.0d0)
        b7 = dt*(1.0d0-gamma)
        b8 = dt*gamma

        w = 6.283185307179586d0/p
        k = w*w
        c = 2.0d0*zeta*w
        keff = b1+b4*c

        !  INITIATION
        para(8) = keff
        rl = [0.0d0,0.0d0,0.0d0]
        rc = [0.0d0,0.0d0,0.0d0]
        para(4:6) = [rl(1),0.0d0,k]
        al = 0.d0

        !  COMPUTATION OF RESPONSE
        do i = 1,n,1
            
            ac = acc(i)
            da = (ac - al)/dble(r)

            do j = 1, r, 1
                ac = al + da*j
        
                feff = ac+(b1*rl(1)+b2*rl(2)+b3*rl(3))&
                +c*(b4*rl(1)+b5*rl(2)+b6*rl(3))
                call solve_balance_equation(x,feff,model,tol,maxiter)

                rc(1) = x
                rc(2) = b4*(rc(1)-rl(1))-b5*rl(2)-b6*rl(3)
                rc(3) = ac-para(5)-c*rc(2)

                rl(1) = rc(1)
                rl(2) = rc(2)
                rl(3) = rc(3)
            end do

            ku(i) = para(13)
            ra(i) = -rc(3)+ac
            rv(i) = -rc(2)
            rd(i) = -rc(1)

            al = ac
        end do
        return
end subroutine rnmknl

subroutine newmark_nl(acc,n,dt0,zeta,P,SPA,SPI,SPV,SPD,SPE,model,tol,maxiter) bind(c)
    !  SUBROUTINE FOR COMPUTATION OF RESPONSE OF SDOF UNDER EARTHQUAKE
    !  acc -----  GROUND ACC
    !  N  -----  LENGTH OF acc
    !  DEL  -----  ORIGINAL TIME INTERVAL
    !  zeta  ----  DAMPING RATIO
    !  PARA  ----  (/K0,FY,K1,XP,FP,KP,FEFF,KEFF/)
    !     K0 -- initial stiffness
    !     FY -- yielding force
    !     K1 -- post-yielding stiffness
    !     XP -- commited displacement of last step
    !     FP -- commited force of last step
    !     KP -- commited stiffness of last step
    !   FEFF -- right-hand-side "force" of the equivalent static balance equation
    !   KEFF -- left-hand-side "stiffness" of the equivalent static balance equation

    integer, intent(in) :: n
    real(8), intent(in) :: dt0, zeta, tol
    real(8), intent(in) :: acc(n), P
    real(8), intent(out) :: SPA,SPV,SPD,SPE
    integer, intent(out) :: SPI
    integer, intent(in) :: maxiter, model

    real(8) :: rc(3),rl(3),ac,al,da,x,aa,rvl,ein
    real(8) :: k,c,w,feff,keff,beta,gamma,dt
    real(8) :: b1, b2, b3, b4, b5, b6, b7, b8
    integer :: i,j,r


    if ( dt0*MPR > P ) then
        r  = ceiling(MPR*dt0/P)
        dt = dt0/r
    else
        r = 1
        dt = dt0
    end if

    beta = 1.d0/4.d0
    gamma = 0.5d0

    b1 = 1.d0/(beta*dt*dt)
    b2 = 1.d0/(beta*dt)
    b3 = 1.d0/(beta*2.0d0)-1.d0
    b4 = gamma/(beta*dt)
    b5 = gamma/beta-1.d0
    b6 = 0.5d0*dt*(gamma/beta-2.0d0)
    b7 = dt*(1.0d0-gamma)
    b8 = dt*gamma

    w = 6.283185307179586d0/p
    k = w*w
    c = 2.0d0*zeta*w
    keff = b1+b4*c

    !  INITIATION
    para(8) = keff
    rl = [0.0d0,0.0d0,0.0d0]
    rc = [0.0d0,0.0d0,0.0d0]
    para(4:6) = [ rl(1),0.0d0,k ]
    al = 0.d0
    rvl = 0.d0

    SPA = 0.d0
    SPV = 0.d0
    SPD = 0.d0
    SPE = 0.d0
    SPI = 0

    !  COMPUTATION OF RESPONSE
    do i = 1,n,1
        
        ac = acc(i)
        da = (ac - al)/dble(r)

        do j = 1, r, 1
            ac = al + da*j
    
            feff = ac+(b1*rl(1)+b2*rl(2)+b3*rl(3))&
            +c*(b4*rl(1)+b5*rl(2)+b6*rl(3))
            call solve_balance_equation(x,feff,model,tol,maxiter)

            rc(1) = x
            rc(2) = b4*(rc(1)-rl(1))-b5*rl(2)-b6*rl(3)
            rc(3) = ac-para(5)-c*rc(2)

            rl(1) = rc(1)
            rl(2) = rc(2)
            rl(3) = rc(3)
        end do

        SPE = SPE+0.5d0*dt0*(rvl*al+rc(2)*ac)
        aa = -rc(3)+ac

        if (abs(aa) > abs(SPA)) then
            SPA=aa
            SPI=i
        endif

        if (abs(rc(2)) > abs(SPV)) then
            SPV=rc(2)
        endif

        if (abs(rc(1)) > abs(SPD)) then
            SPD=rc(1)
        endif

        al = ac
        rvl = rc(2)

    end do

    SPE = sqrt(2.d0*SPE)

    return
end subroutine newmark_nl

subroutine solve_balance_equation(x,feff,model,tol,maxiter) bind(c)

    integer, intent(in) :: maxiter, model
    real(8), intent(in) :: feff, tol
    real(8), intent(out) :: x

    real(8) :: x0, fval, kc

    !     估计变形值
    x0 = para(4)+(feff-para(7))/(para(6)+para(8))
    !     施加荷载
    para(7) = feff
    !     迭代求解变形
    call newton(x,x0,tol,maxiter,model)
    !     更新状态变量(para(4:6))
    call balance_func(fval,kc,x,model,.true.)

    return
end subroutine solve_balance_equation

subroutine newton(x,x0,tol,maxiter,model) bind(c)
    !     牛顿-拉普森法求解非线性方程
    !     输入参数：
    !           x0 : 解的初始估计值
    !           tol : 收敛容差(相对)
    !           maxiter : 最大迭代次数
    !           para : 一维数组, 包含非线性参数及历史状态变量, len(para)=np
    !     输出：
    !           x : 方程的解

    integer, intent(in) :: maxiter,model
    real(8), intent(out) :: x0,x
    real(8), intent(in) :: tol
    
    real(8) :: fval, kc
    integer :: iter
    integer :: status

    !     newton-rapheson method

    kc = 0.0d0
    fval = 0.0d0
    x = x0
    do iter = 1,maxiter,1

        call balance_func(fval,kc,x,model,.false.)

        if (kc .ne. 0.d0) then
            x = x0 - fval/kc
            if (abs(x - x0) < tol) then
                exit
            endif
            x0 = x
        else
            write(*,*) "Derivative was zero. Stop the Iteration."
            x = x0
            exit
        endif

    end do
end subroutine newton

subroutine balance_func(fval,kc,x,model,update) bind(c)
!     已知X, 计算方程F(X)=0中F(X)的值及F'(X)的值KC
!     输入参数：
!           X ---- 变形值
!           PARA,NP ---- 一维数组, 包含非线性参数及历史状态变量,
!                       LEN(PARA)=NP
!           MODEL ---- 非线性滞回模型
!                      0 ---- 双线性(DEFAULT)
!                      1 ---- 无残余变形双线性(卸载指向原点)
!                      2 ---- CLOUGH 退化双线性
!                      3 ---- 武田 退化三线性
!     输入结果：
!           FVAL ---- 平衡力残差
!           KC ---- 当前变形对应的刚度值(切线斜率)
!
    integer, intent(in) :: model
    real(8), intent(in) :: x
    logical, intent(in) :: update
    real(8), intent(out) :: fval, kc

    real(8) :: k0,fy,k1,xp,fp,kp,feff,dx,f,f_try,bup,bdn,uy
    real(8) :: alpha, dxp, xpeak, xmax, xmin, fmax, fmin, ku
    real(8) :: xu, k2, lpx, lpf, lux, luf, xl, xr, fl, fr
    integer :: status

    k0 = para(1)
    fy = para(2)
    k1 = para(3)

    xp = para(4)
    fp = para(5)
    kp = para(6)

    feff = para(7)

    alpha = para(9)
    dxp = para(12)
    ku = para(13)
    lpx = para(14)
    xmax = para(15)
    xmin = para(16)
    lpf = para(22)
    fmax = para(23)
    fmin = para(24)
    
    lux = para(29)
    luf = para(30)
    
    status = int(para(npara))
    
    uy = fy/k0
    dx = x-xp

    select case (model)
    case (1)
        if ( dx == 0.D0 ) then
            kc = kp
            f = fp
        elseif ( x == 0.D0 ) then
            f = 0.D0
            if ( dx > 0.D0 ) then
                ku = fmax/xmax
            else
                ku = fmin/xmin
            kc = ku
            end if
        else
            if ( x>0.D0 ) then
                if ( x <= xmax ) then
                    ku = fmax/xmax
                    f = ku*x
                    kc = ku
                else
                    kc = k1
                    xmax = x
                    f = fy + k1*(x-uy)
                    fmax = f
                    ku = fmax/xmax
                end if
            else
               if ( x >= xmin ) then
                    ku = fmin/xmin
                    f = ku*x
                    kc = ku
                else
                    kc = k1
                    xmin = x
                    f = -fy + k1*(x+uy)
                    fmin = f
                    ku = fmin/xmin
                end if
            end if
        end if
    case (2)
        if ( dx == 0.D0 ) then
            kc = kp
            f = fp
        else
            if ( status /= 0 .and. dx*dxp < 0.D0 ) then
                lux = xp
                luf = fp
            end if

            if ( x > xmax) then
                f = fy + ( x - uy ) * k1
                kc = k1
                status = 1
                xmax = x
                fmax = f
                ku = k0*(uy/xmax)**alpha
                lpx = x
                lpf = f
                lux = x
                luf = f
            elseif ( x < xmin) then
                f = -fy + ( x + uy ) * k1
                kc = k1
                status = 1
                xmin = x
                fmin = f
                ku = k0*(-uy/xmin)**alpha
                lpx = x
                lpf = f
                lux = x
                luf = f
            else if (status == 0) then
                ku = kp
                xu = lux - luf/ku

                if ( xu < lux ) then
                    xl = xu
                    xr = lux
                    fl = 0.D0
                    fr = luf
                else
                    xl = lux
                    xr = xu
                    fl = luf
                    fr = 0.D0
                end if

                if ( x >= xl .and. x<=xr ) then
                    kc = ku
                    f = kc*(x - xu)
                else if ( x > xr ) then
                    if (abs(xmax - xr) < 1.0D-12) then
                        k2 = kp
                    else
                        k2 = (fmax - fr)/(xmax - xr)
                    end if
                    kc = k2
                    f = fr + kc*(x - xr)
                    status = 2
                else
                    if (abs(xmin - xl) < 1.0D-12) then
                        k2 = kp
                    else
                        k2 = (fmin - fl)/(xmin - xl)
                    end if
                    kc = k2
                    f = fl + kc*(x - xl)
                    status = 2
                end if
            else if (status == 1) then
                if ( xp > 0.D0 ) then
                    ku = k0*(uy/xmax)**alpha
                    xu = xmax - fmax/ku

                    if ( x >= xu ) then
                        kc = ku
                        f = kc*(x - xu)
                        status = 0
                    else
                        k2 = fmin/(xmin - xu)
                        kc = k2
                        f = kc*(x - xu)
                        status = 2
                    end if

                    !lux = xp
                    !luf = fp
                else
                    ku = k0*(-uy/xmin)**alpha
                    xu = xmin - fmin/ku

                    if ( x <= xu ) then
                        kc = ku
                        f = kc*(x - xu)
                        status = 0
                    else
                        k2 = fmax/(xmax - xu)
                        kc = k2
                        f = kc*(x - xu)
                        status = 2
                    end if

                    !lux = xp
                    !luf = fp
                end if
            else if (status == 2) then
                if ( luf > 0.D0 ) then
                    ku = k0*(uy/xmax)**alpha
                else
                    ku = k0*(-uy/xmin)**alpha
                end if
                if ( dx*dxp > 0.D0 ) then
                    kc = kp
                    f = fp + kc*dx
                else
                    kc = ku
                    xu = lux - luf/ku
                    f = kc*(x - xu)
                    status = 0

                    !lux = xp
                    !luf = fp
                end if
            end if
        end if
    case default
        f_try = fp+dx*kp
        bup = fy+k1*(x-uy)
        bdn = -fy+k1*(x+uy)
        if ( dx == 0.D0 ) then
            kc = kp
            f = fp
        elseif (status == 0) then
            if (dx > 0.D0) then
                if ( f_try > bup ) then
                    f = bup
                    kc = k1
                    status = 1
                else
                    kc = k0
                    f = f_try
                end if
            else
                if ( f_try < bdn ) then
                    f = bdn
                    kc = k1
                    status = 1
                else
                    kc = k0
                    f = f_try
                end if
            end if
        elseif (status == 1) then
            if (dxp*dx>0) then
                f = f_try
                kc = k1
            else
                kc = k0
                f = fp + kc*dx
                status = 0
                if (f > bup) then
                    f = bup
                    kc = k1
                    status = 1
                elseif (f < bdn) then
                    f = bdn
                    kc = k1
                    status = 1
                end if
            end if
        end if
    end select
    
    if (update) then
        para(4) = x
        para(5) = f
        para(6) = kc
        para(12) = dx
        para(13) = ku
        para(14) = lpx
        para(15) = xmax
        para(16) = xmin
        para(22) = lpf
        para(23) = fmax
        para(24) = fmin
        para(29) = lux
        para(30) = luf
        para(npara) = dble(status)
    end if

    kc = kc+para(8)
    fval = f+para(8)*x-feff
    
    return
end subroutine balance_func

subroutine spmu(acc,n,dt,zeta,P,nP,SPA,SPI,SPV,SPD,SPE,mu,model,rtol,maxiter,uy,rk,alpha) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: spmu
    integer, intent(in) :: n, nP, model, maxiter
    real(8), intent(in) :: mu, dt, zeta, P(nP), rtol, rk, alpha
    real(8), intent(in) :: acc(n)
    real(8), intent(out) :: SPA(nP), SPV(nP), SPD(nP), SPE(nP)
    integer, intent(out) :: SPI(nP)
    real(8), intent(out) :: uy(nP)

    real(8) :: w, k0, k1, fy, mu1, pkd, tol
    integer :: i, j, pki

    do i = 1, nP, 1
        call newmark(acc,n,dt,zeta,P(i),SPA(i),SPI(i),SPV(i),SPD(i),SPE(i))
        uy(i) = abs(SPD(i))/mu
        tol = rtol*abs(SPD(i))

        w = TWO_PI/P(i)
        k0 = w*w
        k1 = k0*rk
        fy = k0*uy(i)
        para = 0.D0
        para(1:3) = [ k0,fy,k1 ]
        para(6) = k0
        para(9) = alpha
        para(13) = k0
        para(15:16) = [ uy,-uy ]
        para(23:24) = [ fy,-fy ]

        do j = 1, 100, 1
            call newmark_nl(acc,n,dt,zeta,P(i),SPA(i),SPI(i),SPV(i),SPD(i),SPE(i),model,tol,maxiter)
            mu1 = abs(SPD(i))/uy(i)

            if (abs(mu1-mu)<1.D-2) exit

            uy(i) = uy(i)*(mu1/mu)**0.5
            fy = k0*uy(i)
            para(2) = fy
            para(4:NPARA) = 0.D0
            para(6) = k0
            para(13) = k0
            para(15:16) = [ uy,-uy ]
            para(23:24) = [ fy,-fy ]
        end do

        SPA(i) = abs(SPA(i))
        SPV(i) = abs(SPV(i))
        SPD(i) = abs(SPD(i))
    end do

    return
end subroutine spmu

subroutine rnl(acc,n,dt,zeta,P,ra,rv,rd,ku,SM,cp) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: rnl
    integer, intent(in) :: n, SM
    real(8), intent(in) :: dt, zeta
    real(8), intent(inout) :: cp(8)
    real(8), intent(in) :: acc(n), P
    real(8), intent(out) :: ra(n),rv(n),rd(n),ku(n)

    integer :: i, model, maxiter
    real(8) :: w, k0, fy, uy, rk, k1, mu, muc, alpha
    real(8) :: pkd, pkf, tol

    mu = cp(1)
    rk = cp(2)
    alpha = cp(3)
    
    model = int(cp(8))

    maxiter = 1000
    tol = 1.D-9

    w = TWO_PI/P
    k0 = w*w
    k1 = k0*rk

    call r(acc,n,dt,zeta,P,ra,rv,rd,ku,2)
    pkd = abs(peak(rd,n))
    uy = pkd/mu
    pkf = pkd*k0
    para = 0.D0
    para(9) = alpha

    if ( SM == 0 ) then
        fy = k0*uy
        para(1:3) = [ k0,fy,k1 ]
        para(6) = k0
        para(13) = k0
        para(15:16) = [ uy,-uy ]
        para(23:24) = [ fy,-fy ]
        do i = 1, 1000, 1
            call rnmknl(acc,n,dt,zeta,P,ra,rv,rd,ku,model,tol,maxiter)
            pkd = abs(peak(rd,n))
            muc = pkd/uy
            if ( abs(muc - mu) < 1.D-2 ) exit
            uy = uy*(muc/mu)**0.5
            fy = k0*uy
            para(2) = fy
            para(4:NPARA) = 0.D0
            para(6) = k0
            para(13) = k0
            para(15:16) = [ uy,-uy ]
            para(23:24) = [ fy,-fy ]
        end do
    else
        fy = pkf/mu
        para(1:3) = [ k0,fy,k1 ]
        para(6) = k0
        para(13) = k0
        para(15:16) = [ uy,-uy ]
        para(23:24) = [ fy,-fy ]
        call rnmknl(acc,n,dt,zeta,P,ra,rv,rd,ku,model,tol,maxiter)
    end if
    
    cp(6) = fy
    cp(7) = fy/k0

    return
end subroutine rnl

subroutine r(acc,n,dt,zeta,P,ra,rv,rd,ku,SM) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: r
    integer, intent(in) :: n, SM
    real(8), intent(in) :: dt, zeta
    real(8), intent(in) :: acc(n), P
    real(8), intent(out) :: ra(n),rv(n),rd(n),ku(n)

    ku = TWO_PI*TWO_PI/(P*P)
    select case (SM)
        case (1)
            call rfreq(acc,n,dt,zeta,P,ra,rv,rd)
        case (2)
            call rnmk(acc,n,dt,zeta,P,ra,rv,rd)
        case (3)
            call rmixed(acc,n,dt,zeta,P,ra,rv,rd)
        case default
            call rmixed(acc,n,dt,zeta,P,ra,rv,rd)
    end select
end subroutine r

subroutine spectrum(acc,n,dt,zeta,P,np,SPA,SPI,SM) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: spectrum
    integer, intent(in) :: n, np, SM
    real(8), intent(in) :: dt,zeta
    real(8), intent(in) :: acc(n), P(np)
    real(8), intent(out) :: SPA(np)
    integer, intent(out) :: SPI(np)

    select case (SM)
        case (1)
            call spafreq(acc,n,dt,zeta,P,nP,SPA,SPI)
        case (2)
            call spanmk(acc,n,dt,zeta,P,nP,SPA,SPI)
        case (3)
            call spamixed(acc,n,dt,zeta,P,nP,SPA,SPI)
        case (4)
            call pspafreq(acc,n,dt,zeta,P,nP,SPA,SPI)
        case (5)
            call pspanmk(acc,n,dt,zeta,P,nP,SPA,SPI)
        case (6)
            call pspamixed(acc,n,dt,zeta,P,nP,SPA,SPI)
        case default
            call spamixed(acc,n,dt,zeta,P,nP,SPA,SPI)
    end select

    SPA = abs(SPA)

    return
end subroutine spectrum

subroutine spectrum_endur(acc,n,dt,zeta,P,np,DI,nD,SPA,SPI,SM) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: spectrum_endur
    integer, intent(in) :: n, np, SM, nD
    integer, intent(in) :: DI(nD)
    real(8), intent(in) :: dt,zeta
    real(8), intent(in) :: acc(n), P(np)
    real(8), intent(out) :: SPA(np*nD)
    integer, intent(out) :: SPI(np*nD)

    real(8), allocatable :: SPA1(:,:)
    integer, allocatable :: SPI1(:,:)

    allocate(SPA1(np,nD))
    allocate(SPI1(np,nD))

    select case (SM)
        case (1)
            call spafreq_endur(acc,n,dt,zeta,P,nP,DI,nD,SPA1,SPI1)
        case (2)
            call spanmk_endur(acc,n,dt,zeta,P,nP,DI,nD,SPA1,SPI1)
        case (3)
            call spamixed_endur(acc,n,dt,zeta,P,nP,DI,nD,SPA1,SPI1)
        case default
            call spamixed_endur(acc,n,dt,zeta,P,nP,DI,nD,SPA1,SPI1)
    end select

    SPA1 = abs(SPA1)

    call flat_double(SPA1,nP,nD,SPA)
    call flat_int(SPI1,nP,nD,SPI)

    deallocate(SPA1)
    deallocate(SPI1)

    return
end subroutine spectrum_endur

subroutine spectrumavd(acc,n,dt,zeta,P,np,SPA,SPI,SPV,SPD,SPE,SM) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: spectrumavd
    integer, intent(in) :: n, np, SM
    real(8), intent(in) :: dt,zeta
    real(8), intent(in) :: acc(n), P(np)
    real(8), intent(out) :: SPA(np), SPV(np), SPD(np), SPE(np)
    integer, intent(out) :: SPI(np)

    select case (SM)
        case (1)
            call spavdfreq(acc,n,dt,zeta,P,nP,SPA,SPI,SPV,SPD,SPE)
        case (2)
            call spavdnmk(acc,n,dt,zeta,P,nP,SPA,SPI,SPV,SPD,SPE)
        case (3)
            call spavdmixed(acc,n,dt,zeta,P,nP,SPA,SPI,SPV,SPD,SPE)
        case default
            call spavdmixed(acc,n,dt,zeta,P,nP,SPA,SPI,SPV,SPD,SPE)
    end select

    SPA = abs(SPA)
    SPV = abs(SPV)
    SPD = abs(SPD)
    
end subroutine spectrumavd

subroutine fitspectrum(acc,n,dt,zeta,P,nP,SPAT,a,tol,mit,fm,kpb) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: fitspectrum
    integer, intent(in) :: n, nP, mit, fm, kpb
    real(8), intent(in) :: dt, zeta, tol
    real(8), intent(in) :: acc(n), P(nP), SPAT(nP)
    real(8), intent(out) :: a(n)

    real(8), allocatable :: EP(:), ESPAT(:)

    if ( fm == 0 ) then
        allocate(EP(nP+2))
        allocate(ESPAT(nP+2))

        EP(2:nP+1) = P
        EP(1) = P(1)*0.5
        EP(nP+2) = P(nP)*1.5
        ESPAT(2:nP+1) = SPAT
        ESPAT(1) = SPAT(1)-(SPAT(2)-SPAT(1))/(P(2)-P(1))*P(1)*0.5d0
        ESPAT(nP+2) = SPAT(nP)+(SPAT(nP)-SPAT(nP-1))/(P(nP)-P(nP-1))*P(nP)*0.5d0

        call fitspectra(acc,n,dt,zeta,EP,nP+2,ESPAT,a,tol,mit,kpb)

        deallocate(EP)
        deallocate(ESPAT)
    else if ( fm == 1 ) then
        call adjustspectra(acc,n,dt,zeta,P,nP,SPAT,a,tol,mit,kpb)
    end if
end subroutine fitspectrum

subroutine initArtWave(a,n,dt,zeta,P,nP,SPT) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: initartwave
    integer, intent(in) :: n, nP
    real(C_DOUBLE), intent(in) :: dt, zeta, P(nP), SPT(nP)
    real(C_DOUBLE), intent(out) :: a(n)

    integer :: Nfft, i, j, k, NPf, IPf1, IPf2
    real(C_DOUBLE) :: phi, Ak, Saw
    real(C_DOUBLE), allocatable :: f(:), Pf(:), SPTf(:)
    complex(C_DOUBLE_COMPLEX), allocatable :: a0(:), af(:)
    
    type(C_PTR) :: plan

    Nfft = nextpow2(n)

    allocate(a0(Nfft))
    allocate(af(Nfft))
    allocate(f(Nfft))
    allocate(Pf(Nfft/2))
    allocate(SPTf(Nfft/2))

    af = (0.D0,0.D0)
    
    plan = fftw_plan_dft_1d(Nfft,af,a0,FFTW_BACKWARD,FFTW_ESTIMATE)

    call fftfreqs(Nfft,1.d0/dt,f)
    Pf(2:Nfft/2) = 1.d0/f(2:Nfft/2)
    Pf(1) = 100.d0*Pf(2)

    call decrfindfirst(Pf,Nfft/2,P(nP),IPf1)
    call decrfindlast(Pf,Nfft/2,P(1), IPf2)
    call decrlininterp(P,SPT,nP,Pf(IPf1:IPf2),SPTf(IPf1:IPf2),IPf2-IPf1+1)

    call random_seed()
    do k = IPf1, IPf2, 1
        call random_number(phi)
        phi = phi*TWO_PI
        Saw = 2.0D0*zeta*SPTf(k)**2.0D0/(pi*2.0D0*pi*f(k)*&
            (-2.0D0*log(-pi*log(0.85D0)/(2.0D0*pi*f(k)*(dt*n-dt)))))
        Ak = 2.0D0*sqrt(Saw*2.0D0*pi*(1.D0/dt)*Nfft/2)
        af(k) = Ak*(cmplx(cos(phi),sin(phi)))
        af(Nfft+2-k) = Ak*(cmplx(cos(phi),-sin(phi)))
    end do

    call fftw_execute_dft(plan,af,a0)
    a = real(a0(1:n))/Nfft
    call fftw_destroy_plan(plan)
    call fftw_cleanup()
    
    deallocate(a0)
    deallocate(af)
    deallocate(f)
    deallocate(Pf)
    deallocate(SPTf)

end subroutine initArtWave

end module eqs
