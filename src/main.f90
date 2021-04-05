!! Regularized Optimization for Hypers-spectral Analysis (ROHSA)
program ROHSA

  use mod_constants
  use mod_start
  use mod_inout
  use mod_rohsa
  use mod_functions
  use mod_array
  use mod_fits
  use mod_convert
  
  implicit none

  logical :: noise           !! if false --> STD map computed by ROHSA with lstd and ustd (if true given by the user)
  logical :: regul           !! if true --> activate regulation
  logical :: descent         !! if true --> activate hierarchical descent to initiate the optimization
  logical :: save_grid       !! save grid of fitted parameters at each step of the multiresolution process
  logical :: lym             !! if true --> activate 2-Gaussian decomposition for Lyman alpha nebula emission
  logical :: init_spec       !! if true --> use params mean spectrum with input
  logical :: init_grid       !! if true --> use fileinit to give the initialization of the last grid

  integer :: n_gauss         !! number of gaussian to fit
  integer :: m               !! number of corrections used in the limited memory matrix by LBFGS-B
  integer :: lstd            !! lower bound to compute the standard deviation map of the cube (if noise .eq. false)
  integer :: ustd            !! upper bound to compute the standrad deviation map of the cube (if noise .eq. false)
  integer :: iprint          !! print option 
  integer :: iprint_init     !! print option init
  integer :: maxiter         !! max iteration for L-BFGS-B alogorithm
  integer :: maxiter_init    !! max iteration for L-BFGS-B alogorithm (init mean spectrum)

  integer :: nl !nl = nline in cat SDC2
  integer :: ns, nf
  integer :: dxy, dv, cx, cy, cv
  integer, dimension(3) :: dim_array, dim_data

  real(xp) :: lambda_amp     !! lambda for amplitude parameter
  real(xp) :: lambda_mu      !! lamnda for mean position parameter
  real(xp) :: lambda_sig     !! lambda for dispersion parameter

  real(xp) :: lambda_var_amp !! lambda for variance amplitude parameter
  real(xp) :: lambda_var_mu  !! lambda for variance mean position parameter
  real(xp) :: lambda_var_sig !! lambda for variance dispersion parameter

  real(xp) :: lambda_lym_sig !! lambda for variance dispersion parameter

  real(xp) :: amp_fact_init  !! times max amplitude of additional Gaussian
  real(xp) :: norm_cube      !! normalize cube to have values around unity
  real(xp) :: sig_init       !! dispersion of additional Gaussian
  real(xp) :: lb_sig_init    !! lower bound sigma init
  real(xp) :: ub_sig_init    !! upper bound sigma init
  real(xp) :: lb_sig         !! lower bound sigma
  real(xp) :: ub_sig         !! upper bound sigma

  character(len=512) :: filename_parameters !! name of the parameters file (default parameters.txt)
  character(len=512) :: filename            !! name of the data file
  character(len=512) :: filename_cat        !! name of the data file
  character(len=512) :: fileout             !! name of the output result
  character(len=512) :: timeout             !! name of the output result
  character(len=512) :: filename_noise      !! name of the file with STD map (if noise .eq. true)
  character(len=512) :: filename_init_spec
  character(len=512) :: fileinit           !! name of the file with init last grid
  character(len=8)   :: init_option !!Init ROHSA with the mean or the std spectrum    

  character(len=512) :: fileout_cat, fileout_cat_rms !! name fitsfile output 

  real(xp) :: start, finish

  real(xp), dimension(:,:,:), allocatable :: data,data2  !! initial fits data 
  real(xp), dimension(:,:,:), allocatable :: data_init   !! initial fits data 
  real(xp), dimension(:,:), allocatable   :: std_map     !! standard deviation map fo the cube is given by the user 
  real(xp), dimension(:), allocatable     :: params_init !! standard deviation map fo the cube is given by the user 
  real(xp), dimension(:,:,:), allocatable :: grid_params !! parameters to optimize at final step (dim of initial cube)

  integer, dimension(:), allocatable  :: x, y, f !!ra, dec, freq
  real(xp), dimension(:), allocatable  :: conf !!confiance

  !Parameters to read fits file
  integer :: stat=0,uni=0,blocksize,naxes(3)
  character(len=80) comment
  character errtext*30
  logical :: undef
  
  real(kind=4), allocatable, dimension(:,:,:) :: array
  integer(kind=4), allocatable, dimension(:) :: fpix,lpix

  real(xp), dimension(:), allocatable :: mean_spect          !! mean spectrum of the observation
  real(xp) :: mean_std=0._xp

  integer :: i, k

  call cpu_time(start)

  !Print header and get filename in argument
  call get_command_argument(1, filename_parameters)
    
  !Default user parameters
  dxy = 8
  dv = 110
  ns = 1
  nf = 2
  n_gauss = 1
  lambda_amp = 1._xp
  lambda_mu = 1._xp
  lambda_sig = 1._xp
  lambda_var_amp = 0._xp
  lambda_var_mu = 0._xp
  lambda_var_sig = 0._xp
  lambda_lym_sig = 0._xp
  amp_fact_init = 2._xp/3._xp
  norm_cube = 1._xp
  sig_init = 5._xp
  lb_sig_init = 1._xp
  ub_sig_init = 100._xp
  lb_sig = 1._xp
  ub_sig = 100._xp
  maxiter_init = 15000
  maxiter = 800
  m = 10
  noise = .false.
  regul = .true.
  descent = .false.
  lstd = 1; ustd = 20
  init_option = "mean"
  iprint = -1
  iprint_init = -1
  save_grid = .true.
  lym = .false.
  init_grid = .false.
  init_spec = .false.
 
  !Read parameters
  call read_parameters(filename_parameters, filename, filename_cat, fileout, timeout, filename_noise, filename_init_spec, &
       n_gauss, lambda_amp, lambda_mu, lambda_sig, lambda_var_amp, lambda_var_mu, lambda_var_sig, lambda_lym_sig, &
       amp_fact_init, sig_init, lb_sig_init, ub_sig_init, lb_sig, ub_sig, init_option, maxiter_init, maxiter, &
       m, noise, regul, descent, lstd, ustd, iprint, iprint_init, save_grid, lym, init_grid, fileinit, init_spec, &
       dv, dxy, nl, ns, nf, norm_cube)
  
  !Call header
  call header()  

  !Load cat Benoit 
  print*, "Read catalogue Benoit SDC2"
  allocate(x(nl),y(nl),f(nl),conf(nl))
  call read_cat(filename_cat,x,y,f,conf,nl)
  
  print*, "filename = '",trim(filename),"'"

  !Print params
  print*, "fileout = '",trim(fileout),"'"
  print*, "timeout = '",trim(timeout),"'"
  
  print*, " "
  print*, "______Parameters_____"
  print*, "n_gauss = ", n_gauss
  
  print*, "lambda_amp = ", lambda_amp
  print*, "lambda_mu = ", lambda_mu
  print*, "lambda_sig = ", lambda_sig
  print*, "lambda_var_sig = ", lambda_var_sig
  print*, "amp_fact_init = ", amp_fact_init
  print*, "norm_cube = ", norm_cube
  print*, "sig_init = ", sig_init
  print*, "lb_sig_init = ", lb_sig_init
  print*, "ub_sig_init = ", ub_sig_init
  print*, "lb_sig = ", lb_sig
  print*, "ub_sig = ", ub_sig
  print*, "init_option = ", init_option
  print*, "maxiter_init = ", maxiter_init
  print*, "maxiter = ", maxiter
  print*, "lstd = ", lstd
  print*, "ustd = ", ustd
  print*, "noise = ", noise
  print*, "regul = ", regul
  print*, "descent = ", descent
  print*, "save_grid = ", save_grid
  print*, "init_spec = ", init_spec  
  print*, ""
  
  dim_array(1) = 2*dxy
  dim_array(2) = 2*dxy
  dim_array(3) = 2*dv

  !Check if nf < nl
  if (nf .ge. nl) stop ("nf must be lower than nl")

  !Read dim full cube SDC2
  call ftgiou(uni,stat)
  if (stat .gt. 0) then
     call ftgerr(stat,errtext)
     print *,'FITSIO Error Status =',stat,': ',errtext
     stop
  end if
  call ftdkopn(uni,filename,0,blocksize,stat)      
  call ftgkyj(uni,'NAXIS1',naxes(1),comment,stat)
  call ftgkyj(uni,'NAXIS2',naxes(2),comment,stat)
  call ftgkyj(uni,'NAXIS3',naxes(3),comment,stat)
  write(*,*) 'Data is',naxes(1),'by',naxes(2),'by',naxes(3)
  call ftclos(uni,stat)

  !Start loop
  do k=ns, nf
     !Allocate memory
     allocate(array(dim_array(1),dim_array(2),dim_array(3)))
     allocate(data(dim_array(3),dim_array(2),dim_array(1)))  
     allocate(mean_spect(dim_array(3)),params_init(3))
     allocate(grid_params(3*n_gauss,dim_array(2),dim_array(1))) !include noise map
     allocate(std_map(2*dxy,2*dxy))

     !Load data
     print*,"Read FITS file data cube"
     allocate(fpix(3), lpix(3))
     
     print*, "Position in cube = ", x(k),y(k),f(k)
     cx=x(k)+1; cy=y(k)+1; cv=f(k)+1
     
     !Check that cutout is in full cube
     if (cv < dv) goto 18
     if (cx < dxy) goto 18
     if (cy < dxy) goto 18
     if (cv > naxes(3)-dv) goto 18
     if (cx > naxes(1)-dxy) goto 18
     if (cy > naxes(2)-dxy) goto 18

     fpix(1) = cx-dxy
     lpix(1) = cx+dxy-1
     fpix(2) = cy-dxy
     lpix(2) = cy+dxy-1
     fpix(3) = cv-dv
     lpix(3) = cv+dv-1
          
     call readp_fits(filename,array,fpix,lpix)
     call roll_fits(array,data)
     dim_data = shape(data)

     !Compute std map noise
     call set_stdmap(std_map, data, lstd, ustd)
     ! mean_std = mean_2D(std_map,2*dxy,2*dxy)
     
     !Normalize cube
     data = data * norm_cube
     
     !Init Gaussian mean spectrum
     if (init_spec .eqv. .true.) then
        call mean_spectrum(data, mean_spect, dim_data(1), dim_data(2), dim_data(3))
        params_init(1) = amp_fact_init * maxval(mean_spect)
        params_init(2) = dv
        params_init(3) = sig_init
     end if
     
     !Call ROHSA subroutine
     call main_rohsa(grid_params, data, std_map, fileout, timeout, n_gauss, lambda_amp, lambda_mu, lambda_sig, &
          lambda_var_amp, lambda_var_mu, lambda_var_sig, lambda_lym_sig, amp_fact_init, sig_init, lb_sig_init, &
          ub_sig_init, lb_sig, ub_sig, maxiter_init, maxiter, m, noise, regul, descent, lstd, ustd, init_option, &
          iprint, iprint_init, save_grid, lym, init_grid, fileinit, data_init, params_init, init_spec)  
     
     !Unnomalize amplitude
     do i=1, n_gauss
        grid_params(1+3*(i-1),:,:) = grid_params(1+3*(i-1),:,:) / norm_cube
     end do
     
     !Write fits file
     fileout_cat = fileout(:len(trim(fileout))-5)//"_"//trim(str(k))//".fits"
     fileout_cat_rms = fileout(:len(trim(fileout))-5)//"_"//trim(str(k))//"_rms.fits"
     call writefits3D(fileout_cat,real(grid_params,kind=4),3*n_gauss,dim_data(2),dim_data(3))
     call writefits2D(fileout_cat_rms,real(std_map,kind=4),2*dxy,2*dxy)
     
18   deallocate(array,data,params_init,grid_params,std_map)
     deallocate(mean_spect)
     deallocate(fpix,lpix)

  end do

  !Call ender
  call ender()
  call cpu_time(finish)
  print '("Time = ",f6.3," seconds.")',finish-start
   
end program ROHSA
