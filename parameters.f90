MODULE parameters

   IMPLICIT NONE

    integer, parameter :: n1=512, n2=512, n3=512
    integer, parameter :: npe=128

    integer, parameter :: n1d=n1+2, n2d=n2, n3d=n3
    integer, parameter :: n3h0=n3/npe, n3h1=n3/npe+2, n3h2=n3/npe+4
    integer, parameter :: n3d1=n3d+2*npe                !for transposed f, all z-lev + 1-lev halos stacked

    integer, parameter :: izbot1=2,izbot2=3
    integer, parameter :: iztop1=n3h1-1,iztop2=n3h2-2
    

    integer, parameter :: ktx = n1/2,  kty =n2/2
    integer, parameter :: iktx= ktx+1, ikty=n2, iktyp=n2/npe

    double complex :: i = (0.,1.)
    double precision, parameter :: twopi=4.D0*asin(1.D0)

    double precision, parameter :: dom_x = 222000                            !Horizontal domain size (in m)
    double precision, parameter :: dom_z = 3000                              !Vertical   domain size (in m)
    double precision, parameter :: L1=twopi, L2=twopi, L3=twopi               !Domain size
    double precision, parameter :: dx=L1/n1,dy=L2/n2,dz=L3/n3                 !Cell dimensions  

    real, parameter :: ktrunc_x = twopi/L1 * float(n1)/3.           ! dimensional truncation wavenumber (x)
    real, parameter :: ktrunc_z = twopi/L3 * float(n3)/3.           ! dimensional truncation wavenumber (x)


    !Tags to specify run!
    !-------------------!

    !Gaussian wave initial condition
    double precision, parameter :: delta_a = 50.
    double precision, parameter :: xi_a = dom_z/(L3*delta_a)

    !C's and N's
    integer, parameter :: n_one = 1, n_two = 2
    double precision, parameter :: c_one = 1./( n_one*n_one - n_one*n_two*tanh(twopi*n_one)/tanh(twopi*n_two) )
    double precision, parameter :: c_two = 1./( n_two*n_two - n_one*n_two*tanh(twopi*n_two)/tanh(twopi*n_one) )

    integer, parameter :: barotropize = 0       !1: Waves only feel the effects of a barotropized flow;  0: waves and flow feel the same streamfunction (regular setup)
    integer, parameter :: bt_level = n3         !Level at which the barotropic streamfunction is defined (n3: top, 1: bottom, etc.)

    integer, parameter :: fixed_flow = 0        !1: Skip the psi-inversion steps
    integer, parameter :: passive_scalar = 0    !1: Set A and refraction to 0 and skip the LA -> A inversion. BR and BI become two (independent) passive scalars.
    integer, parameter :: ybj_plus = 1                  !1: B is L+A and A is recovered from B like psi is recovered from q (exception of the 1/4 factor). 0: Regular YBJ equation   
    
    integer, parameter :: no_waves = 0                  !1: Wave part ignored.
    integer, parameter :: no_feedback = 0               !1: Wave do not feedback on the flow; ): they do
    integer, parameter :: eady = 1                      !1: Eady version: add a bunch of terms
    integer, parameter :: eady_bnd = 1                  !1: Eady version: include the boundary terms (set NOT to zero only for testing purposes)

    integer, parameter :: no_dispersion=0
    integer, parameter :: no_refraction=0
    integer, parameter :: linear=0                      !1: set the nonlinear terms (advection) to 0. 
    integer, parameter :: init_wageo=0                  !1: Initialize wk with Ro*wak

    integer, parameter :: zero_aveB=1                   !1: Set B=LA vertical average to zero

    integer :: dealiasing=1         ! 1: Dealias, 0: don't. May not work though...

    !Should eventually plot both energies
!    integer, parameter :: plot_energy=1      !Use 1: energy_linear (equivalent to boussinesq including variable density, 2: energy_lipps)

    integer, parameter :: restoring_wind = 1              !1: Restore wind. 0: do not.
    double precision, parameter :: tau = 10*(3600*24)     !Dimensional wind-restoring time scale in seconds (days * 3600s/h 24h/d)

    !Initial structure!
    !-----------------!

    !Eady only
    integer, parameter :: ave_k=10              !Average wavenumber                                                                                          
    real, parameter ::    var_k=10.              !Variance of of the gaussian in wavenumbers                                                                                          
    double precision, parameter :: psi_0=0.1!1.     


    integer, parameter :: generic=1 
    integer, parameter :: init_vertical_structure=1
    integer, parameter :: linear_vert_structure=0          !1: LINEAR psi as in SB2013      2: kz=1 for all k's       OTHERWISE: kz consistent with QG.  
    integer, dimension(1) :: seed = 2     ! Seed for the random number generator                                                                                                    
    integer, parameter :: initial_k = 5                    !Peak of k-spec
    integer, parameter :: enveloppe = 0                    !1: Enveloppe allowing b=0 at the boundaries
    double precision, parameter :: z_env   = twopi/8!twopi/3.!twopi/8         !Center of the tanh enveloppe
    double precision, parameter :: sig_env = twopi/24!twopi/6.!twopi/24      !Width  of the tanh enveloppe
    double precision, parameter :: z0  = L3!L3/2           !Exponential N: ~exp[(z-z0)*N2_scale]        !Middle of the domain / Position of the tropopause (between 0 and L3)

    
    !Normalization at the tropopause!
    !-------------------------------!
    
    integer, parameter :: norm_trop = 1                !Normalize (1) (or not: 0) with the RMS value of U at the tropopause (overrides normalize=1, normalization from total energy)                                                        
    integer, parameter :: trop_height = n3/2         !Position (in stag) where we want U_RMS to be computed for normalization                                                      
    double precision, parameter :: URMS = 1.D0

    !Or normalization from total energy!
    !----------------------------------!

    integer, parameter :: norm_energy=1      !Use 1: energy_linear (equivalent to boussinesq including variable density, 2: energy_lipps to normalize fields.)
    integer, parameter :: normalize=0
    real, parameter :: e_init=0.175!0.007            !0.175 for  generic psi and 0.007 for TS-type of init give u'~1.
    real, parameter :: k_init=(2./3.)*e_init            
    real, parameter :: p_init=e_init-k_init    


    !Base-state!
    !----------!

    integer, parameter :: tropopause=1, exponential=2, constant_N=3, skewed_gaussian=4
    integer, parameter :: stratification = skewed_gaussian

    !Stratification = tropopause!
    integer, parameter :: fraction=128                   !If h#=150m, then fraction=133.333333~128
    double precision :: H_N = L3/fraction                          !Caracteristic length scale of N^2 for the TANH prof. (1/alpha...)
    double precision, parameter :: N_2_trop = 0.0001 !(2.*grav/10000.)*(1-t0_bot)/(1+t0_bot)             !BV frequency at the tropopause + in the tropsphere
    double precision, parameter :: N_2_stra = 0.0004                    !BV frequency in the stratosphere
    double precision, parameter :: gamma_N1=(sqrt(N_2_stra)-sqrt(N_2_trop))/(sqrt(N_2_stra)+sqrt(N_2_trop))       !This is alpha for N~1+alpha tanh(z/h)

    !Stratification = exponential!
    double precision, parameter :: N2_scale = dom_z/(twopi*450)   !N^2 ~ exp(N2_scale*(z-z0)), thus xi=H/h = 4000/(2pi*800) = 5/2pi 
!    double precision, parameter :: N0  =  sqrt(2e-5)             !Actual N is s^-1, not squared.  If ExpEady==1 ==> N0 = Nmax. 

    !Stratification = skewed gaussian!
    double precision, parameter :: N0 = 0.001550529072004        !Surface value of the fitted N2
    double precision, parameter :: N02_sg = 0.537713935783168
    double precision, parameter :: N12_sg = 2.684198470106461
    double precision, parameter :: sigma_sg = 0.648457170048730
    double precision, parameter :: z0_sg = 6.121537923499139
    double precision, parameter :: alpha_sg = -5.338431587899242
    double precision, parameter :: Xi = 0.155309488603754   !1./int_N2_nd        !Nondimensional parameter in front of the vQy term in Eady: H Ns2/ int(N^2) for the skewed gaussian

    !For the Eady case
    double precision :: U_mean(n3h0)                        !Eady: base-state velocity profile
    double precision :: Theta_y                             !Eady: base-state meriodional potential temperature gradient (constant, nondimensionalized)



   ! USEFUL INDEX !                                                                                                                          
   ! ------------ !                                                                                                                         

    integer :: ikx,iky,ikyp,izh0,izh1,izh2,izth   
    integer :: ix,iy,iz,izs
    integer :: kx,ky,kh2
    integer :: jj
    double precision :: x,y,z,zs

    ! USEFUL ARRAYS !                                                                                                                                                 
    ! ------------- !                                                                                                   
                                        
    integer, save :: kxa(iktx),kya(ikty)
    integer, save :: L(iktx,ikty)
    integer, save :: zath(n3)

    double precision, save :: xa(n1),ya(n2)
    double precision, save :: za(n3) ,zah2(n3h2) ,zah1(n3h1) ,zah0(n3h0)
    double precision, save :: zas(n3),zash2(n3h2),zash1(n3h1),zash0(n3h0)  !staggered version zasX(iz)=zaX(iz) + dz/2

    double precision, save :: r_1(n3h2),r_2(n3h2),r_3(n3h2)  !z-dependent r coefficients (r_1,2 unstag, while r_3 is stag)
    double precision, save :: r_1st(n3),r_2st(n3)            !Staggered and transposed versions of r_1 and r_2 for the elliptic equation.   (This could be a single value, not an array...) 
    double precision, save :: r_3u(n3h2),r_5u(n3h2)          !Unstaggered versions of r_3 and r_5 for the omega-equation verification.   
    double precision, save :: r_1ut(n3),r_2ut(n3)            !Unstaggered and transposed versions of r_1 and r_2 for the omega equation.  
    double precision, save :: r_3t(n3)                       !Transposed verion of r_3   (contains all n3 z-levels) for the pressure solver (still stag) 
    double precision, save :: r_3ut(n3)                      !Unstag version of r_3t
    double precision, save :: r_5ut(n3)                      
    double precision, save :: rho_st(n3)                     !Transposed verion of rho_s (contains all n3 z-levels) for diagnostics         (still stag) 
    double precision, save :: rho_ut(n3)                     !Transposed verion of rho_u (contains all n3 z-levels) for omega eqn (unstag)
    double precision, save :: a_ell_t(n3),b_ell_t(n3)        !coefficients of the elliptic equation for psi (LHQG)  --- transposed (for elliptic.f90)
    double precision, save :: a_ell_ut(n3)                   !coefficient of omega eqn --- transposed and UNstag - only computed for smooth_N for now
    double precision, save :: a_helm(n3),b_helm(n3)          !coefficients of the elliptic equation for phi ( QG )  --- transposed (for elliptic.f90)
    double precision, save :: a_ell(n3h2),b_ell(n3h2)        !coefficients of the elliptic equation for psi  --- not transposed (for setting the initial q in QG, and recover q from psi in LH)
    double precision, save :: a_ell_u(n3h2),b_ell_u(n3h2)    !coefficients of the elliptic equation for psi  --- unstaggered
    double precision, save :: rho_s(n3h2),rho_u(n3h2)        !Staggered and unstaggered versions of rho_0 (BS density)
    double precision, save :: r_1s(n3h2),r_2s(n3h2)          !Staggered version of r_1,2 (necessary to plot PV and initialize q) +  also conditions of integrability...
    double precision, save :: pi_0(n3h2)                     !For computing E_lh, Exner function's base-state (unstaggered) 
    double precision, save :: pi_0s(n3h2)                    !Staggered version of pi_0 (useful in two_exp BS
    double precision, save :: pi_0st(n3)                     !Staggered version of transposed pi_0 (useful in two_exp BS

    double precision, save :: eigen_vectors(n3,n3)               !Ouput = Eigenvectors(z^s_i,mode #)                                                                                     
    double precision, save :: eigen_values(n3)                   !Input = values on the diagonal. Output: eigenvalues in ascending order

    double precision, save :: czero(n3)                      !Temporary variable to estimate WPE creation: initial Az = C. Compute during the first call of A_solver_ybj_plus



    !I choose:
    !x N=0.03 s^-1 (to get N_troposphere ~ 0.01 right)
    !f=0.0001 s^-1 (real value of earth)
    !H=20 000m / 2pi
    !L=H/Ar 
    !U=1m/s 
    ! => Fr = pi/3
    ! => Ro = pi Ar
    !Keep in mind that: N/f = 1/Ar * (Ro/Fr).


    !Primary parameters!
    !------------------!

    double precision, parameter :: H_scale=dom_z/L3                  !Actual H in m ( z_real = H z' where z' in [0:L3]  is the nondim z.)
    double precision, parameter :: L_scale=dom_x/L1                  !Actual L in m ( x_real = L x' where x' in [0:2pi] is the nondim x.)
    double precision, parameter :: cor=1.2419D-04                    !Actual f = 0.0001 s^-1 (real value of planet Earth)
    double precision, parameter :: U_scale=0.01D0!0.025                       !Actual U in m/s (u_real = U u' where u' is the nondim velocity ur implemented in the code)
    double precision, parameter :: Uw_scale=0.1D0                       !Characteristic magnitude of wave velocity (wave counterpart to U_scale for flow)
    double precision, parameter :: Ar2 = (H_scale/L_scale)**2                                   !(1./64.)**2!(1./10.)**2 !0.01     !Aspect ratio squared = (H/L)^2     
    double precision, parameter :: Ro  = U_scale/(cor*L_scale)                                  !Rossby number  U/fL
    double precision, parameter :: Fr  = U_scale/(N0*H_scale)                                   !Froude number  U/N(z0)H
    double precision, parameter :: W2F = (Uw_scale/U_scale)**2                                  ! wave to flow velocity magnitude squared
    double precision, parameter :: Bu  = Fr*Fr/(Ro*Ro)                                          ! (Fr/Ro)^2 = Burger number 

    double precision, parameter :: delta_E = 6.D0                                                !Depth of the Ekman layer: 63 m
    double precision, parameter :: Ek  = delta_E/(Ro*H_scale)                                   !Ekman term = delta_E/(Ro H)

    double precision, parameter :: tau_wind = U_scale*tau/L_scale                          !Nondimensional wind-restoring timescale

    !Timestepping!
    !------------!

    real :: time=0.
    integer :: iter=0
    integer :: itermax=1000000000
    real :: maxtime=1000                      
    double precision, parameter :: delt= 0.002*dx!Ro/20.   !0.5*Bu*Ro/(2.*ktrunc_x*ktrunc_x)   !0.01*dx   !0.5*Bu*Ro/(2.*ktrunc_x*ktrunc_x) !0.25/ktrunc_x !0.5*Bu*Ro/(2.*ktrunc_x*ktrunc_x) 
    double precision, parameter :: gamma=1e-3                                  !Robert filter parameter


    !------------------------------!
    !--- Dissipation parameters ---!
    !------------------------------!

    !Assumes dissipation operator takes the form [ nuh1X*nabla^(2*ilap1X) + nuh2X*nabla^(2*ilap2X) ]. Suffix w is acting on waves.

    double precision, parameter :: coeff1  = 0.01
    double precision, parameter :: coeff2  = 10.
    double precision, parameter :: coeff1w = 0.
    double precision, parameter :: coeff2w = 10.

!    double precision, parameter :: coeff1  = 0.!0.01
!    double precision, parameter :: coeff2  = 0.!10.
!    double precision, parameter :: coeff1w = 0.!0.
!    double precision, parameter :: coeff2w = 0.!10.

    integer, parameter :: ilap1  = 2
    integer, parameter :: ilap2  = 6
    integer, parameter :: ilap1w = 2
    integer, parameter :: ilap2w = 6
    
    double precision, parameter :: nuh1   =  coeff1  * (64./(1.*n1)) **(4./3.) * (3./n1)**(2*(ilap1 -1))   !Dissipation operator 1, flow        
    double precision, parameter :: nuh2   =  coeff2  * (64./(1.*n1)) **(4./3.) * (3./n1)**(2*(ilap2 -1))   !Dissipation operator 2, flow             
    double precision, parameter :: nuh1w  =  coeff1w * (64./(1.*n1)) **(4./3.) * (3./n1)**(2*(ilap1w-1))   !Dissipation operator 1, wave             
    double precision, parameter :: nuh2w  =  coeff2w * (64./(1.*n1)) **(4./3.) * (3./n1)**(2*(ilap2w-1))   !Dissipation operator 2, wave             

    !Output!
    !------!

    integer, parameter :: out_etot   = 1, freq_etot   = INT(0.1/delt)!50!346!n3/64!n3!64!n3!50*n3/64      !Total energy                                                    
    integer, parameter :: out_we     = 1, freq_we     = INT(0.1/delt)!50!346!n3/64!n3!64!n3!50*n3/64      !Total energy                                                    
    integer, parameter :: out_conv   = 1, freq_conv   = freq_we      !Conversion terms in the potential energy equation.
    integer, parameter :: out_hspec  = 1, freq_hspec  = 1*freq_etot!n3/64!n3!freq_etot*10     !Horizontal energy spectrum at various heights 
    integer, parameter :: out_hspecw = 1, freq_hspecw = 1*freq_etot!n3/64!n3!freq_etot*10     !Horizontal energy spectrum at various heights 
    integer, parameter :: out_hg     = 0                 !Output geostrophic horizontal spectrum as well?
    integer, parameter :: out_vspec  = 0, freq_vspec =  freq_hspec
    integer, parameter :: out_vbuoy  = 0, freq_vbuoy =  freq_hspec
    integer, parameter :: out_vbuoyr = 0, freq_vbuoyr=  freq_etot
    integer, parameter :: out_ens    = 0, freq_ens   =  3*n3!freq_etot*10
    integer, parameter :: out_pv     = 0, freq_pv    =  3*n3!freq_etot*10

    integer, parameter :: out_ez     = 1, freq_ez    =  freq_etot        !E(z)  (freq has to be a multiple of that of etot) 
    integer, parameter :: out_wz     = 1, freq_wz    =  freq_we          !WE(z) (freq has to be a multiple of that of we)
    integer, parameter :: out_wshear = 1                                 !Calculate wave vertical shear
    integer, parameter :: out_rotz   = 0, freq_rotz  =  freq_etot 
    integer, parameter :: out_ensz   = 0, freq_ensz  =  3*n3!freq_ens
    integer, parameter :: out_pvz    = 0, freq_pvz   =  freq_pv
    integer, parameter :: out_cond   = 0, freq_cond  =  5*freq_etot!*10        !Plot the conditions of integrability of the balance equations.
    integer, parameter :: out_grow   = 0, freq_grow  =  5*freq_etot!*10        !Plot the conditions of integrability of the balance equations.
    integer, parameter :: out_omega  = 0, freq_omega =  5*freq_etot!*10        !Compute the QG ageotrophic vertical velocity wak and war
    integer, parameter :: out_condwz = 0, freq_condwz=  freq_omega!*10        !Plot the w_z condition (requires out_omega = 1)
    integer, parameter :: out_cont   = 0, freq_cont  =  freq_etot!*10        !Plot the anelastic divergence (should be 0 because of the proj method)



    

    !For conditions:
    double precision :: jump_region_width = 5.

    !For vspec
    integer, parameter :: num_couples=n1 !Should be well enough to have a good estimate of the vertical spectrum.
    integer, save :: x0(num_couples) 
    integer, save :: y0(num_couples) 
    integer, parameter :: variance_spectrum=1               !1: Plots variance spectrum (no rho(z) or other z-factors, just fields squared). 0: regular spectra

    integer, parameter :: parabolic_L_peak = 1              !1: Refines L_peak by using npt (odd number>=3) to fit a parabola.  0: Use the peak simply (causes discontinuous L_peak)
    integer, parameter :: npt = 5

    !For slices                                                                                                                     
    integer, parameter :: stag=1,unstag=2
    integer, parameter :: where_bz=unstag
    integer, parameter :: num_spec = 10

    integer, parameter :: height(num_spec)=[1, n3/8,  n3/4,  3*n3/8, n3/2, 5*n3/8, 3*n3/4, 7*n3/8,  n3-2 , n3]
    !                                       0   1      2       3      4      5       6        7      8      9    

    !Slices
    integer, parameter :: max_slices = 999     
    integer, parameter :: nfields  = 7         !Don't forget to change tag_slice_xz(nfields) accordingly in "mpi.f90"
    integer, parameter :: nfieldsw = 6         !Don't forget to change tag_slice_xz(nfields) accordingly in "mpi.f90"
    integer :: count_slice(nfields) = 0       !number of slices
    integer :: count_slicew(nfieldsw) = 0       !number of slices

    integer, parameter :: nvslices =3     
    integer :: yval(nvslices)=[n2,n2-n2/16,n2-n2/8]  !n2/2  

    integer :: hlvl(nfields)=[2,2,1,1,2,1,1]                                   
    integer :: hlvlw(nfieldsw)=[0,0,0,0,0,0]                                   

    integer, parameter :: bot_height = n3-34!1
    integer, parameter :: mid_height = n3-17!n3/2
    integer, parameter :: top_height = n3!-9 !n3-1

    integer, parameter :: out_slab = 0, freq_slab = 1
    integer, parameter :: slab_mype   = npe/2-1 
    integer :: count_eta = 0
    
                                              !halo levels (u=2,zz=1...)                                                                                                                                                     
    integer :: id_field                       !dummy index to differenciate fields plotted  

    integer, parameter :: out_slice   = 1, freq_slice =  10*freq_etot
    integer, parameter :: out_slicew  = 1, freq_slicew=  10*freq_etot
    integer, parameter :: out_eta     = 0, freq_eta   =  freq_hspec
    integer, parameter :: out_tspec   = 0

    !Restart
    integer :: count_restart = 0                                 !when dumping: restart file number 
    integer, parameter :: dump = 0, freq_dump = freq_slice*10    !dump = 1 means you dump, every "freq_dump" timestep
    integer, parameter :: restart = 1                            !restart = 1 start from file
    integer, parameter :: restart_no = 15                         !Restart file number (from 0 to 99)
    character(len = 64), parameter :: floc='../../dE60_dt0.01_512_7/output/'   !Location of the restart file (when restarting only: dumping in local output/ folder)


    !Filtering of A modes
    integer, parameter :: filter_A=1, freq_filter_A=1!*freq_etot
    integer, parameter :: print_A=1, freq_print_A=1*freq_we
    integer :: count_A=0
    double precision, parameter :: YBJ_criterion =3! 100000.           !Tolerate modes with (Nkh/fkz)^2 < YBJ_criterion.



END MODULE parameters
