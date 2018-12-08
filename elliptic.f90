MODULE elliptic

USE parameters
USE fft
USE mpi

IMPLICIT NONE

CONTAINS
  

    SUBROUTINE helmholtzdouble(p,ft,b_bot,b_top)  

      double complex, dimension(iktx,ikty,n3h1),  intent(out) :: p    !Pressure in usual z-parallelization                                                                       
      double complex, dimension(iktx,n3d1, iktyp)             :: pt   !Transposed (ky-parallelization) pressure                                                   
      double complex, dimension(iktx,n3, iktyp), intent(in) :: ft   !Transposed (ky-parallelization) right-hand side               
      double complex, dimension(iktx,ikty), intent(in) :: b_bot,b_top   !Top and bottom buoyancy terms             

      integer :: step1,step2                               !Intermediate steps in calculation of indices
      integer :: info                                      ! Returns error for sgtsv                                                                                           

      double precision :: d(n3),dl(n3-1),du(n3-1)          !diagonal value   
      double precision :: br(n3), bi(n3)                 !real and imaginary parts of rhs                                                                      

                                                                                                      
      DO ikx=1,iktx
         kx=kxa(ikx)
         DO ikyp=1,iktyp
            iky=ikyp+iktyp*mype
            ky=kya(iky)

            kh2=kx*kx + ky*ky

            if(kh2/=0 .and. L(ikx,iky)==1 ) then
               
               do iz=1,n3
                  
                  br(iz)=dz*dz*DBLE(  ft(ikx,iz,ikyp) )  
                  bi(iz)=dz*dz*DIMAG( ft(ikx,iz,ikyp) )  

               end do

               !--- For nonzero w* at the top/bot, we must add an extra to the RHS ---!
               !----------------------------------------------------------------------!

               iz=1
               
               br(iz) = br(iz) + (a_helm(iz) - 0.5*b_helm(iz)*dz)* DBLE(b_bot(ikx,iky))*dz
               bi(iz) = bi(iz) + (a_helm(iz) - 0.5*b_helm(iz)*dz)*DIMAG(b_bot(ikx,iky))*dz

               iz=n3
               
               br(iz) = br(iz) - (a_helm(iz) + 0.5*b_helm(iz)*dz)* DBLE(b_top(ikx,iky))*dz
               bi(iz) = bi(iz) - (a_helm(iz) + 0.5*b_helm(iz)*dz)*DIMAG(b_top(ikx,iky))*dz

               !----------------------------------------------------------------------!

 

               do iz=1,n3
                  if(iz==1) then
                     d(iz) = -a_helm(iz)-0.5*b_helm(iz)*dz-kh2*dz*dz
                    du(iz) =  a_helm(iz)+0.5*b_helm(iz)*dz
                  elseif(iz==n3) then
                     d(iz) = -a_helm(iz)+0.5*b_helm(iz)*dz-kh2*dz*dz
                  dl(iz-1) =  a_helm(iz)-0.5*b_helm(iz)*dz
                  else !1<iz<n3                                                                                                                                                                                                             
                     d(iz) = -2*a_helm(iz)-kh2*dz*dz
                    du(iz) =  a_helm(iz)+0.5*b_helm(iz)*dz
                  dl(iz-1) =  a_helm(iz)-0.5*b_helm(iz)*dz
                  end if
               end do

               
               call DGTSV( n3, 1, dl, d, du, br, n3, info )
               if(info/=0) write(*,*) "problem in helmdouble", info
               
               do iz=1,n3
                  if(iz==1) then
                     d(iz) = -a_helm(iz)-0.5*b_helm(iz)*dz-kh2*dz*dz
                    du(iz) =  a_helm(iz)+0.5*b_helm(iz)*dz
                  elseif(iz==n3) then
                     d(iz) = -a_helm(iz)+0.5*b_helm(iz)*dz-kh2*dz*dz
                  dl(iz-1) =  a_helm(iz)-0.5*b_helm(iz)*dz
                  else !1<iz<n3                                                                                                                                                                                                             
                     d(iz) = -2*a_helm(iz)-kh2*dz*dz
                    du(iz) =  a_helm(iz)+0.5*b_helm(iz)*dz
                  dl(iz-1) =  a_helm(iz)-0.5*b_helm(iz)*dz
                  end if
               end do
               
               call DGTSV( n3, 1, dl, d, du, bi, n3, info )
               if(info/=0) write(*,*) "problem in helmdouble (imag)", info


               !*** Put the solution in pt ***!    
               
               do izth=1,n3d1
                  
                  step1=(izth-1)/n3h1   !mype                                                       
                  step2=izth-n3h1*step1 !position in mype                                           
                  iz=step1*n3h0+step2-1 !-1 factor since h1 fields                                  
                  
                  if(izth/=1 .and. izth/=n3d1) then
                     pt(ikx,izth,ikyp)=br(iz)+i*bi(iz)
                  else
                     pt(ikx,izth,ikyp)=(0.,0.)
                  end if
                  
               end do
      
            else
            
               !set kh=0 modes to zero...                                                                                                                                    

                DO izth=1,n3d1
                   pt(ikx,izth,ikyp)=(0.,0.)
                END DO

             
             end if

         END DO
      END DO


      

      



      !*********** Transposition to z-parallelized p ***************!                                                                                                               

      call mpitranspose(pt,iktx,n3d1,iktyp,p,ikty,n3h1)


    END SUBROUTINE HELMHOLTZDOUBLE






    SUBROUTINE psi_solver(psik,qt)

      ! This subroutines solves the elliptic equation a_ell(z) d^2 psi/dz^2 + b_ell(z) d psi/dz - kh2 psi = q for the streamfunction psi. 

      double complex, dimension(iktx,ikty,n3h1),  intent(out) :: psik   !Pressure in usual z-parallelization                                                                       
      double complex, dimension(iktx,n3d1, iktyp)             :: psit   !Transposed (ky-parallelization) pressure                                                   
      double complex, dimension(iktx,n3, iktyp), intent(in) :: qt       !Transposed (ky-parallelization) right-hand side               

      integer :: step1,step2                               !Intermediate steps in calculation of indices
      integer :: info                                      ! Returns error for sgtsv                                                                                           

      double precision :: d(n3),dl(n3-1),du(n3-1)          !diagonal value   
      double precision :: br(n3), bi(n3)                 !real and imaginary parts of rhs                                                                      

                                                                                                                                                                          
      DO ikx=1,iktx
         kx=kxa(ikx)
         DO ikyp=1,iktyp
            iky=ikyp+iktyp*mype
            ky=kya(iky)

            kh2=kx*kx + ky*ky

            if(kh2/=0 .and. L(ikx,iky)==1 ) then
               
               do iz=1,n3
                  
                  br(iz)=dz*dz*DBLE(  qt(ikx,iz,ikyp) )  
                  bi(iz)=dz*dz*DIMAG( qt(ikx,iz,ikyp) )  

               end do

               do iz=1,n3
                  if(iz==1) then                     
                     d(iz) = -( rho_ut(iz)*a_ell_ut(iz)/rho_st(iz) + kh2*dz*dz )
                    du(iz) = rho_ut(iz)*a_ell_ut(iz)/rho_st(iz)
                  elseif(iz==n3) then
                     d(iz) = -( rho_ut(iz-1)*a_ell_ut(iz-1)/rho_st(iz) + kh2*dz*dz )
                  dl(iz-1) = rho_ut(iz-1)*a_ell_ut(iz-1)/rho_st(iz)
                  else !1<iz<n3
                     d(iz) = -( (rho_ut(iz)*a_ell_ut(iz) + rho_ut(iz-1)*a_ell_ut(iz-1))/rho_st(iz) + kh2*dz*dz )
                    du(iz) = rho_ut(iz)*a_ell_ut(iz)/rho_st(iz)
                  dl(iz-1) = rho_ut(iz-1)*a_ell_ut(iz-1)/rho_st(iz)
                  end if
               end do
               
               call DGTSV( n3, 1, dl, d, du, br, n3, info )
               if(info/=0) write(*,*) "problem in helmdouble", info
               
               do iz=1,n3
                  if(iz==1) then                     
                     d(iz) = -( rho_ut(iz)*a_ell_ut(iz)/rho_st(iz) + kh2*dz*dz )
                    du(iz) = rho_ut(iz)*a_ell_ut(iz)/rho_st(iz)
                  elseif(iz==n3) then
                     d(iz) = -( rho_ut(iz-1)*a_ell_ut(iz-1)/rho_st(iz) + kh2*dz*dz )
                  dl(iz-1) = rho_ut(iz-1)*a_ell_ut(iz-1)/rho_st(iz)
                  else !1<iz<n3
                     d(iz) = -( (rho_ut(iz)*a_ell_ut(iz) + rho_ut(iz-1)*a_ell_ut(iz-1))/rho_st(iz) + kh2*dz*dz )
                    du(iz) = rho_ut(iz)*a_ell_ut(iz)/rho_st(iz)
                  dl(iz-1) = rho_ut(iz-1)*a_ell_ut(iz-1)/rho_st(iz)
                  end if
               end do

               call DGTSV( n3, 1, dl, d, du, bi, n3, info )
               if(info/=0) write(*,*) "problem in helmdouble (imag)", info
               
               !*** Put the solution in pt ***!    
               
               do izth=1,n3d1
                  
                  step1=(izth-1)/n3h1   !mype                                                       
                  step2=izth-n3h1*step1 !position in mype                                           
                  iz=step1*n3h0+step2-1 !-1 factor since h1 fields                                  
                  
                  if(izth/=1 .and. izth/=n3d1) then
                     psit(ikx,izth,ikyp)=br(iz)+i*bi(iz)
                  else
                     psit(ikx,izth,ikyp)=(0.,0.)
                  end if
                  
               end do
               

            else
            
               !set kh=0 modes to zero...                                                                                                                                    

                DO izth=1,n3d1
                   psit(ikx,izth,ikyp)=(0.,0.)
                END DO

             
             end if

         END DO
      END DO

      !*********** Transposition to z-parallelized p ***************!                                                                                                               

      call mpitranspose(psit,iktx,n3d1,iktyp,psik,ikty,n3h1)


    END SUBROUTINE PSI_SOLVER
    






    SUBROUTINE omega_equation(wak,qt)

      ! This subroutines solves the elliptic equation a_ell(z) d^2 psi/dz^2 + b_ell(z) d psi/dz - kh2 psi = q for the streamfunction psi. 

      double complex, dimension(iktx,ikty,n3h1),  intent(out) :: wak   !vertical velocity in usual z-parallelization                                   
      double complex, dimension(iktx,n3d1, iktyp)             :: wat   !Transposed (ky-parallelization) wak                                                   
      double complex, dimension(iktx,n3, iktyp), intent(in) :: qt      !Transposed (ky-parallelization) right-hand side               

      integer :: step1,step2                               !Intermediate steps in calculation of indices
      integer :: info                                      ! Returns error for sgtsv                                                                                           

      double precision :: d(n3-1),dl(n3-2),du(n3-2)          !diagonal value   
      double precision :: br(n3-1), bi(n3-1)                 !real and imaginary parts of rhs                                                                      

                                                                                                                                                                          
      DO ikx=1,iktx
         kx=kxa(ikx)
         DO ikyp=1,iktyp
            iky=ikyp+iktyp*mype
            ky=kya(iky)

            kh2=kx*kx + ky*ky

            if(kh2/=0 .and. L(ikx,iky)==1 ) then
               
               do iz=1,n3-1
                  
                  br(iz)=dz*dz*DBLE(  qt(ikx,iz,ikyp) )  
                  bi(iz)=dz*dz*DIMAG( qt(ikx,iz,ikyp) )  

               end do

               do iz=1,n3-1
                  if(iz==1) then                     
                     d(iz) = -(rho_ut(iz)*(rho_st(iz+1)+rho_st(iz))/(rho_st(iz+1)*rho_st(iz))  + dz*dz*kh2/a_ell_ut(iz))
!-(rho_ut(iz)/(rho_st(iz+1)+rho_st(iz)) + dz*dz*kh2/a_ell_ut(iz))      WRONG
                    du(iz) = rho_ut(iz+1)/rho_st(iz+1)
                  elseif(iz==(n3-1)) then
                     d(iz) = -(rho_ut(iz)*(rho_st(iz+1)+rho_st(iz))/(rho_st(iz+1)*rho_st(iz))  + dz*dz*kh2/a_ell_ut(iz))
                  dl(iz-1) = rho_ut(iz-1)/rho_st(iz)
                  else !1<iz<n3-1
                     d(iz) = -(rho_ut(iz)*(rho_st(iz+1)+rho_st(iz))/(rho_st(iz+1)*rho_st(iz))  + dz*dz*kh2/a_ell_ut(iz))
                    du(iz) = rho_ut(iz+1)/rho_st(iz+1)
                  dl(iz-1) = rho_ut(iz-1)/rho_st(iz)
                  end if
               end do
               
               call DGTSV( n3-1, 1, dl, d, du, br, n3-1, info )
               if(info/=0) write(*,*) "problem in helmdouble", info

               do iz=1,n3-1
                  if(iz==1) then                     
                     d(iz) = -(rho_ut(iz)*(rho_st(iz+1)+rho_st(iz))/(rho_st(iz+1)*rho_st(iz))  + dz*dz*kh2/a_ell_ut(iz))
                    du(iz) = rho_ut(iz+1)/rho_st(iz+1)
                  elseif(iz==(n3-1)) then
                     d(iz) = -(rho_ut(iz)*(rho_st(iz+1)+rho_st(iz))/(rho_st(iz+1)*rho_st(iz))  + dz*dz*kh2/a_ell_ut(iz))
                  dl(iz-1) = rho_ut(iz-1)/rho_st(iz)
                  else !1<iz<n3-1
                     d(iz) = -(rho_ut(iz)*(rho_st(iz+1)+rho_st(iz))/(rho_st(iz+1)*rho_st(iz))  + dz*dz*kh2/a_ell_ut(iz))
                    du(iz) = rho_ut(iz+1)/rho_st(iz+1)
                  dl(iz-1) = rho_ut(iz-1)/rho_st(iz)
                  end if
               end do
               
               call DGTSV( n3-1, 1, dl, d, du, bi, n3-1, info )
               if(info/=0) write(*,*) "problem in helmdouble (imag)", info
               
               
               !*** Put the solution in pt ***!    
               
               do izth=1,n3d1
                  
                  step1=(izth-1)/n3h1   !mype                                                       
                  step2=izth-n3h1*step1 !position in mype                                           
                  iz=step1*n3h0+step2-1 !-1 factor since h1 fields                                  
                  
                  if(izth/=1 .and. izth/=n3d1 .and. iz/=n3) then    !Difference with psi_solver, if iz=n3, w=0. I didn't even solve the equation for that point...
                     wat(ikx,izth,ikyp)=br(iz)+i*bi(iz)
                  else
                     wat(ikx,izth,ikyp)=(0.,0.)
                  end if
                  
               end do
               

            else
            
               !set kh=0 modes to zero...                                                                                                                                    

                DO izth=1,n3d1
                   wat(ikx,izth,ikyp)=(0.,0.)
                END DO

             
             end if

         END DO
      END DO

      !*********** Transposition to z-parallelized p ***************!                                                                                                               

      call mpitranspose(wat,iktx,n3d1,iktyp,wak,ikty,n3h1)


    END SUBROUTINE OMEGA_EQUATION




    SUBROUTINE compute_eigen

      !This subroutine computes the eigenvalues and vectors for A and LA. LAPACK actually solves for L' = - (N0 H_scale / cor)^2 L,                                
      !Where the minus sign implies that the eigenvalues are positive and thus the eigenvectors are sorted from the gravest to highest mode.          


      !Eigenvalues output of LAPACK == Em, then                                                                                                   
      !Dimensional RDR = N0 H_scale / sqrt(Em) f    (Rossby deformation radius: LA == - (1/RDR)^2 A)                                     
      !Dimensional kappa = 1/RDR                    (Where LA == - kappa^2 A)                                                                

      double precision :: ds(n3-1)          !diagonal and sub/super diagonal values                                                      

      double precision :: WORK(2*n3-2)          !diagonal and sub/super diagonal values for DSTEV                                                   
      integer :: info                                      ! Returns error for LAPACK                                                   
            
      !Validate orthonormality of the modes                                                                                             
      integer :: iz1,iz2
      double precision :: leftover

      integer :: m


      !Initial condition as a function of z and m
      double precision :: CRz(n3),BRz(n3),ARz(n3)
      double precision :: CRm(n3),BRm(n3),ARm(n3)

      double precision :: z

      double precision :: ave_A,ave_B,ave_C
      double precision :: xi_a,xi_b,xi_c
      double precision :: delta_a,delta_b,delta_c

      double precision :: cons

      !Define matrix L'                                                                                                                            
      !Compute center diagonal:
      do iz=1,n3
         if(iz==1) then
            eigen_values(iz) = -(1./r_2ut(iz))
         elseif(iz==n3) then
            eigen_values(iz) = -(1./r_2ut(iz-1))
         else !1<iz<n3
            eigen_values(iz) = -(1./r_2ut(iz) + 1./r_2ut(iz-1))
         end if
      end do

      !Compute lower and upper diagonals: ds_i = 1/N^2(z^u_i) from i = 1 to N-1
      do iz=1,n3-1
         ds(iz) = 1./r_2ut(iz)
      end do

      !Give LAPACK -L = - d/dz'( (1/N')^2 d/dz') so that eigenvalues are positive and thus eigenvectors are sorted from gravest to highest                         
      eigen_values=-eigen_values/(dz*dz)
      ds=-ds/(dz*dz)

      !Calculate the eigenvectors and eigenvalues using LAPACK                                                                    
      call DSTEV( 'V', n3, eigen_values, ds, eigen_vectors, n3, WORK, INFO )       !DSTEV( JOBZ, N, D, E, Z, LDZ, WORK, INFO )                                                
      if(info/=0) write(*,*) "problem in compute_eigen", info

      !LAPACK solved the problem -LA = E A, our actual eigen values are kappa = sqrt(E) ---- The original problem is LA = - kappa^2 A                                         
      eigen_values = sqrt(abs(eigen_values))

      !Reinforce the first (constant) mode:                                                                                                                                   
!      eigen_values(1) = 0.
!      do iz=1,n3-1
!         eigen_vectors(iz,1) = 1/sqrt(1.*n3)
!      end do





      !********************************************************!                                                                                                            
      !* Print the first baroclinic Rossby deformation radius *!                                                                                                        
      !* the eigen values and eigen vectors and validate them *!                                                                                                           
      !********************************************************!                                                                                                                  


      if(mype==0) then

         write(*,*) "Rossby deformation radius (km) = ",N0*H_scale/(cor*eigen_values(2)*1000)

         !Print eigenvalues and eigenvectors of first 3 modes
         open (unit = 154679, file = "eigen.dat")
         do iz=1,n3
            write(154679,"(1x,E12.5,1x,E12.5,1x,E12.5,1x,E12.5,1x,E12.5)") zas(iz),eigen_values(iz),eigen_vectors(iz,1),eigen_vectors(iz,2),eigen_vectors(iz,3)
         end do

         !Test: verify that the modes are orthonormal
         do iz=1,n3
            do iz1=1,n3

               leftover=0.
               do iz2=1,n3
                  leftover = leftover + eigen_vectors(iz2,iz1)*eigen_vectors(iz2,iz)
               end do

               if(abs(leftover)>1e-14 .and. iz1/=iz) write(*,*) "Problem with orthonormality: modes (iz1,iz) have leftover=",leftover,iz1,iz

            end do
         end do

      end if

      !********************************************************************************************************************!
      !* Print the legality matrix: modes satisfying (or not) the YBJ criterion N lambda_m / f lambda_hor < YBJ_criterion *!
      !********************************************************************************************************************!

      if(mype==0) then 
         open (unit=3334,file='legality.dat',action="write",status="replace")

         do m = 1,n3
            do ikx = 0,ktx

               if(m == 1 .or. ikx <= sqrt(YBJ_criterion*Bu)*eigen_values(m)) then    !The mode is legal for the chosen YBJ_criterion 
                  write(unit=3334,fmt=3334) 1.
               else                                                                  !The mode is not legal
                  write(unit=3334,fmt=3334) 0.
               end if

               !Superbly unefficient, but just for test
               if(m == 2) then
                  if(ikx <= sqrt(YBJ_criterion*Bu)*eigen_values(m)) then 
                     write(*,*) "Mode kh = ",ikx," is legal for the first baroclinic mode"
                  else
                     write(*,*) "Mode kh = ",ikx," is not legal for the first baroclinic mode"
                  end if
               end if

            end do
            write(unit=3334,fmt=*) '           '
         end do
         
3334     format(1x,F10.3,1x)
         close (unit=3334)

      end if


      !*********************************************************************!
      !* Print the initial condition for waves and its (vertical) spectrum *!
      !*********************************************************************!
      if(mype==0) then


         ARm = 0.
         BRm = 0.
         CRm = 0.

         !Regenerate the initial condition (the cosh one)

         delta_a = 200.
         delta_b = 100.
         delta_c = 50.

         xi_a = H_scale/delta_a
         xi_b = H_scale/delta_b
         xi_c = H_scale/delta_c

         cons = 2./sqrt(twopi/2.)

         ave_A = 0.
         ave_B = 0.
         ave_C = 0.

         do iz=1,n3

            z = zas(iz)

            ARz(iz) = cons*xi_a*exp(-xi_a*(z-twopi)**2)
            BRz(iz) = cons*xi_b*exp(-xi_b*(z-twopi)**2)
            CRz(iz) = cons*xi_c*exp(-xi_c*(z-twopi)**2)

            ave_A=ave_A + ARz(iz)
            ave_B=ave_B + BRz(iz)
            ave_C=ave_C + CRz(iz)

         end do
         
         !Remove the horizontal average from B
         ARz=ARz-ave_A/n3
         BRz=BRz-ave_B/n3
         CRz=CRz-ave_C/n3

         !Print the initial condition
         open (unit=6969,file='init.dat',action="write",status="replace")
         do iz=1,n3
            write(unit=6969,fmt=3335) zas(iz)*H_scale,ARz(iz),BRz(iz),CRz(iz)
         enddo
         
3335     format(1x,E12.5,1x,E12.5,1x,E12.5,1x,E12.5,1x)

         close (unit=6969)


         !Calculate the vertical coefficients (for both A and B)
         do m=1,n3
            do iz=1,n3

               ARm(m) = ARm(m) + ARz(iz)*eigen_vectors(iz,m)
               BRm(m) = BRm(m) + BRz(iz)*eigen_vectors(iz,m)
               CRm(m) = CRm(m) + CRz(iz)*eigen_vectors(iz,m)

            end do
         end do
         
         
         !Print the(dimensional WKE spectrum using B and A
         open (unit=696969,file='init_spec.dat',action="write",status="replace")
         do m=1,n3
            write(unit=696969,fmt=3336) m-1 ,   ARm(m)*ARm(m)/(2.*n3),   BRm(m)*BRm(m)/(2.*n3),   CRm(m)*CRm(m)/(2.*n3)  
         enddo
         
3336     format(1x,I3,1x,E12.5,1x,E12.5,1x,E12.5,1x)
         close (unit=696969)
      end if
         
    end SUBROUTINE compute_eigen



END MODULE elliptic
