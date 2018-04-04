PROGRAM main

  USE parameters
  USE mpi
  USE fft
  USE init
  USE derivatives
  USE elliptic
  USE diagnostics
  USE files

  !********************** Declaring variables *****************************!

  double precision, dimension(n1d,n2d,n3h2)   :: ur,vr,wr,br                     !Velocity and potential temperature fields (r-space)
  double complex,   dimension(iktx,ikty,n3h2) :: uk,vk,wk,bk                     !Velocity and potential temperature fields (k-space)

  double precision, dimension(n1d,n2d,n3h1)   :: war       !1st order vertical velocity as computed in QG (w^1)
  double complex,   dimension(iktx,ikty,n3h1) :: wak

  double precision, dimension(n1d,n2d,n3h1)   :: qr                  
  double complex,   dimension(iktx,ikty,n3h1) :: qk        
  double complex,   dimension(iktx,ikty,n3h1) :: qok          !Old q 
  double complex,   dimension(iktx,ikty,n3h1) :: qtempk        

  !**** B = LA, and both A and B are decomposed into their real and imag parts (ex.: A = AR + iAI)
  double complex,   dimension(iktx,ikty,n3h0) :: BRk, BIk, ARk, AIk
  double precision, dimension(n1d,n2d,n3h0)   :: BRr, BIr, ARr, AIr
  double complex,   dimension(iktx,ikty,n3h0) :: BRok, BIok            !B at the old time step
  double complex,   dimension(iktx,ikty,n3h0) :: BRtempk, BItempk      !B before filering

  !**** n = nonlinear advection term J(psi,B) **** r = refractive term ~ B*vort
  double complex,   dimension(iktx,ikty,n3h0) :: nBRk, nBIk, rBRk, rBIk
  double precision, dimension(n1d,n2d,n3h0)   :: nBRr, nBIr, rBRr, rBIr

  !**** qw, the wave-averaged feedback onto QGPV ****!
  double complex,   dimension(iktx,ikty,n3h0) :: qwk
  double precision, dimension(n1d,n2d,n3h0)   :: qwr

  double complex,   dimension(iktx,ikty,n3h0) :: dqk         !dissipation
  double complex,   dimension(iktx,ikty,n3h1) :: psik        !pressure, and rhs of pressure equation!
  double precision, dimension(n1d,n2d,n3h1)   :: psir
  double complex,   dimension(iktx,ikty,n3h1) :: psi_old     !For computing w...

  double complex,   dimension(iktx,ikty,n3h0) :: rhs         !RHS of elliptic equation (n3h0 version of q at n+1)
  double precision, dimension(n1d,n2d,n3h0)   :: rhsr

  double complex, dimension(iktx,n3, iktyp) :: qt            !Transposed (ky-parallelization) right-hand side   
  double complex, dimension(iktx,n3, iktyp) :: BRkt          !Transposed (ky-parallelization) BRk (this array can most likely be recycled)
  double complex, dimension(iktx,n3, iktyp) :: BIkt          !Transposed (ky-parallelization) BIk (this array can most likely be recycled)

  double precision, dimension(n1d,n2d,n3h0)   :: nqr                  
  double complex,   dimension(iktx,ikty,n3h0) :: nqk        

  double complex, dimension(iktx,ikty,2) :: sigma    !Vertial integral of A(kx,ky), 1=real part, 2=imag part


  equivalence(ur,uk)
  equivalence(vr,vk)
  equivalence(wr,wk)
  equivalence(br,bk)
  equivalence(war,wak)

  equivalence(rhsr,rhs)

  equivalence(psir,psik)
  equivalence(qr,qk)
  equivalence(qwr,qwk)
  equivalence(nqr,nqk)

  equivalence(BRr,BRk)
  equivalence(BIr,BIk)
  equivalence(ARr,ARk)
  equivalence(AIr,AIk)

  equivalence(nBRr,nBRk)
  equivalence(nBIr,nBIk)
  equivalence(rBRr,rBRk)
  equivalence(rBIr,rBIk)

  double precision, dimension(n1d,n2d) :: array2dr
  double complex,   dimension(iktx,ikty) :: array2di

  double precision, dimension(n3)   :: fr_even,fk_even
  double precision, dimension(n3-1) :: fr_odd ,fk_odd

  equivalence(fr_even,fk_even)
  equivalence(fr_odd ,fk_odd )
  equivalence(array2dr,array2di)

  !For implicit dissipation
  double precision :: diss             ! nu_H * kH**(2*ilap) delt

  !Rotational part of u for slice...                                                                                                                                                                                                         
  double complex, dimension(iktx,ikty,n3h1) :: u_rot
  double precision, dimension(n1d,n2d,n3h1) :: u_rotr

  equivalence(u_rotr,u_rot)

  !For the comprehensive test only.
  double complex,   dimension(iktx,ikty,n3h0) :: FbRk, FpRk, FaRk, FtRk    !Forcing terms, real part
  double precision, dimension(n1d,n2d,n3h0)   :: FbRr, FpRr, FaRr

  double complex,   dimension(iktx,ikty,n3h0) :: FbIk, FpIk, FaIk, FtIk    !Forcing terms, imaginary part
  double precision, dimension(n1d,n2d,n3h0)   :: FbIr, FpIr, FaIr

  double precision, dimension(n1d,n2d,n3h0)   :: BRtr, BItr, ARtr, AItr    !True solution

  double complex,   dimension(iktx,ikty,n3h1) :: psitk       !Exact solution for psi 
  double precision, dimension(n1d,n2d,n3h1)   :: psitr   

  double precision ::   error_r,L1_local_r,L2_local_r,Li_local_r,L1_global_r,L2_global_r,Li_global_r
  double precision ::   error_i,L1_local_i,L2_local_i,Li_local_i,L1_global_i,L2_global_i,Li_global_i

  character(len = 32) :: fnamer,fnamei                !future file names                                                                                          

  equivalence(FbRr,FbRk)
  equivalence(FpRr,FpRk)
  equivalence(FaRr,FaRk)

  equivalence(FbIr,FbIk)
  equivalence(FpIr,FpIk)
  equivalence(FaIr,FaIk)

  equivalence(psitr,psitk)

  !********************** Initializing... *******************************!


  iter=0

  call initialize_mpi
  call init_files
  call initialize_fftw(array2dr,array2di,fr_even,fk_even,fr_odd,fk_odd)
  call init_arrays
  call init_base_state


  !Initialize the test!
  !*******************!

  !First generate the exact solution and the associated forcing term
  call generate_fields_stag(ARtr,n3h0,BRtr,n3h0,psitr,n3h1) 
!  call generate_fields_stag2(AItr,n3h0,BItr,n3h0,war,n3h1) 
!  call generate_fields_stag3(FbRr,n3h0,FpRr,n3h0,FaRr,n3h0) 
!  call generate_fields_stag4(FbIr,n3h0,FpIr,n3h0,FaIr,n3h0) 

  psir = psitr

  BRr = BRtr
  BIr = BItr
  ARr = ARtr
  AIr = AItr
 
  !Move to Fourier space.
  call fft_r2c(psir,psik,n3h1)

  call fft_r2c(ARr,ARk,n3h0)
  call fft_r2c(AIr,AIk,n3h0)

  call fft_r2c(BRr,BRk,n3h0)
  call fft_r2c(BIr,BIk,n3h0)

  call fft_r2c(FbRr,FbRk,n3h0)
  call fft_r2c(FpRr,FpRk,n3h0)
  call fft_r2c(FaRr,FaRk,n3h0)

  call fft_r2c(FbIr,FbIk,n3h0)
  call fft_r2c(FpIr,FpIk,n3h0)
  call fft_r2c(FaIr,FaIk,n3h0)

  !Set QG terms to 0.
   qk = (0.D0,0.D0)
  dqk=(0.D0,0.D0)

  !Initialize the other fields
  call compute_velo(uk,vk,wk,bk,psik)
  call generate_halo(uk,vk,wk,bk)

  !Initialize error file!
  if(mype==0) then
     write (fnamer, "(A3,I1,A1,I1)") "br_",vres,"_",tres
     write (fnamei, "(A3,I1,A1,I1)") "bi_",vres,"_",tres
     open (unit=154673,file=fnamer,action="write",status="replace")
     open (unit=154674,file=fnamei,action="write",status="replace")
  end if


 !************************************************************************!
 !*** 1st time timestep using the projection method with Forward Euler ***!
 !************************************************************************!
 
 time=delt
 if(itermax>0) then
! if(mype==0) write(*,*) "First time step"

 iter=1
 
 call convol_waqg(nqk,nBRk,nBIk,nqr,nBRr,nBIr,uk,vk,qk,BRk,BIk,ur,vr,qr,BRr,BIr)
 call refraction_waqg(rBRk,rBIk,rBRr,rBIr,BRk,BIk,psik,BRr,BIr,psir)

  qok = qk
 BRok = BRk
 BIok = BIk


 !Compute the forcing!
 if(forcing==1) then
    do izh0=1,n3h0
       do iky=1,ikty
          do ikx=1,iktx
             if (L(ikx,iky).eq.1) then
                FtRk(ikx,iky,izh0) = FbRk(ikx,iky,izh0)*sin(a_t*(time-delt)) + FpRk(ikx,iky,izh0)*cos(a_t*(time-delt)) + FaRk(ikx,iky,izh0)*cos(a_t*(time-delt))
                FtIk(ikx,iky,izh0) = FbIk(ikx,iky,izh0)*sin(a_t*(time-delt)) + FpIk(ikx,iky,izh0)*cos(a_t*(time-delt)) + FaIk(ikx,iky,izh0)*cos(a_t*(time-delt))
             else
                FtRk(ikx,iky,izh0) = (0.D0,0.D0)
                FtIk(ikx,iky,izh0) = (0.D0,0.D0)
             endif
          enddo
       enddo
    enddo
 else
    FtRk(ikx,iky,izh0) = (0.D0,0.D0)
    FtIk(ikx,iky,izh0) = (0.D0,0.D0)
 end if

 !Compute q^1 and B^1 with Forward Euler  
 do izh0=1,n3h0
    izh1=izh0+1
    do iky=1,ikty
       ky = kya(iky)
       do ikx=1,iktx
          kx = kxa(ikx)
          kh2=kx*kx+ky*ky
          diss = nuh*delt*(1.*kh2)**(1.*ilap)              !This does not work !!!!! diss = nuh*(kh2**ilap)*delt 
          if (L(ikx,iky).eq.1) then
             qk(ikx,iky,izh1) = (  qok(ikx,iky,izh1) - delt* nqk(ikx,iky,izh0)  + delt*dqk(ikx,iky,izh0) )*exp(-diss)
   BRk(ikx,iky,izh0) = ( BRok(ikx,iky,izh0) + delt*FtRk(ikx,iky,izh0)  - delt*nBRk(ikx,iky,izh0)  - delt*(0.5/(Bu*Ro))*kh2*AIk(ikx,iky,izh0) + delt*0.5*rBIk(ikx,iky,izh0) )*exp(-diss)
   BIk(ikx,iky,izh0) = ( BIok(ikx,iky,izh0) + delt*FtIk(ikx,iky,izh0)  - delt*nBIk(ikx,iky,izh0)  + delt*(0.5/(Bu*Ro))*kh2*ARk(ikx,iky,izh0) - delt*0.5*rBRk(ikx,iky,izh0) )*exp(-diss)
          else
             qk(ikx,iky,izh1) = (0.D0,0.D0)
            BRk(ikx,iky,izh0) = (0.D0,0.D0)
            BIk(ikx,iky,izh0) = (0.D0,0.D0)
          endif
       enddo
    enddo
 enddo


 ! --- Recover A from B --- !

 call compute_sigma(sigma,nBRk, nBIk, rBRk, rBIk)              !Compute the sum of A
 call mpitranspose(BRk,iktx,ikty,n3h0,BRkt,n3,iktyp)           !Transpose BR to iky-parallelized space 
 call mpitranspose(BIk,iktx,ikty,n3h0,BIkt,n3,iktyp)           !Transpose BK to iky-parallelized space 
 call compute_A(ARk,AIK,BRkt,BIkt,sigma)                       !Compute A!

 ! ------------------------ !

end if



 !********************************************************************************!
 !*** Subsequent timesteps using the projection method + leapfrog timestepping ***!
 !********************************************************************************!


!  if(mype==0) write(*,*) "Subsequent time steps"
  do iter=2,itermax

!     if(mype==0)  cputime=etime(tarray1)
     
     time=iter*delt

     call convol_waqg(nqk,nBRk,nBIk,nqr,nBRr,nBIr,uk,vk,qk,BRk,BIk,ur,vr,qr,BRr,BIr)
     call refraction_waqg(rBRk,rBIk,rBRr,rBIr,BRk,BIk,psik,BRr,BIr,psir)

     !Compute the forcing!
     if(forcing==1) then
        do izh0=1,n3h0
           do iky=1,ikty
              do ikx=1,iktx
                 if (L(ikx,iky).eq.1) then
                    FtRk(ikx,iky,izh0) = FbRk(ikx,iky,izh0)*sin(a_t*(time-delt)) + FpRk(ikx,iky,izh0)*cos(a_t*(time-delt)) + FaRk(ikx,iky,izh0)*cos(a_t*(time-delt))
                    FtIk(ikx,iky,izh0) = FbIk(ikx,iky,izh0)*sin(a_t*(time-delt)) + FpIk(ikx,iky,izh0)*cos(a_t*(time-delt)) + FaIk(ikx,iky,izh0)*cos(a_t*(time-delt))
                 else
                    FtRk(ikx,iky,izh0) = (0.D0,0.D0)
                    FtIk(ikx,iky,izh0) = (0.D0,0.D0)
                 endif
              enddo
           enddo
        enddo
     else
        FtRk(ikx,iky,izh0) = (0.D0,0.D0)
        FtIk(ikx,iky,izh0) = (0.D0,0.D0)
     end if
     
     !Compute q^n+1 and B^n+1 using leap-frog
     do izh0=1,n3h0
        izh1=izh0+1
        do iky=1,ikty
           ky = kya(iky)
           do ikx=1,iktx
              kx = kxa(ikx)
              kh2=kx*kx+ky*ky
              diss = nuh*delt*(1.*kh2)**(1.*ilap)
              if (L(ikx,iky).eq.1) then
                 qtempk(ikx,iky,izh1) =  qok(ikx,iky,izh1)*exp(-2*diss) - 2*delt*nqk(ikx,iky,izh0)*exp(-diss)  + 2*delt*dqk(ikx,iky,izh0)*exp(-2*diss)
  BRtempk(ikx,iky,izh0) = BRok(ikx,iky,izh0)*exp(-2*diss) - 2*delt*(nBRk(ikx,iky,izh0) - FtRk(ikx,iky,izh0) + (0.5/(Bu*Ro))*kh2*AIk(ikx,iky,izh0) - 0.5*rBIk(ikx,iky,izh0) )*exp(-diss)
  BItempk(ikx,iky,izh0) = BIok(ikx,iky,izh0)*exp(-2*diss) - 2*delt*(nBIk(ikx,iky,izh0) - FtIk(ikx,iky,izh0) - (0.5/(Bu*Ro))*kh2*ARk(ikx,iky,izh0) + 0.5*rBRk(ikx,iky,izh0) )*exp(-diss)
              else
                 qtempk(ikx,iky,izh1) = (0.D0,0.D0)
                BRtempk(ikx,iky,izh0) = (0.D0,0.D0)
                BItempk(ikx,iky,izh0) = (0.D0,0.D0)
              endif
           enddo
        enddo
     enddo


     !Apply Robert-Asselin filter to damp the leap-frog computational mode
     do izh0=1,n3h0
        izh1=izh0+1
        do iky=1,ikty
           do ikx=1,iktx
              if (L(ikx,iky).eq.1) then
                 qok(ikx,iky,izh1) =  qk(ikx,iky,izh1) + gamma * (  qok(ikx,iky,izh1) - 2 *  qk(ikx,iky,izh1) +  qtempk(ikx,iky,izh1) )
                BRok(ikx,iky,izh0) = BRk(ikx,iky,izh0) + gamma * ( BRok(ikx,iky,izh0) - 2 * BRk(ikx,iky,izh0) + BRtempk(ikx,iky,izh0) )
                BIok(ikx,iky,izh0) = BIk(ikx,iky,izh0) + gamma * ( BIok(ikx,iky,izh0) - 2 * BIk(ikx,iky,izh0) + BItempk(ikx,iky,izh0) )
              else
                 qok(ikx,iky,izh1) = (0.D0,0.D0)
                BRok(ikx,iky,izh0) = (0.D0,0.D0)
                BIok(ikx,iky,izh0) = (0.D0,0.D0)
              endif
           enddo
        enddo
     enddo

!Overwrite the new field uk with u^{n+1} 
 qk =  qtempk
BRk = BRtempk
BIk = BItempk


 ! --- Recover A from B --- !                                                                                                                                 

 call compute_sigma(sigma,nBRk, nBIk, rBRk, rBIk)              !Compute the sum of A                                                                                    
 call mpitranspose(BRk,iktx,ikty,n3h0,BRkt,n3,iktyp)           !Transpose BR to iky-parallelized space                                                                   
 call mpitranspose(BIk,iktx,ikty,n3h0,BIkt,n3,iktyp)           !Transpose BK to iky-parallelized space                                                                  
 call compute_A(ARk,AIK,BRkt,BIkt,sigma)                       !Compute A!                                                                                  

 ! ------------------------ !       



 !*** Compute the error ***!
 !-------------------------!

 if(out_slice ==1 .and. mod(iter,freq_slice)==0) then
    
    !Print slices of the solution!
    call fft_c2r(BRk,BRr,n3h0)
    call fft_c2r(BIk,BIr,n3h0)
!    do id_field=8,nfields                                            
       !call slices(uk,vk,wk,bk,wak,u_rot,ur,vr,wr,br,war,u_rotr,psir,psitr*cos(a_t*time),id_field)
       !      call slices(uk,vk,wk,bk,wak,u_rot,ur,vr,wr,br,war,u_rotr,qr,qtr*cos(a_t*time),id_field)
!    end do
    
    !Compute the error on psi!
    error_r   =0.
    L1_local_r=0.
    L2_local_r=0.
    Li_local_r=0.
    L1_global_r=0.
    L2_global_r=0.
    Li_global_r=0.

    error_i   =0.
    L1_local_i=0.
    L2_local_i=0.
    Li_local_i=0.
    L1_global_i=0.
    L2_global_i=0.
    Li_global_i=0.
    
    do izh0=1,n3h0
       izh1=izh0+1
       do ix=1,n1
          do iy=1,n2
             
             error_r    = ABS( BRtr(ix,iy,izh0)*cos(a_t*time) - BRr(ix,iy,izh0) )
             L1_local_r = L1_local_r + error_r
             L2_local_r = L2_local_r + error_r**2
             if(error_r > Li_local_r) Li_local_r = error_r

             error_i    = ABS( BItr(ix,iy,izh0)*cos(a_t*time) - BIr(ix,iy,izh0) )
             L1_local_i = L1_local_i + error_i
             L2_local_i = L2_local_i + error_i**2
             if(error_i > Li_local_i) Li_local_i = error_i
             
          end do
       end do
    end do
    
    call mpi_reduce(L1_local_r,L1_global_r, 1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    call mpi_reduce(L2_local_r,L2_global_r, 1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    call mpi_reduce(Li_local_r,Li_global_r, 1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD,ierror)
   
    call mpi_reduce(L1_local_i,L1_global_i, 1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    call mpi_reduce(L2_local_i,L2_global_i, 1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    call mpi_reduce(Li_local_i,Li_global_i, 1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD,ierror)

    if (mype==0) then
       write(154673,"(E12.5,E12.5,E12.5,E12.5)") time,L1_global_r/(n1*n2*n3),sqrt(L2_global_r/(n1*n2*n3)),Li_global_r
       write(154674,"(E12.5,E12.5,E12.5,E12.5)") time,L1_global_i/(n1*n2*n3),sqrt(L2_global_i/(n1*n2*n3)),Li_global_i
    end if
    
    call fft_r2c(BRr,BRk,n3h0)
    call fft_r2c(BIr,BIk,n3h0)    
 end if
 
 if(time>maxtime) EXIT
end do !End loop

 
!************ Terminating processes **********************!                                                                                                                         

  call kill_fftw                                                                                                                                              
  call kill_mpi                                                                                                                                  
 
END PROGRAM main
