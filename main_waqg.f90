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

  !**** C = Az and is decomposed into real and imag parts (ex.: C = CR + iCI) even though in Fourier-space both CRk and CIk are complex
  double complex,   dimension(iktx,ikty,n3h0) :: CRk, CIk

  !**** n = nonlinear advection term J(psi,B) **** r = refractive term ~ B*vort
  double complex,   dimension(iktx,ikty,n3h0) :: nBRk, nBIk, rBRk, rBIk
  double precision, dimension(n1d,n2d,n3h0)   :: nBRr, nBIr, rBRr, rBIr

  !**** qw, the wave-averaged feedback onto QGPV ****!
  double complex,   dimension(iktx,ikty,n3h0) :: qwk
  double precision, dimension(n1d,n2d,n3h0)   :: qwr

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

  double complex,   dimension(iktx,ikty,n3h0) :: FtRk, FtIk        !YBJ Forcing term (or any other term on the right-hand side of the YBJ equation, e.g. mean-flow advection of LA) Must be set to zero if no forcing present)
  double precision, dimension(n1d,n2d,n3h0)   :: FtRr, FtIr             

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

  equivalence(FtRr,FtRk)
  equivalence(FtIr,FtIk)

  double precision, dimension(n1d,n2d) :: array2dr
  double complex,   dimension(iktx,ikty) :: array2di

  double precision, dimension(n3)   :: fr_even,fk_even

  equivalence(fr_even,fk_even)
  equivalence(fr_odd ,fk_odd )
  equivalence(array2dr,array2di)

  !Integrating factor (includes both dissipation and mean-flow advection in the Eady problem (must be complex in the latter case). One for the flow, one for the waves (w)
  double complex :: int_factor,int_factor_w

  !Rotational part of u for slice...                                                                                      
  double complex, dimension(iktx,ikty,n3h1) :: u_rot
  double precision, dimension(n1d,n2d,n3h1) :: u_rotr

  equivalence(u_rotr,u_rot)

  !Test the potential energy balance by directly diagnosing d/dt B.
  double complex,   dimension(iktx,ikty,n3h0) :: dBRk, dBIk
  double precision, dimension(n1d,n2d,n3h0)   :: dBRr, dBIr

  equivalence(dBRr,dBRk)
  equivalence(dBIr,dBIk)

  !Barotropic version of the streamfunction
  double complex,   dimension(iktx,ikty,n3h1) :: psik_bt        !pressure, and rhs of pressure equation!
  double precision, dimension(n1d,n2d,n3h1)   :: psir_bt

  equivalence(psir_bt,psik_bt)

  !********************** Initializing... *******************************!

  call initialize_mpi
  call init_files
  call initialize_fftw(array2dr,array2di,fr_even,fk_even)
  call init_arrays
  call init_base_state
  if(mype==0)  call validate_run

  if(restart==1) then
     call read_restart(psik)                                                                                                                                                            
     if(npe > 1) call generate_halo_q(psik)                                                                                                                                             
  else
     call init_eady(psik,psir)
  end if

  call init_q(qk,psik)
  call compute_velo(uk,vk,wk,bk,psik)
  if(npe > 1) call generate_halo(uk,vk,wk,bk)
  if(npe > 1) call generate_halo_q(qk) 
 
 psi_old = psik 
     qok = qk 

  !**Initialize a storm**!
  call generate_fields_stag(wr,n3h2,ARr,n3h0,BRr,n3h0) 
  call fft_r2c(ARr,ARk,n3h0)
  call fft_r2c(BRr,BRk,n3h0)
  AIk = (0.D0,0.D0)
  BIk = (0.D0,0.D0)
  CRk = (0.D0,0.D0)
  CIk = (0.D0,0.D0)
  !----------------------!

  


 !Initial diagnostics!
 !*******************!

 !Compute war/wak if desired                                                                                                                                                     
 if(out_omega==1)  then
    call omega_eqn_rhs(rhs,rhsr,psik)
    call mpitranspose(rhs,iktx,ikty,n3h0,qt,n3,iktyp)
    call omega_equation(wak,qt)
 end if

 if(out_etot ==1) call diag_zentrum(uk,vk,wk,bk,wak,psik,u_rot)

 if(out_we   ==1) call wave_energy(BRk,BIk,CRk,CIk)

 do id_field=1,nfields                                            
    if(out_slice ==1) call slices(uk,vk,bk,psik,qk,ur,vr,br,psir,qr,id_field)
 end do

 do id_field=1,nfieldsw                                            
    if(out_slicew ==1) call slices_waves(BRk,BIk,BRr,BIr,CRk,CIk,id_field)
 end do

 do iz=1,num_spec
    if(out_hspecw ==1) call hspec_waves(BRk,BIk,CRk,CIk,iz)
 end do

 if(dump==1) call dump_restart(psik)


 !************************************************************************!
 !*** 1st time timestep using the projection method with Forward Euler ***!
 !************************************************************************!
 
 time=delt
 if(itermax>0) then
 iter=1

 
 if(no_waves == 1) then
    call convol_q(nqk,nqr,uk,vk,qk,ur,vr,qr)
    nBRk=(0.D0,0.D0)
    nBIk=(0.D0,0.D0)
    rBRk = (0.D0,0.D0)
    rBIk = (0.D0,0.D0)
    ARk = (0.D0,0.D0)
    AIk = (0.D0,0.D0)
 else if(no_waves == 0 .and. barotropize == 1) then
    call convol_q(nqk,nqr,uk,vk,qk,ur,vr,qr)

    !Compute the u and v from the barotropized streamfunction (to be used in convol_waves)                                                                                          
    call barotropize_psi(psik,psik_bt)
    call compute_velo(uk,vk,wk,bk,psik_bt)
    if(npe > 1) call generate_halo(uk,vk,wk,bk)

    if(zero_aveB==1) call sumB(BRk,BIk)                                               !Resets the vertical sum of B to zero                                                         
    call convol_waves(nBRk,nBIk,nBRr,nBIr,uk,vk,BRk,BIk,ur,vr,BRr,BIr)
    call refraction_waqg(rBRk,rBIk,rBRr,rBIr,BRk,BIk,psik_bt,BRr,BIr,psir_bt)
 else
    if(zero_aveB==1) call sumB(BRk,BIk)                                               !Resets the vertical sum of B to zero
    call convol_waqg(nqk,nBRk,nBIk,nqr,nBRr,nBIr,uk,vk,qk,BRk,BIk,ur,vr,qr,BRr,BIr)
    call refraction_waqg(rBRk,rBIk,rBRr,rBIr,BRk,BIk,psik,BRr,BIr,psir)
 end if

 if(linear==1) then
    nqk=(0.D0,0.D0)
   nBRk=(0.D0,0.D0)
   nBIk=(0.D0,0.D0)
 end if

 if(no_dispersion==1) then
    ARk=(0.D0,0.D0)
    AIk=(0.D0,0.D0)
 end if

 if(no_refraction==1) then
    rBRk = (0.D0,0.D0)
    rBIk = (0.D0,0.D0)
 end if

 if(passive_scalar==1) then
    ARk = (0.D0,0.D0)
    AIk = (0.D0,0.D0)
   rBRk = (0.D0,0.D0)
   rBIk = (0.D0,0.D0)
end if

  qok = qk
 BRok = BRk
 BIok = BIk

 !Compute q^1 and B^1 with Forward Euler  
 do izh0=1,n3h0
    izh1=izh0+1
    do iky=1,ikty
       ky = kya(iky)
       do ikx=1,iktx
          kx = kxa(ikx)
          kh2=kx*kx+ky*ky

          if(eady==1 .and. expeady == 0 .and. barotropize == 0) then      !Integrating factor includes the terms Uqx and ULAx
             int_factor   = delt* (i*kx*zash0(izh0)                    +  nuh1 *((1.*kx)**(2.*ilap1 ) + (1.*ky)**(2.*ilap1 )) + nuh2 *((1.*kx)**(2.*ilap2 ) + (1.*ky)**(2.*ilap2 ))  )
             int_factor_w = delt* (i*kx*zash0(izh0)                    +  nuh1w*((1.*kx)**(2.*ilap1w) + (1.*ky)**(2.*ilap1w)) + nuh2w*((1.*kx)**(2.*ilap2w) + (1.*ky)**(2.*ilap2w))  )
          elseif(eady==1 .and. expeady == 1 .and. barotropize == 0) then      !Integrating factor includes the terms U exp() qx and U exp() LAx
             int_factor   = delt* (i*kx*exp(N2_scale*(zash0(izh0)-z0)) +  nuh1 *((1.*kx)**(2.*ilap1 ) + (1.*ky)**(2.*ilap1 )) + nuh2 *((1.*kx)**(2.*ilap2 ) + (1.*ky)**(2.*ilap2 ))  )
             int_factor_w = delt* (i*kx*exp(N2_scale*(zash0(izh0)-z0)) +  nuh1w*((1.*kx)**(2.*ilap1w) + (1.*ky)**(2.*ilap1w)) + nuh2w*((1.*kx)**(2.*ilap2w) + (1.*ky)**(2.*ilap2w))  )
          else if(eady==1 .and. expeady == 0 .and. barotropize == 1) then      !Integrating factor includes the term Uqx but NOT the U LAx term
             int_factor   = delt* (i*kx*zash0(izh0)                    +  nuh1 *((1.*kx)**(2.*ilap1 ) + (1.*ky)**(2.*ilap1 )) + nuh2 *((1.*kx)**(2.*ilap2 ) + (1.*ky)**(2.*ilap2 ))  )
             int_factor_w = delt* (                                       nuh1w*((1.*kx)**(2.*ilap1w) + (1.*ky)**(2.*ilap1w)) + nuh2w*((1.*kx)**(2.*ilap2w) + (1.*ky)**(2.*ilap2w))  )
          elseif(eady==1 .and. expeady == 1 .and. barotropize == 1) then      !Integrating factor includes the terms U exp() qx but not U exp() LAx
             int_factor   = delt* (i*kx*exp(N2_scale*(zash0(izh0)-z0)) +  nuh1 *((1.*kx)**(2.*ilap1 ) + (1.*ky)**(2.*ilap1 )) + nuh2 *((1.*kx)**(2.*ilap2 ) + (1.*ky)**(2.*ilap2 ))  )
             int_factor_w = delt* (                                       nuh1w*((1.*kx)**(2.*ilap1w) + (1.*ky)**(2.*ilap1w)) + nuh2w*((1.*kx)**(2.*ilap2w) + (1.*ky)**(2.*ilap2w))  )
          else
             int_factor   = delt* (                                       nuh1 *((1.*kx)**(2.*ilap1 ) + (1.*ky)**(2.*ilap1 )) + nuh2 *((1.*kx)**(2.*ilap2 ) + (1.*ky)**(2.*ilap2 ))  )
             int_factor_w = delt* (                                       nuh1w*((1.*kx)**(2.*ilap1w) + (1.*ky)**(2.*ilap1w)) + nuh2w*((1.*kx)**(2.*ilap2w) + (1.*ky)**(2.*ilap2w))  )
          end if

          if (L(ikx,iky).eq.1) then
             qk(ikx,iky,izh1) = (  qok(ikx,iky,izh1) - delt* nqk(ikx,iky,izh0) )*exp(-int_factor)
            BRk(ikx,iky,izh0) = ( BRok(ikx,iky,izh0) - delt*nBRk(ikx,iky,izh0)  - delt*(0.5/(Bu*Ro))*kh2*AIk(ikx,iky,izh0) + delt*0.5*rBIk(ikx,iky,izh0) )*exp(-int_factor_w)
            BIk(ikx,iky,izh0) = ( BIok(ikx,iky,izh0) - delt*nBIk(ikx,iky,izh0)  + delt*(0.5/(Bu*Ro))*kh2*ARk(ikx,iky,izh0) - delt*0.5*rBRk(ikx,iky,izh0) )*exp(-int_factor_w)

            if(eady == 1 .and. eady_bnd == 1) then  !Add bottom and top boundary terms
               !Bottom boundary: add vQy and the Ekman term                
               if(mype == 0 .and. izh0 == 1) then
                  if(expeady==1) then    !Expeady => extra H/h factor in the temperature term                                                                                        
                     qk(ikx,iky,izh1) = qk(ikx,iky,izh1) + delt*(1./dz)*(i*kx*N2_scale*Bu*psik(ikx,iky,izh1) + (1.*kh2)*Ek*psi_old(ikx,iky,izh1) )*exp(-int_factor)
                  else
                     qk(ikx,iky,izh1) = qk(ikx,iky,izh1) + delt*(1./dz)*(i*kx*         Bu*psik(ikx,iky,izh1) + (1.*kh2)*Ek*psi_old(ikx,iky,izh1) )*exp(-int_factor)
                  end if
               end if

               !Top Boundary: add vQy            
               if(mype == (npe-1) .and. izh0 == n3h0) then
                  if(expeady==1) then    !Expeady => extra H/h factor in the temperature term                        
                     qk(ikx,iky,izh1) = qk(ikx,iky,izh1) - delt*(1./dz)*(i*kx*N2_scale*Bu)*psik(ikx,iky,izh1)*exp(-int_factor)
                  else
                     qk(ikx,iky,izh1) = qk(ikx,iky,izh1) - delt*(1./dz)*(i*kx*         Bu)*psik(ikx,iky,izh1)*exp(-int_factor)
                  end if
               end if
            end if

          else
             qk(ikx,iky,izh1) = (0.D0,0.D0)
            BRk(ikx,iky,izh0) = (0.D0,0.D0)
            BIk(ikx,iky,izh0) = (0.D0,0.D0)
          endif



       enddo
    enddo
 enddo


 !Generate halo for q
 if(npe > 1) call generate_halo_q(qk)

if(fixed_flow==0) then
 ! --- Recover the streamfunction --- !

 !Keep the old version of psi for the stability of the Ekman term
 psi_old = psik 

 if(no_feedback == 1 .or. no_waves == 1) then
    qwk = (0.D0,0.D0)
 else
    call compute_qw(qwk,BRk,BIk,qwr,BRr,BIr)           ! Compute qw
 end if

 do izh0=1,n3h0                                     ! Compute q* = q - qw
    izh1=izh0+1
    do iky=1,ikty
       do ikx=1,iktx
          if (L(ikx,iky).eq.1) then
             qwk(ikx,iky,izh0)=  qk(ikx,iky,izh1) - qwk(ikx,iky,izh0)
          endif
       enddo
    enddo
 enddo

 call mpitranspose(qwk,iktx,ikty,n3h0,qt,n3,iktyp)  !Transpose q*                                                                                      
 call psi_solver(psik,qt)                           !Solve the QGPV equation L(phi)=q*, assuming psi_z = 0 at top/bot (homogeneous problem)                    

 ! ----------------------------------- !
end if

 !Compute the corresponding u,v,w and t (u and v to be used in convol)                                                                                    
 call compute_velo(uk,vk,wk,bk,psik)
 if(npe > 1) call generate_halo(uk,vk,wk,bk)

end if



 !********************************************************************************!
 !*** Subsequent timesteps using the projection method + leapfrog timestepping ***!
 !********************************************************************************!


  do iter=2,itermax
     
     time=iter*delt

     if(no_waves == 1) then
        call convol_q(nqk,nqr,uk,vk,qk,ur,vr,qr)
        nBRk=(0.D0,0.D0)
        nBIk=(0.D0,0.D0)
        rBRk = (0.D0,0.D0)
        rBIk = (0.D0,0.D0)
        ARk = (0.D0,0.D0)
        AIk = (0.D0,0.D0)
     else if(no_waves == 0 .and. barotropize == 1) then
        call convol_q(nqk,nqr,uk,vk,qk,ur,vr,qr)
        
        !Compute the u and v from the barotropized streamfunction (to be used in convol_waves)                                                                         
        call barotropize_psi(psik,psik_bt)
        call compute_velo(uk,vk,wk,bk,psik_bt)
        if(npe > 1) call generate_halo(uk,vk,wk,bk)
        
        if(zero_aveB==1) call sumB(BRk,BIk)                           !Resets the vertical sum of B to zero                                                                         
        call convol_waves(nBRk,nBIk,nBRr,nBIr,uk,vk,BRk,BIk,ur,vr,BRr,BIr)
        call refraction_waqg(rBRk,rBIk,rBRr,rBIr,BRk,BIk,psik_bt,BRr,BIr,psir_bt)
     else
        if(zero_aveB==1) call sumB(BRk,BIk)                           !Resets the vertical sum of B to zero
        call convol_waqg(nqk,nBRk,nBIk,nqr,nBRr,nBIr,uk,vk,qk,BRk,BIk,ur,vr,qr,BRr,BIr)
        call refraction_waqg(rBRk,rBIk,rBRr,rBIr,BRk,BIk,psik,BRr,BIr,psir)
     end if
     
     if(linear==1) then
        nqk=(0.D0,0.D0)
        nBRk=(0.D0,0.D0)
        nBIk=(0.D0,0.D0)
     end if
     
     if(no_dispersion==1) then
        ARk=(0.D0,0.D0)
        AIk=(0.D0,0.D0)
     end if
     
     if(no_refraction==1) then
        rBRk = (0.D0,0.D0)
        rBIk = (0.D0,0.D0)
     end if

     if(passive_scalar==1) then
        ARk = (0.D0,0.D0)
        AIk = (0.D0,0.D0)
       rBRk = (0.D0,0.D0)
       rBIk = (0.D0,0.D0)
     end if

     if(passive_scalar/=1 .and. no_dispersion/=1 .and. no_waves/=1) then     !This was moved to the right place on October 1st, 2018 (see diary on the same day for details)
        ! --- Recover A from B --- !
        if(eady==1 .and. barotropize == 0) then !Include the mean-flow advection in the right-hand side, J(Psi,LA) --> i*k_x*z*LA
           do izh0=1,n3h0
              do iky=1,ikty
                 do ikx=1,iktx
                    kx = kxa(ikx)

                    if(expeady==1) then
                       if (L(ikx,iky).eq.1) then
                          FtRk(ikx,iky,izh0)=i*kx*exp(N2_scale*(zash0(izh0)-z0))*BRk(ikx,iky,izh0)
                          FtIk(ikx,iky,izh0)=i*kx*exp(N2_scale*(zash0(izh0)-z0))*BIk(ikx,iky,izh0)
                       end if
                    else   !Regular linear shear
                       if (L(ikx,iky).eq.1) then
                          FtRk(ikx,iky,izh0)=i*kx*zash0(izh0)*BRk(ikx,iky,izh0)
                          FtIk(ikx,iky,izh0)=i*kx*zash0(izh0)*BIk(ikx,iky,izh0)
                       end if
                    end if

                 end do
              end do
           end do
        else
           FtRk=(0.D0,0.D0)
           FtIk=(0.D0,0.D0)
        end if

        call compute_sigma(sigma,nBRk, nBIk, rBRk, rBIk, FtRk, FtIk)  !Compute the sum of A
        call mpitranspose(BRk,iktx,ikty,n3h0,BRkt,n3,iktyp)           !Transpose BR to iky-parallelized space 
        call mpitranspose(BIk,iktx,ikty,n3h0,BIkt,n3,iktyp)           !Transpose BK to iky-parallelized space 
        call compute_A(ARk,AIK,BRkt,BIkt,CRk,CIK,sigma)               !Compute A!
        
        ! ------------------------ !
     end if




     !Compute q^n+1 and B^n+1 using leap-frog
     do izh0=1,n3h0
        izh1=izh0+1
        do iky=1,ikty
           ky = kya(iky)
           do ikx=1,iktx
              kx = kxa(ikx)
              kh2=kx*kx+ky*ky

              if(eady==1 .and. expeady == 0 .and. barotropize == 0) then      !Integrating factor includes the terms Uqx and ULAx
                 int_factor   = delt* (i*kx*zash0(izh0)                    +  nuh1 *((1.*kx)**(2.*ilap1 ) + (1.*ky)**(2.*ilap1 )) + nuh2 *((1.*kx)**(2.*ilap2 ) + (1.*ky)**(2.*ilap2 ))  )
                 int_factor_w = delt* (i*kx*zash0(izh0)                    +  nuh1w*((1.*kx)**(2.*ilap1w) + (1.*ky)**(2.*ilap1w)) + nuh2w*((1.*kx)**(2.*ilap2w) + (1.*ky)**(2.*ilap2w))  )
              elseif(eady==1 .and. expeady == 1 .and. barotropize == 0) then      !Integrating factor includes the terms U exp() qx and U exp() LAx
                 int_factor   = delt* (i*kx*exp(N2_scale*(zash0(izh0)-z0)) +  nuh1 *((1.*kx)**(2.*ilap1 ) + (1.*ky)**(2.*ilap1 )) + nuh2 *((1.*kx)**(2.*ilap2 ) + (1.*ky)**(2.*ilap2 ))  )
                 int_factor_w = delt* (i*kx*exp(N2_scale*(zash0(izh0)-z0)) +  nuh1w*((1.*kx)**(2.*ilap1w) + (1.*ky)**(2.*ilap1w)) + nuh2w*((1.*kx)**(2.*ilap2w) + (1.*ky)**(2.*ilap2w))  )
              else if(eady==1 .and. expeady == 0 .and. barotropize == 1) then      !Integrating factor includes the term Uqx but NOT the U LAx term
                 int_factor   = delt* (i*kx*zash0(izh0)                    +  nuh1 *((1.*kx)**(2.*ilap1 ) + (1.*ky)**(2.*ilap1 )) + nuh2 *((1.*kx)**(2.*ilap2 ) + (1.*ky)**(2.*ilap2 ))  )
                 int_factor_w = delt* (                                       nuh1w*((1.*kx)**(2.*ilap1w) + (1.*ky)**(2.*ilap1w)) + nuh2w*((1.*kx)**(2.*ilap2w) + (1.*ky)**(2.*ilap2w))  )
              elseif(eady==1 .and. expeady == 1 .and. barotropize == 1) then      !Integrating factor includes the terms U exp() qx but not U exp() LAx
                 int_factor   = delt* (i*kx*exp(N2_scale*(zash0(izh0)-z0)) +  nuh1 *((1.*kx)**(2.*ilap1 ) + (1.*ky)**(2.*ilap1 )) + nuh2 *((1.*kx)**(2.*ilap2 ) + (1.*ky)**(2.*ilap2 ))  )
                 int_factor_w = delt* (                                       nuh1w*((1.*kx)**(2.*ilap1w) + (1.*ky)**(2.*ilap1w)) + nuh2w*((1.*kx)**(2.*ilap2w) + (1.*ky)**(2.*ilap2w))  )
              else
                 int_factor   = delt* (                                       nuh1 *((1.*kx)**(2.*ilap1 ) + (1.*ky)**(2.*ilap1 )) + nuh2 *((1.*kx)**(2.*ilap2 ) + (1.*ky)**(2.*ilap2 ))  )
                 int_factor_w = delt* (                                       nuh1w*((1.*kx)**(2.*ilap1w) + (1.*ky)**(2.*ilap1w)) + nuh2w*((1.*kx)**(2.*ilap2w) + (1.*ky)**(2.*ilap2w))  )
              end if


              if (L(ikx,iky).eq.1) then
                 qtempk(ikx,iky,izh1) =  qok(ikx,iky,izh1)*exp(-2*int_factor) - 2*delt*nqk(ikx,iky,izh0)*exp(-int_factor) 
                BRtempk(ikx,iky,izh0) = BRok(ikx,iky,izh0)*exp(-2*int_factor_w) - 2*delt*(nBRk(ikx,iky,izh0) + (0.5/(Bu*Ro))*kh2*AIk(ikx,iky,izh0) - 0.5*rBIk(ikx,iky,izh0) )*exp(-int_factor_w)
                BItempk(ikx,iky,izh0) = BIok(ikx,iky,izh0)*exp(-2*int_factor_w) - 2*delt*(nBIk(ikx,iky,izh0) - (0.5/(Bu*Ro))*kh2*ARk(ikx,iky,izh0) + 0.5*rBRk(ikx,iky,izh0) )*exp(-int_factor_w)
                
                if(eady == 1 .and. eady_bnd == 1) then  !Add bottom and top boundary terms      
                   !Bottom boundary: add vQy and the Ekman term                                                                                                             
                   if(mype == 0 .and. izh0 == 1) then
                      if(expeady==1) then    !Expeady => extra H/h factor in the temperature term
                         qtempk(ikx,iky,izh1) = qtempk(ikx,iky,izh1) + 2.*delt*(1./dz)*i*kx*N2_scale*Bu*psik(ikx,iky,izh1)*exp(-int_factor) + 2.*delt*(1./dz)*(1.*kh2)*Ek*psi_old(ikx,iky,izh1)*exp(-2*int_factor)
                      else
                         qtempk(ikx,iky,izh1) = qtempk(ikx,iky,izh1) + 2.*delt*(1./dz)*i*kx*         Bu*psik(ikx,iky,izh1)*exp(-int_factor) + 2.*delt*(1./dz)*(1.*kh2)*Ek*psi_old(ikx,iky,izh1)*exp(-2*int_factor)
                      end if
                   end if

                   !Top Boundary: add vQy                                                                                                                                        
                   if(mype == (npe-1) .and. izh0 == n3h0) then
                      if(expeady==1) then    !Expeady => extra H/h factor in the temperature term 
                         qtempk(ikx,iky,izh1) = qtempk(ikx,iky,izh1) - 2.*delt*(1./dz)*(i*kx*N2_scale*Bu)*psik(ikx,iky,izh1)*exp(-int_factor)
                      else
                         qtempk(ikx,iky,izh1) = qtempk(ikx,iky,izh1) - 2.*delt*(1./dz)*(i*kx*         Bu)*psik(ikx,iky,izh1)*exp(-int_factor)
                      end if
                   end if
                end if

              else
                 qtempk(ikx,iky,izh1) = (0.D0,0.D0)
                BRtempk(ikx,iky,izh0) = (0.D0,0.D0)
                BItempk(ikx,iky,izh0) = (0.D0,0.D0)
              endif

           enddo
        enddo
     enddo

     !Compute d/dt B at the same time that the right-hand side is computed (prior to filtering)
     dBRk = (BRtempk - BRok)/(2.*delt)
     dBIk = (BItempk - BIok)/(2.*delt)
     !Switch signs for the advection due to the base-state, because it is interpreted as Forcing in the right-hand side in the conversion subroutines.
     if(eady==1) then
        FtRk = -FtRk
        FtIk = -FtIk
     end if
     if(out_conv ==1 .and. mod(iter,freq_conv )==0)  then
        call wke_conversion(BRk, BIk, BRr, BIr, FtRk, FtIk, FtRr, FtIr, dBRk, dBIk, dBRr, dBIr)
        call we_conversion(ARk, AIk, BRk, BIk, CRk, CIk,  dBRk, dBIk, nBRk, nBIk, rBRk, rBIk, FtRk, FtIk, dBRr, dBIr, nBRr, nBIr, rBRr, rBIr, FtRr, FtIr)
     end if

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

 !Generate halo for q
 if(npe > 1) call generate_halo_q(qk)
 if(npe > 1) call generate_halo_q(qok)
 

if(fixed_flow==0) then
 ! --- Recover the streamfunction --- !                                                                                                                   

 !Keep the old version of psi for the stability of the Ekman term
 psi_old = psik 

 if(no_feedback == 1 .or. no_waves == 1) then
    qwk = (0.D0,0.D0)
 else
    call compute_qw(qwk,BRk,BIk,qwr,BRr,BIr)           ! Compute qw                                                                                                                    
 end if

 do izh0=1,n3h0                                     ! Compute q* = q - qw                                                                                 
    izh1=izh0+1
    do iky=1,ikty
       do ikx=1,iktx
          if (L(ikx,iky).eq.1) then
             qwk(ikx,iky,izh0)=  qk(ikx,iky,izh1) - qwk(ikx,iky,izh0)
          endif
       enddo
    enddo
 enddo

 call mpitranspose(qwk,iktx,ikty,n3h0,qt,n3,iktyp)  !Transpose rhs -> ft                                                                            
 call psi_solver(psik,qt)                           !Solve the pressure equation laplacian(phi)=f                                                              

 ! ----------------------------------- !  
end if

 !Compute the corresponding u,v,w and t 
 call compute_velo(uk,vk,wk,bk,psik)
 if(npe > 1) call generate_halo(uk,vk,wk,bk) 


 !*** Diagnostics ***!
 !-------------------!

 !Compute w if desired
 if(out_omega==1 .and. (mod(iter,freq_omega) ==0))  then
    call omega_eqn_rhs(rhs,rhsr,psik)
    call mpitranspose(rhs,iktx,ikty,n3h0,qt,n3,iktyp)
    call omega_equation(wak,qt)
    if(npe > 1) call generate_halo_q(wak)
 end if
 
if(out_etot ==1 .and. mod(iter,freq_etot )==0) call diag_zentrum(uk,vk,wk,bk,wak,psik,u_rot)

do id_field=1,nfields
   if(out_slice ==1 .and. mod(iter,freq_slice)==0 .and. count_slice(id_field)<max_slices) call slices(uk,vk,bk,psik,qk,ur,vr,br,psir,qr,id_field)
end do


do id_field=1,nfieldsw
   if(out_slicew ==1 .and. mod(iter,freq_slicew)==0 .and. count_slicew(id_field)<max_slices) call slices_waves(BRk,BIk,BRr,BIr,CRk,CIk,id_field)
end do

do iz=1,num_spec
   if(out_hspecw ==1  .and. mod(iter,freq_hspecw)==0 ) call hspec_waves(BRk,BIk,CRk,CIk,iz)
end do


if(dump==1 .and. mod(iter,freq_dump)==0) call dump_restart(psik)

if(out_we ==1   .and. mod(iter,freq_we   )==0)  call wave_energy(BRk,BIk,CRk,CIk)
 

if(time>maxtime) EXIT
end do !End loop

!************ Terminating processes **********************!                                                                                                                         

  call kill_fftw                                                                                                                                              
  call kill_mpi                                                                                                                                  
 
END PROGRAM main
