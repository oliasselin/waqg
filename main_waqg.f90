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

  equivalence(fr_even,fk_even)
  equivalence(fr_odd ,fk_odd )
  equivalence(array2dr,array2di)

  !For implicit dissipation
  double complex :: diss             ! nu_H * kH**(2*ilap) delt

  !Rotational part of u for slice...                                                                                                                                                                                                         
  double complex, dimension(iktx,ikty,n3h1) :: u_rot
  double precision, dimension(n1d,n2d,n3h1) :: u_rotr

  equivalence(u_rotr,u_rot)

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

  !For the stability test!
  call generate_fields_stag(BRr,n3h0,BIr,n3h0,wr,n3h2) 
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
 else
    call convol_waqg(nqk,nBRk,nBIk,nqr,nBRr,nBIr,uk,vk,qk,BRk,BIk,ur,vr,qr,BRr,BIr)
    call refraction_waqg(rBRk,rBIk,rBRr,rBIr,BRk,BIk,psik,BRr,BIr,psir)
 end if

 !Compute dissipation 
 call dissipation_q_nv(dqk,qok)
 
 if(inviscid==1) then
    dqk=(0.D0,0.D0)
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

  qok = qk
 BRok = BRk
 BIok = BIk

 if(passive_scalar==1) then
    ARk = (0.D0,0.D0)
    AIk = (0.D0,0.D0)
   rBRk = (0.D0,0.D0)
   rBIk = (0.D0,0.D0)
end if

 !Compute q^1 and B^1 with Forward Euler  
 do izh0=1,n3h0
    izh1=izh0+1
    do iky=1,ikty
       ky = kya(iky)
       do ikx=1,iktx
          kx = kxa(ikx)
          kh2=kx*kx+ky*ky

          if(eady==1) then
             diss = delt* (i*kx*zash0(izh0) +  nuh*((1.*kx)**(2.*ilap) + (1.*ky)**(2.*ilap))  )     !Integrating factor includes the term Uqx
          else
             diss = nuh*delt*((1.*kx)**(2.*ilap) + (1.*ky)**(2.*ilap))          
          end if

          if (L(ikx,iky).eq.1) then
             qk(ikx,iky,izh1) = (  qok(ikx,iky,izh1) - delt* nqk(ikx,iky,izh0)  + delt*dqk(ikx,iky,izh0) )*exp(-diss)
            BRk(ikx,iky,izh0) = ( BRok(ikx,iky,izh0) - delt*nBRk(ikx,iky,izh0)  - delt*(0.5/(Bu*Ro))*kh2*AIk(ikx,iky,izh0) + delt*0.5*rBIk(ikx,iky,izh0) )*exp(-diss)
            BIk(ikx,iky,izh0) = ( BIok(ikx,iky,izh0) - delt*nBIk(ikx,iky,izh0)  + delt*(0.5/(Bu*Ro))*kh2*ARk(ikx,iky,izh0) - delt*0.5*rBRk(ikx,iky,izh0) )*exp(-diss)

            if(eady == 1 .and. eady_bnd == 1) then  !Add bottom and top boundary terms
               !Bottom boundary: add vQy and the Ekman term                
               if(mype == 0 .and. izh0 == 1) then
                  qk(ikx,iky,izh1) = qk(ikx,iky,izh1) + delt*(1./dz)*(i*kx*Bu*psik(ikx,iky,izh1) + (1.*kh2)*Ek*psi_old(ikx,iky,izh1) )*exp(-diss)
               end if

               !Top Boundary: add vQy            
               if(mype == (npe-1) .and. izh0 == n3h0) then
                  qk(ikx,iky,izh1) = qk(ikx,iky,izh1) - delt*(1./dz)*(i*kx*Bu)*psik(ikx,iky,izh1)*exp(-diss)
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

 if(no_waves == 1) then
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

if(passive_scalar==0 .and. no_waves/=1) then
 ! --- Recover A from B --- !

 if(zero_aveB==1) call sumB(BRk,BIk)                           !Resets the vertical sum of B to zero

 call compute_sigma(sigma,nBRk, nBIk, rBRk, rBIk)              !Compute the sum of A
 call mpitranspose(BRk,iktx,ikty,n3h0,BRkt,n3,iktyp)           !Transpose BR to iky-parallelized space 
 call mpitranspose(BIk,iktx,ikty,n3h0,BIkt,n3,iktyp)           !Transpose BK to iky-parallelized space 
 call compute_A(ARk,AIK,BRkt,BIkt,CRk,CIK,sigma)               !Compute A!

 ! ------------------------ !
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
     else
        call convol_waqg(nqk,nBRk,nBIk,nqr,nBRr,nBIr,uk,vk,qk,BRk,BIk,ur,vr,qr,BRr,BIr)
        call refraction_waqg(rBRk,rBIk,rBRr,rBIr,BRk,BIk,psik,BRr,BIr,psir)
     end if
 
     !Compute dissipation from qok
     call dissipation_q_nv(dqk,qok)

     if(inviscid==1) then
        dqk=(0.D0,0.D0)
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

     !Compute q^n+1 and B^n+1 using leap-frog
     do izh0=1,n3h0
        izh1=izh0+1
        do iky=1,ikty
           ky = kya(iky)
           do ikx=1,iktx
              kx = kxa(ikx)
              kh2=kx*kx+ky*ky

              if(eady==1) then
                 diss = delt* (i*kx*zash0(izh0) +  nuh*((1.*kx)**(2.*ilap) + (1.*ky)**(2.*ilap)) )     !Integrating factor includes the term Uqx
              else
                 diss = nuh*delt*((1.*kx)**(2.*ilap) + (1.*ky)**(2.*ilap))          
              end if

              if (L(ikx,iky).eq.1) then
                 qtempk(ikx,iky,izh1) =  qok(ikx,iky,izh1)*exp(-2*diss) - 2*delt*nqk(ikx,iky,izh0)*exp(-diss)  + 2*delt*dqk(ikx,iky,izh0)*exp(-2*diss)
                BRtempk(ikx,iky,izh0) = BRok(ikx,iky,izh0)*exp(-2*diss) - 2*delt*(nBRk(ikx,iky,izh0) + (0.5/(Bu*Ro))*kh2*AIk(ikx,iky,izh0) - 0.5*rBIk(ikx,iky,izh0) )*exp(-diss)
                BItempk(ikx,iky,izh0) = BIok(ikx,iky,izh0)*exp(-2*diss) - 2*delt*(nBIk(ikx,iky,izh0) - (0.5/(Bu*Ro))*kh2*ARk(ikx,iky,izh0) + 0.5*rBRk(ikx,iky,izh0) )*exp(-diss)
                
                if(eady == 1 .and. eady_bnd == 1) then  !Add bottom and top boundary terms      
                   !Bottom boundary: add vQy and the Ekman term                                                                                                             
                   if(mype == 0 .and. izh0 == 1) then
                      qtempk(ikx,iky,izh1) = qtempk(ikx,iky,izh1) + 2.*delt*(1./dz)*i*kx*Bu*psik(ikx,iky,izh1)*exp(-diss) + 2.*delt*(1./dz)*(1.*kh2)*Ek*psi_old(ikx,iky,izh1)*exp(-2*diss)
                   end if

                   !Top Boundary: add vQy                                                                                                                                        
                   if(mype == (npe-1) .and. izh0 == n3h0) then
                      qtempk(ikx,iky,izh1) = qtempk(ikx,iky,izh1) - 2.*delt*(1./dz)*(i*kx*Bu)*psik(ikx,iky,izh1)*exp(-diss)
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

 if(no_waves == 1) then
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


if(passive_scalar==0 .and. no_waves/=1) then
 ! --- Recover A from B --- !                                                                                                                                 

 if(zero_aveB==1) call sumB(BRk,BIk)                           !Resets the vertical sum of B to zero

 call compute_sigma(sigma,nBRk, nBIk, rBRk, rBIk)              !Compute the sum of A                                                                                    
 call mpitranspose(BRk,iktx,ikty,n3h0,BRkt,n3,iktyp)           !Transpose BR to iky-parallelized space                                                                   
 call mpitranspose(BIk,iktx,ikty,n3h0,BIkt,n3,iktyp)           !Transpose BK to iky-parallelized space                                                                  
 call compute_A(ARk,AIK,BRkt,BIkt,CRk,CIK,sigma)               !Compute A!                                                                                                               

 ! ------------------------ !       
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

if(dump==1 .and. mod(iter,freq_dump)==0) call dump_restart(psik)

if(out_we ==1   .and. mod(iter,freq_we   )==0)  call wave_energy(BRk,BIk,CRk,CIk)
 

if(time>maxtime) EXIT
end do !End loop

!************ Terminating processes **********************!                                                                                                                         

  call kill_fftw                                                                                                                                              
  call kill_mpi                                                                                                                                  
 
END PROGRAM main
