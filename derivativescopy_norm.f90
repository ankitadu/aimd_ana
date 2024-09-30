program derivatives_norm  
      integer :: nFrame_ho, iFrame_ho
      integer :: nFrame_hh, iFrame_hh
      real(kind(0.0d0)) :: CN(1:1000), Peak(1:1000), nF
      real(kind(0.0d0)) :: ENER(1:1000)
      real(kind(0.0d0)) :: roo, dU_deps_oo, S_dU_deps_oo, dU_drmin_oo, S_dU_drmin_oo 
      real(kind(0.0d0)) :: eps_oo, eps_ow, rmin_oo    
      real(kind(0.0d0)) :: roh, dU_deps_oh, S_dU_deps_oh, dU_drmin_oh, S_dU_drmin_oh 
      real(kind(0.0d0)) :: eps_oh, eps_hw, rmin_oh  
      real(kind(0.0d0)) :: rho, dU_deps_ho, S_dU_deps_ho, dU_drmin_ho, S_dU_drmin_ho 
      real(kind(0.0d0)) :: eps_ho, rmin_ho   
      real(kind(0.0d0)) :: rhh, dU_deps_hh, S_dU_deps_hh, dU_drmin_hh, S_dU_drmin_hh 
      real(kind(0.0d0)) :: eps_hh, rmin_hh, Int_ener
      real(kind(0.0d0)) :: a_CN, a_CNh, a_Peak, a_Peakh, a_Int, a_Inth
      real(kind(0.0d0)) :: CN_i, Peak_i, CN_j, Peak_j, Int_ener_i, Int_ener_j
      real(kind(0.0d0)) :: S_dU_deps_owat, CN_dU_deps_owat, SCN_dU_deps_owat, a_dU_deps_owat, a_CN_dU_deps_owat, dCN_depso
      real(kind(0.0d0)) :: Int_dU_deps_owat, SInt_dU_desp_owat, Int_dU_drmin_owat, SInt_dU_drmin_owat
      real(kind(0.0d0)) :: a_Int_dU_desp_owat, a_Int_dU_drmin_owat, dInt_deps, dInt_drmin
      real(kind(0.0d0)) :: Int_dU_depsh_owat, SInt_dU_desph_owat, Int_dU_drminh_owat, SInt_dU_drminh_owat
      real(kind(0.0d0)) :: a_Int_dU_desph_owat, a_Int_dU_drminh_owat, dInt_depsh, dInt_drminh
      real(kind(0.0d0)) :: S_dU_drmin_owat, CN_dU_drmin_owat, SCN_dU_drmin_owat, a_dU_drmin_owat, a_CN_dU_drmin_owat, dCN_drmin
      real(kind(0.0d0)) :: Peak_dU_deps_owat, SPeak_dU_deps_owat, a_Peak_dU_deps_owat, dPeak_deps
      real(kind(0.0d0)) :: Peak_dU_drmin_owat, SPeak_dU_drmin_owat, a_Peak_dU_drmin_owat, dPeak_drmin
      real(kind(0.0d0)) :: CN_k, Peak_k, CN_l, Peak_l, Int_ener_k, Int_ener_l
      real(kind(0.0d0)) :: S_dU_depsh_owat, CN_dU_depsh_owat, SCN_dU_depsh_owat, a_dU_depsh_owat, a_CN_dU_depsh_owat 
      real(kind(0.0d0)) :: S_dU_drminh_owat, CN_dU_drminh_owat, SCN_dU_drminh_owat, a_dU_drminh_owat, a_CN_dU_drminh_owat
      real(kind(0.0d0)) :: Peak_dU_depsh_owat, SPeak_dU_depsh_owat, a_Peak_dU_depsh_owat, dPeak_depsh
      real(kind(0.0d0)) :: Peak_dU_drminh_owat, SPeak_dU_drminh_owat, a_Peak_dU_drminh_owat, dPeak_drminh
      real(kind(0.0d0)) :: dCN_drminh, dCN_depsh
      real(kind(0.0d0)) :: SS_dU_deps_owat, SS_dU_depsh_owat, SS_dU_drmin_owat, SS_dU_drminh_owat
      real(kind(0.0d0)) :: err_CN, err_Peak, S_err, chi_sq, dchi_sq, alpha, new_eps, new_rmin, new_epsh, new_rminh
      real(kind(0.0d0)) :: dchi_sq_depso, dchi_sq_drmino, dchi_sq_depsh, dchi_sq_drminh, alpha_epso, alpha_drmino
      real(kind(0.0d0)) :: epso_norm, rmino_norm, epsh_norm, rminh_norm
      real(kind(0.0d0)) :: epso_unorm, rmino_unorm, epsh_unorm, rminh_unorm
      character*20 :: dumch                        !dumch read in can be 20 characters long                               
      
     
      

!get energy and orientation data 

      n= 600
      open(unit = 13, file = 'cn_peak_ener.dat')
      do i = 1, n
         read(13,*) nF, CN(i), Peak(i), ENER(i)
      end do
      i = 0

!calculate dE/dlam for OH-OW atom pairs------------------------------------------ 

      open(unit = 10, file = 'natomdistoow.txt')
      open(unit = 11, file = 'gradientoow.txt')

      read(10,*) !Skip the header                     !similarly 10 references the input file
      nFrame_oo = 1
      nWat = 0
      S_dU_deps_oo = 0d0
      S_dU_drmin_oo = 0d0


12    read(10,*,end=13) dumch ! Read the first argument of the line. At the end of the file goto line 13


      if (dumch .eq. "Parameter:") then ! One frame is finished. Write the results of the frame in the output file.
         i = i + 1
         Int_ener = EXP(-(ENER(i)/0.594))
         write(11,*) nFrame_oo, CN(i), Peak(i), Int_ener, S_dU_deps_oo, S_dU_drmin_oo 

         nWat = 0
         S_dU_deps_oo = 0d0
         S_dU_drmin_oo = 0d0
         nFrame_oo = nFrame_oo + 1



      else ! The line is a data line
         BACKSPACE 10 !Go back one line and read data of the line again.
         read(10,*) roo       !deta structure of the excel file I sent

         nWat = nWat + 1
         eps_oo = sqrt(0.024)
         eps_ow = sqrt(0.1521)
         rmin_oo = (3.965 + 3.5364)
         dU_deps_oo = (1/(2*eps_oo))*eps_ow* ((((rmin_oo)/(2*roo))**12) - (2*((rmin_oo)/(2*roo))**6))
         dU_drmin_oo = ((3*eps_oo*eps_ow)*((rmin_oo**11)/(1024*(roo**12)))) - ((6*eps_oo*eps_ow)* ((rmin_oo**5)/(32*(roo**6))))
         S_dU_deps_oo = S_dU_deps_oo + dU_deps_oo
         S_dU_drmin_oo = S_dU_drmin_oo + dU_drmin_oo
         



         endif
         goto 12

13       i = i + 1
         write(11,*) nFrame_oo, CN(i), Peak(i), Int_ener, S_dU_deps_oo, S_dU_drmin_oo 

         close(10)
         close(11)

!calculate dE/dlam for OH-HW atom pairs----------------------------------- 

      open(unit = 16, file = 'natomdistohw.txt')
      open(unit = 17, file = 'gradientohw.txt')

      
      read(16,*) !Skip the header                     !similarly 10 references the input file
      nFrame_oh = 1
      nWat = 0
      S_dU_deps_oh = 0d0
      S_dU_drmin_oh = 0d0
      


      i = 0

14    read(16,*,end=15) dumch ! Read the first argument of the line. At the end of the file goto line 13

      if (dumch .eq. "Parameter:") then ! One frame is finished. Write the results of the frame in the output file.
         i = i + 1
         Int_ener = EXP(-(ENER(i)/0.594))

         write(17,*) nFrame_oh, CN(i), Peak(i), Int_ener, S_dU_deps_oh, S_dU_drmin_oh 
         nWat = 0
         S_dU_deps_oh = 0d0
         S_dU_drmin_oh = 0d0
         nFrame_oh = nFrame_oh + 1


      else ! The line is a data line
         BACKSPACE 16 !Go back one line and read data of the line again.
         read(16,*) roh       !deta structure of the excel file I sent

         nWat = nWat + 1
         eps_oh = sqrt(0.024)
         eps_hw = sqrt(0.046)
         rmin_oh = (3.965 + 0.449)
         dU_deps_oh = (1/(2*eps_oh))*eps_hw* ((((rmin_oh)/(2*roh))**12) - (2*((rmin_oh)/(2*roh))**6))
         dU_drmin_oh = ((3*eps_oh*eps_hw)*((rmin_oh**11)/(1024*(roh**12)))) - ((6*eps_oh*eps_hw)* ((rmin_oh**5)/(32*(roh**6))))
         S_dU_deps_oh = S_dU_deps_oh + dU_deps_oh
         S_dU_drmin_oh = S_dU_drmin_oh + dU_drmin_oh


         endif
         goto 14
                                
15       i = i + 1
         write(17,*) nFrame_oh, CN(i), Peak(i), Int_ener, S_dU_deps_oh, S_dU_drmin_oh 

         close(16)
         close(17)

!calculate dE/dlam for HO-OW atom pairs-------------------------------------------- 

      open(unit = 18, file = 'natomdistHOW.txt')
      open(unit = 19, file = 'gradientHOW.txt')

      
      read(18,*) !Skip the header                     !similarly 10 references the input file
      nFrame_ho = 1
      nWat = 0
      S_dU_deps_ho = 0d0
      S_dU_drmin_ho = 0d0
      


      i = 0

16    read(18,*,end=17) dumch ! Read the first argument of the line. At the end of the file goto line 13

      if (dumch .eq. "Parameter:") then ! One frame is finished. Write the results of the frame in the output file.
         i = i + 1
         Int_ener = EXP(-(ENER(i)/0.594))

         write(19,*) nFrame_ho, CN(i), Peak(i), Int_ener, S_dU_deps_ho, S_dU_drmin_ho 
         nWat = 0
         S_dU_deps_ho = 0d0
         S_dU_drmin_ho = 0d0
         nFrame_ho = nFrame_ho + 1


      else ! The line is a data line
         BACKSPACE 18 !Go back one line and read data of the line again.
         read(18,*) rho       !deta structure of the excel file I sent

         nWat = nWat + 1
         eps_ho = sqrt(0.046)
         eps_ow = sqrt(0.1521)
         rmin_ho = (3.5364 + 0.449)
         dU_deps_ho = (1/(2*eps_ho))*eps_ow* ((((rmin_ho)/(2*rho))**12) - (2*((rmin_ho)/(2*rho))**6))
         dU_drmin_ho = ((3*eps_ho*eps_ow)*((rmin_ho**11)/(1024*(rho**12)))) - ((6*eps_ho*eps_ow)* ((rmin_ho**5)/(32*(rho**6))))
         S_dU_deps_ho = S_dU_deps_ho + dU_deps_ho
         S_dU_drmin_ho = S_dU_drmin_ho + dU_drmin_ho


         endif
         goto 16

17       i = i + 1
         write(19,*) nFrame_ho, CN(i), Peak(i), Int_ener, S_dU_deps_ho, S_dU_drmin_ho 

         close(18)
         close(19)
    

!calculate dE/dlam for HO-HW atom pairs-------------------------------------------- 

      open(unit = 20, file = 'natomdistHHW.txt')
      open(unit = 21, file = 'gradientHHW.txt')

      
      read(20,*) !Skip the header                     !similarly 10 references the input file
      nFrame_hh = 1
      nWat = 0
      S_dU_deps_hh = 0d0
      S_dU_drmin_hh = 0d0
      


      i = 0

18    read(20,*,end=19) dumch ! Read the first argument of the line. At the end of the file goto line 13

      if (dumch .eq. "Parameter:") then ! One frame is finished. Write the results of the frame in the output file.
         i = i + 1
         Int_ener = EXP(-(ENER(i)/0.594))

         write(21,*) nFrame_hh, CN(i), Peak(i), Int_ener, S_dU_deps_hh, S_dU_drmin_hh 
         nWat = 0
         S_dU_deps_hh = 0d0
         S_dU_drmin_hh = 0d0
         nFrame_hh = nFrame_hh + 1


      else ! The line is a data line
         BACKSPACE 20 !Go back one line and read data of the line again.
         read(20,*) rhh       !deta structure of the excel file I sent

         nWat = nWat + 1
         eps_hh = sqrt(0.046)
         eps_hw = sqrt(0.046)
         rmin_hh = (0.449 + 0.449)
         dU_deps_hh = (1/(2*eps_hh))*eps_hw* ((((rmin_hh)/(2*rhh))**12) - (2*((rmin_hh)/(2*rhh))**6))
         dU_drmin_hh = ((3*eps_hh*eps_hw)*((rmin_hh**11)/(1024*(rhh**12)))) - ((6*eps_hh*eps_hw)* ((rmin_hh**5)/(32*(rhh**6))))
         S_dU_deps_hh = S_dU_deps_hh + dU_deps_hh
         S_dU_drmin_hh = S_dU_drmin_hh + dU_drmin_hh


         endif
         goto 18

19       i = i + 1
         write(21,*) nFrame_oh, CN(i), Peak(i), Int_ener, S_dU_deps_hh, S_dU_drmin_hh 

         close(20)
         close(21)

!All average calculations-----------------------------------------------------------------
! dCN/deps (for O-wat)------------------------------------   
      
         open(unit = 22, file = 'gradientoow.txt')
         open(unit = 23, file = 'gradientohw.txt')
         a_CN = 0d0
         a_Peak = 0d0
         do j = 1, 599
            read(22,*) nFrame_oo, CN_i, Peak_i, Int_ener_i, S_dU_deps_oo, S_dU_drmin_oo 
            read(23,*) nFrame_oh, CN_j, Peak_j, Int_ener_j, S_dU_deps_oh, S_dU_drmin_oh 
                        
                        
                             a_CN = a_CN + CN_i
                             a_Peak = a_Peak + Peak_i
                             a_Int = a_Int + Int_ener_i

                             S_dU_deps_owat = S_dU_deps_oo + S_dU_deps_oh
                             SS_dU_deps_owat = SS_dU_deps_owat + S_dU_deps_owat

                             CN_dU_deps_owat = CN_i * S_dU_deps_owat
                             Peak_dU_deps_owat = Peak_i * S_dU_deps_owat
                             Int_dU_deps_owat = Int_ener_i * S_dU_deps_owat


                             SCN_dU_deps_owat = SCN_dU_deps_owat + CN_dU_deps_owat
                             SPeak_dU_deps_owat = SPeak_dU_deps_owat + Peak_dU_deps_owat
                             SInt_dU_desp_owat = SInt_dU_desp_owat + Int_dU_deps_owat

                             S_dU_drmin_owat = S_dU_drmin_oo  + S_dU_drmin_oh
                             SS_dU_drmin_owat = SS_dU_drmin_owat + S_dU_drmin_owat

                             CN_dU_drmin_owat = CN_i * S_dU_drmin_owat
                             Peak_dU_drmin_owat = Peak_i * S_dU_drmin_owat
                             Int_dU_drmin_owat = Int_ener_i * S_dU_drmin_owat

                             SCN_dU_drmin_owat = SCN_dU_drmin_owat+ CN_dU_drmin_owat
                             SPeak_dU_drmin_owat = SPeak_dU_drmin_owat + Peak_dU_drmin_owat
                             SInt_dU_drmin_owat = SInt_dU_drmin_owat + Int_dU_drmin_owat 




                             
                             

         enddo

         a_CN = a_CN/599
         a_Peak = a_Peak/599
         a_Int = a_Int/599

print*, 'avg CN=', a_CN
print*, 'avg_Peak=', a_Peak
print*, 'avg_FE=', a_INT

         a_dU_deps_owat = SS_dU_deps_owat/600
         a_dU_drmin_owat = SS_dU_drmin_owat/600

         a_CN_dU_deps_owat = SCN_dU_deps_owat/600
         a_CN_dU_drmin_owat = SCN_dU_drmin_owat/600

         a_Peak_dU_deps_owat = SPeak_dU_deps_owat/600
         a_Peak_dU_drmin_owat = SPeak_dU_drmin_owat/600

         a_Int_dU_desp_owat = SInt_dU_desp_owat/600
         a_Int_dU_drmin_owat = SInt_dU_drmin_owat/600
         


         dCN_depso = -1.69*(a_CN_dU_deps_owat - (a_CN * a_dU_deps_owat))
         dCN_drmin = -1.69*(a_CN_dU_drmin_owat - (a_CN * a_dU_drmin_owat))
         dPeak_deps = -1.69*(a_Peak_dU_deps_owat - (a_Peak * a_dU_deps_owat))
         dPeak_drmin = -1.69*(a_Peak_dU_drmin_owat - (a_Peak * a_dU_drmin_owat))
         dInt_deps = -0.594*(1/a_Int)*(-1.69*(a_Int_dU_desp_owat - (a_Int * a_dU_deps_owat)))
         dInt_drmin = -0.594*(1/a_Int)*(-1.69*(a_Int_dU_drmin_owat - (a_Int * a_dU_drmin_owat)))


        close(22)
        close(23)


         open(unit = 24, file = 'gradientHOW.txt')
         open(unit = 25, file = 'gradientHHW.txt')
         a_CNh= 0d0
         a_Peakh = 0d0
         do k = 1, 600
            read(24,*) nFrame_ho, CN_k, Peak_k, Int_ener_k, S_dU_deps_ho, S_dU_drmin_ho 
            read(25,*) nFrame_hh, CN_l, Peak_l, Int_ener_l, S_dU_deps_hh, S_dU_drmin_hh 
                        
                        
                             a_CNh = a_CNh + CN_k
                             a_Peakh = a_Peakh + Peak_k
                             a_Inth = a_Inth + Int_ener_k

                             S_dU_depsh_owat = S_dU_deps_ho + S_dU_deps_hh
                             SS_dU_depsh_owat = SS_dU_depsh_owat + S_dU_depsh_owat

                             CN_dU_depsh_owat = CN_k * S_dU_depsh_owat
                             Peak_dU_depsh_owat = Peak_k * S_dU_depsh_owat
                             Int_dU_depsh_owat = Int_ener_k * S_dU_depsh_owat

                             SCN_dU_depsh_owat = SCN_dU_depsh_owat + CN_dU_depsh_owat
                             SPeak_dU_depsh_owat = SPeak_dU_depsh_owat + Peak_dU_depsh_owat
                             SInt_dU_desph_owat = SInt_dU_desph_owat + Int_dU_depsh_owat
    

                             S_dU_drminh_owat = S_dU_drmin_ho  + S_dU_drmin_hh
                             SS_dU_drminh_owat = SS_dU_drminh_owat + S_dU_drminh_owat

                             CN_dU_drminh_owat = CN_k * S_dU_drminh_owat
                             Peak_dU_drminh_owat = Peak_k * S_dU_drminh_owat
                             Int_dU_drminh_owat = Int_ener_k * S_dU_drminh_owat

                             SCN_dU_drminh_owat = SCN_dU_drminh_owat+ CN_dU_drminh_owat
                             SPeak_dU_drminh_owat = SPeak_dU_drminh_owat + Peak_dU_drminh_owat
                             SInt_dU_drminh_owat = SInt_dU_drminh_owat + Int_dU_drminh_owat




                             
                             

         enddo

         a_CNh = a_CNh/600
         a_Peakh = a_Peakh/600

         a_dU_depsh_owat = SS_dU_depsh_owat/600
         a_dU_drminh_owat = SS_dU_drminh_owat/600

         a_CN_dU_depsh_owat = SCN_dU_depsh_owat/600
         a_CN_dU_drminh_owat = SCN_dU_drminh_owat/600

         a_Peak_dU_depsh_owat = SPeak_dU_depsh_owat/600
         a_Peak_dU_drminh_owat = SPeak_dU_drminh_owat/600
         
         a_Int_dU_desph_owat = SInt_dU_desph_owat/600
         a_Int_dU_drminh_owat = SInt_dU_drminh_owat/600


         dCN_depsh = -1.69*(a_CN_dU_depsh_owat - (a_CNh * a_dU_depsh_owat))
         dCN_drminh = -1.69*(a_CN_dU_drminh_owat - (a_CNh * a_dU_drminh_owat))
         dPeak_depsh = -1.69*(a_Peak_dU_depsh_owat - (a_Peakh * a_dU_depsh_owat))
         dPeak_drminh = -1.69*(a_Peak_dU_drminh_owat - (a_Peakh * a_dU_drminh_owat))
         dInt_depsh = -0.594*(1/a_Int)*(-1.69*(a_Int_dU_desph_owat - (a_Int * a_dU_depsh_owat)))
         dInt_drminh = -0.594*(1/a_Int)*(-1.69*(a_Int_dU_drminh_owat - (a_Int * a_dU_drminh_owat)))

        close(22)
        close(23)

!running a quick check 
!print*, 'dCN/dEpso=', dCN_depso
!print*, 'dCN/drmino=', dCN_drmin
!print*, 'dPeak/dEpso=', dPeak_deps
!print*, 'dPeak/drmino=',  dPeak_drmin 
!print*, 'dCN/dEpsh=', dCN_depsh
!print*, 'dCN/drminh=', dCN_drminh
!print*, 'dPeak/dEpsh=', dPeak_depsh
!print*, 'dPeak/drminh=', dPeak_drminh 
!print*, 'dInt/deps =', dInt_deps
!print*, 'dInt/drmin=', dInt_drmin
!print*, 'dInt/depsh =', dInt_depsh
!print*, 'dInt/drminh =', dInt_drminh


!Normalizing---------------------------------------------------

epso_norm = (0.024 - 0.015)/(0.070 - 0.015)
rmino_norm = (3.965 - 3.000)/(4.000 - 3.000)
epsh_norm = (0.046 - 0.040)/(0.050 - 0.040)
rminh_norm = (0.449 - 0.4000)/(3.9000 - 0.4000)

!Applying the steepest descent method---------------------------------------------------------

!Calculate chi^2 and its deriv

err_CN = (a_CN - 4.5)
err_Peak = (a_Peak - 4)
!S_err = (err_CN) + (err_Peak)
!chi_sq = ((a_CN - 4.5)**2) + ((a_Peak- 4)**2)
!dchi_sq = (2*S_err*dCN_depso ) + (2*S_err*dCN_drmin ) + (2*S_err*dPeak_deps ) + (2*S_err*dPeak_drmin ) + (2*S_err*dCN_depsh ) &
!          + (2*S_err*dCN_drminh ) + (2*S_err*dPeak_depsh ) + (2*S_err*dPeak_drminh )

print*, 'error CN=', err_CN
print*, 'error_Peak=', err_Peak

dchi_sq_depso = (2*err_CN*dCN_depso)+(2*err_Peak*dPeak_deps)
dchi_sq_drmino = (2*err_CN*dCN_drmin)+(2*err_Peak*dPeak_drmin)
dchi_sq_depsh = (2*err_CN*dCN_depsh)+(2*err_Peak*dPeak_depsh)
dchi_sq_drminh = (2*err_CN*dCN_drminh)+(2*err_Peak*dPeak_drminh)

print*, 'dchi_sq/depso =', dchi_sq_depso
print*, 'dchi_sq/drmino =', dchi_sq_drmino 
print*, 'dchi_sq/depsh =', dchi_sq_depsh 
print*, 'dchi_sq/drminh =', dchi_sq_drminh

!Calculate alpha for both parameters 

alpha = 0.0002

!compute new parameter for oxygen 

new_eps = epso_norm - dchi_sq_depso*alpha*(0.070 - 0.015)
new_rmin = rmino_norm - dchi_sq_drmino*alpha*(4.000 - 3.000)
new_epsh = epsh_norm - dchi_sq_depsh*alpha*(0.050 - 0.040)
new_rminh = rminh_norm - dchi_sq_drminh*alpha*(3.9000 - 0.4000)


!Un-Normalize---------------------------------------------------

epso_unorm = new_eps * (0.07 - 0.015) + 0.015
rmino_unorm = new_rmin * (4.000 - 3.000) + 3.000
epsh_unorm = new_epsh * (0.050 - 0.040) + 0.040
rminh_unorm = new_rminh * (3.9000 - 0.4000) + 0.4000

!----------------------------------------------------------------

open(unit = 15, file = 'param_update.txt')
write(21,*) epso_unorm, rmino_unorm, epsh_unorm, rminh_unorm
      
end program derivatives_norm
