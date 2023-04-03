1) modifications for using gate soundings

a)include the line in the begin of the module:
 USE module_gate
 
b) Include the lines below for setting the parameters for GATE 
   check if the names of temp and mix ratio are correct (t2d,q2d,tn, QO)
!- begin: for GATE soundings-------------------------------------------
    
     if(use_gate) then  
       do i=its,itf
         do k=kts,ktf
           p2d (i,k) = ppres(jl,k)
           t2d (i,k) = ptemp(jl,k)+273.15
           q2d (i,k) = pq(jl,k)/1000.
           us  (i,k) = pu(jl,k)
           vs  (i,k) = pv(jl,k)
           omeg(i,k,:)=pvervel(jl,k) 
           phil(i,k) = pgeo(jl,k)*g   !geo
           po  (i,k) = p2d (i,k)
           zo  (i,k) = phil(i,k)/g    !meters
           rhoi(i,k) =1.e2*p2d(i,k)/(rgas*t2d(i,k))	 
          enddo
          ter11(i)  = phil(i,1)/g  ! phil is given in g*h.
          psur(i)   = p2d(i,1)                
	  tsur(i)   = t2d(i,1)
          gswi(i,1) = 0.0 
          kpbli(i)  = 4                      
          zws  (i)  = 1.0 ! wstar 
          do k=kts,ktf
           tn  (i,k) = t2d(i,k) + dt *(zadvt(jl,k)+zqr(jl,k))/86400.
           QO  (i,k) = q2d(i,k) + dt * zadvq(jl,k)
          enddo
	  pbl (i)  = zo(i,kpbli(i))
       enddo
     endif     
!- end:   for GATE soundings-------------------------------------------


c) After the lines:
 
 !--- DEEP CONVECTION
  IF(IDEEP_GF ==1 )then
  
 include the if command
 
    "if(.not. use_gate) then ! GATE" 
 
 and before the line "call cup_gf", include the line
    "endif" 
 
 
d) Just before the line    "END SUBROUTINE CUP_gf" include the lines"
 
 !- begin: for GATE soundings-------------------------------------------
  if(use_gate) then
    do i=its,itf
     if(ierr(i).eq.0) then
      !- 2-d section
      do k=kts,ktop(i)

       nvar=1 
       call setgradsvar(i,k,nvar,xmb(i)*zuo(i,k),"mup",'m flux up (kg/s/m^2)')

       nvar=2 
       call setgradsvar(i,k,nvar,-edto(i)*xmb(i)*zdo(i,k),"mdn" ,'m flux down (kg/s/m^2)')

       nvar=3 
       call setgradsvar(i,k,nvar,dellah(i,k),"delh" ,'dellah')

       nvar=4 
       call setgradsvar(i,k,nvar,dellaq(i,k)*86400.*xl/cp*xmb(i), "dellq" ,'dellaq')

       nvar=5 
       call setgradsvar(i,k,nvar,dellaqc(i,k)*86400.*xl/cp*xmb(i),"dellqc" ,'dellaqc')

       !outvars(i,k,7)=up_massentro(i,k)
       !outvars(i,k,8)=up_massdetro(i,k)
      enddo
	   
 endif
!- end  : for GATE soundings-------------------------------------------

e) include the routine below in the module_cu_gf90   
  
  SUBROUTINE SETGRADSVAR(i,k,nvar,f,name1,name2)
     implicit none
     integer, intent(in) :: nvar,i,k
     real, intent(in) :: f
     character*(*) :: name1,name2
     
     cupout(nvar)%varp(i,k)= f
     cupout(nvar)%varn(1)=name1 
     cupout(nvar)%varn(2)=name2 

   END SUBROUTINE SETGRADSVAR
