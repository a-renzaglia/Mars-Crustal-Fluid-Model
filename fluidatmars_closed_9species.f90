Program fluidatmars_closed
Implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This code has incorporated continuity and momentum equations to a quasi-1D plasma flow model. Quasi-1D because it takes 3D equations but only looks along magnetic flux tubes.
!! Current setup is for L-shell of 1.1, but should be able to be outfitted for multiple L-shell values. Necessary changes: total length of system, delta(s), s vs lambda input file (#25) for new flux tubes
!! Advances through time using 2-step Lax Wendroff method (explained below)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Inputs/Outputs
!! #25: Inputs s [km] vs lambda [radians]
!! #26: Inputs neutral density [cm^-3] vs altitude [km]
!! #27: Inputs initial values for w1(NGrid),w2(NGrid),w3(NGrid) (also outputs final values to same file) Allows for taking results from run 1 and immediately plug into run 2 as inputs
!! #35: Outputs subroutine results. Allows us to check that they are running correctly
!! #45: Outputs density, velocity, and temperature of system along the flux tubes, at multiple time intervals. The main output from the code that we need.
!! #101: Output for general checks on the code
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Fluid Dynamics
!! w1 = A*n, f1 = u*w1 = A*u*n
!! w2 = f1 = A*u*n
!! f2 = u*w2+CSS0*w1 = A*u*u*n+(kb/m*(T_e+T_i))*A*n
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 2 step Lax Wendroff
!! moves w1, w2 one step forward in time
!! first find half time step forward, half spatial step forward/backward
!! find fhp1/fhm1/fhp2/fhp2 from whp1/whm1/whp2/whm2
!! then find full time step forward, in same location
!!
!! w1m(I)=0.5*(w1(I-1)+w1(I))-0.5*DTDS*(F1(I)-F1(I-1))+0.25*(R1(jj)+R1(j))*DELT
!! w1p(I)=0.5*(w1(I+1)+w1(I))-0.5*DTDS*(F1(I+1)-F1(I))+0.25*(R1(jj)+R1(j))*DELT
!! w1(I)=w1(I)-DTDS*(fhp1(I)-fhm1(I))+R1(jj)*DELT
!!
!! similar equations for w2
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Other notes
!! We increase number of steps on each end of the flux tube (by 10 steps each). Allows us to deal with boundary discontinuities.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Revised for Mars Crustal Fields
!!two magnetic dipoles, 100km below surface, 500km apart
!!Case 1: Max Altitude of field line at 400km
!!Case 2: Max Altitude of field line at 600km
!!Case 3: Max Altitude of field line at 800km
!!!!!!!!!!!!!!!!!!!!!!!!!!!


!! introduce variables used in code
integer n
parameter(n=125)           !magnitude of s vs lambda input file
!real orig_lam(n), toss(n), orig_s(n)
integer i,j,k,NGrid1,jj
integer NIon,NNeut,NGrid,NReact,NSteps
parameter(NIon=9)		!number of ion species
parameter(NNeut=10)		!number of neutral species
parameter(NGrid=125)		!number of gridpoints
parameter(NReact=1)		!number of ion-neutral chemical reactions
parameter(NSteps=1000000)	!number of time steps
integer m,h,mm
parameter(m=500)        !Magnitude of input files (primprod, neutral density, electron temps)
parameter(mm=901)	!Magnitude of ion temps input file
real orig_s(n), x(n), y(n), alt(n), toss(n)
real og_z(m),og_ProdRate(m,NIon),ProdRate(NGrid,NIon)
real or_z(m),or_NeutDen(m,NNeut),NeutDen(NGrid,NNeut)
!real og_Npr(m),og_N2pr(m),og_CO2pr(m),og_COpr(m),og_Opr(m),og_Cpr(m),og_Hepr(m)
!real og_NOpr(m),og_O2pr(m),og_Arpr(m),og_Hpr(m)
!real og_Ndoubpr(m),og_CO2doubpr(m),og_COdoubpr(m),og_Odoubpr(m),og_total_PR(m)
!real Npr(NGrid),O2pr(NGrid),CO2pr(NGrid),COpr(NGrid),Opr(NGrid),Cpr(NGrid)
!real Hepr(NGrid),NOpr(NGrid),O2pr(NGrid),Arpr(NGrid),Hpr(NGrid)
!real Ndoubpr(NGrid),CO2doubpr(NGrid),COdoubpr(NGrid),Odoubpr(NGrid),total_PR(NGrid)
!real orz(m),ornden(m)
!real or_N2(m),or_CO2(m),or_CO(m)
!real or_O(m),or_He(m),or_total_nden(m)
!real N2(NGrid),CO2(NGrid),CO(NGrid),O(NGrid),He(NGrid),total_nden(NGrid)
real DELT,DTDS,TME,den,uu,temp,Rmars,Rdipole
real DNN(NSteps,NGrid,NIon),UUU(NSteps,NGrid,NIon),NElec(NGrid),elecden(NSteps,NGrid)
real w1(NGrid,NIon),w1m(NGrid,NIon),w1p(NGrid,NIon),w1new(NGrid,NIon)
real w2(NGrid,NIon),w2m(NGrid,NIon),w2p(NGrid,NIon),w2new(NGrid,NIon)
real f1(NGrid,NIon),f2(NGrid,NIon),fhm1(NGrid,NIon),fhm2(NGrid,NIon),fhp1(NGrid,NIon),fhp2(NGrid,NIon)
real A(NGrid),GPAR(NGrid),lambda(NGrid),R(NGrid),z(NGrid),s(NGrid),nden(NGrid),Atld(NGrid)
real orig_te(m),orig_alt(m),orig_ti(mm),orig_alt2(mm)
real t_i(NGrid),t_e(NGrid),CSS0(NGrid,NIon),mass(NIon)
real R1(NGrid,NIon),R2(NGrid,NIon)
real prod(NGrid,NIon),loss(NGrid,NIon),grav(NGrid,NIon),coll(NGrid,NIon),epres(NGrid,NIon),mload(NGrid,NIon)
real RHS1,RHS2
real L,SMAX,DELS,A0,DDEN,DVEL,infinity
real lam_deg(NGrid)
real alpha(NGrid,NIon),prod_rate(NGrid,NIon),den_i(NGrid,NIon),uu_i(NGrid,NIon),elecden_i(NGrid)

doubleprecision newtime, time(NSteps)

!!!!!!!!!!!!!!!!!!
!! Notes on variables defined
!! n is length of s vs lambda input file; orig_lam(n)/orig_s(n) are lambda/s values from that file; toss(n) is used to get rid of initial index
!! si is the number of spatial steps in the system; NSteps is the number of time steps in the system
!! Mnden is length of the nuetral density vs altitude file; orz(mnden)/ornden(mnden) are altitude (z)/neutral density values from that file
!! delta(t); dt/ds; time; ion density; ion velocity; ion temperature; radius of Mars (core to surface)
!! time (function of # of steps); output ion density; output ion velocity; output ion temperature
!! establishing wN(NGrid) variables from fluid dynamics (M/P notations stand for minus/plus (for 2-step LW))
!! establishing fN(NGrid) variables from fluid dynamics (hm/hp notations stand for half-minus/half-plus (for 2-step LW))
!! subroutine outputs: area; external accel(centripetal - gravity); lambda; radius; altitude(above cloud tops); distance along flux tube; neutral H2 density
!! subroutine output for electron and neutral temperatures; speed of sound squared for electrons
!! dw3ds is approximation of dw3/ds=(w3m-w3p)/(2*delta(s))
!! R1 is RHS for continuity eqn; R2 is RHS for momentum eqn; R3 is RHS for energy eqn
!! various constants: L-shell value; max length of flux tube; delta(s); initial area of flux tube, damping term (for final calculations of wN terms); infinity term
!!!!!!!!!!!!!!!!!!!!!

!! open necessary files
open(25, file='inputs/sgrid_case1_400km.txt')
open(26, file='inputs/mars_neutraldensity_lowH.txt')
open(27, file='inputs/primprod_transport_lowH.txt')
open(28, file='inputs/electemp.txt')
open(29, file='inputs/iontemp.txt')
open(30, file='init_conditions_time.txt')
open(31, file='init_conditions_w1.txt')
open(32, file='init_conditions_w2.txt')
open(35, file='diagnostics/test_subroutines.txt')
open(36, file='diagnostics/check_neutraldensities.txt')
open(37, file='diagnostics/check_productionrates.txt')
open(38, file='diagnostics/check_alphavalues.txt')
open(40, file='outputs/o2+_densities.txt')
open(41, file='outputs/o+_densities.txt')
open(42, file='outputs/co2+_densities.txt')
open(43, file='outputs/co+_densities.txt')
open(44, file='outputs/no+_densities.txt')
open(45, file='outputs/n2+_densities.txt')
open(46, file='outputs/oh+_densities.txt')
open(47, file='outputs/h+_densities.txt')
open(48, file='outputs/hco+_densities.txt')
open(50, file='outputs/o2+_velocities.txt')
open(51, file='outputs/o+_velocities.txt')
open(52, file='outputs/co2+_velocities.txt')
open(53, file='outputs/co+_velocities.txt')
open(54, file='outputs/no+_velocities.txt')
open(55, file='outputs/n2+_velocities.txt')
open(56, file='outputs/oh+_velocities.txt')
open(57, file='outputs/h+_velocities.txt')
open(58, file='outputs/hco+_velocities.txt')
open(60, file='outputs/e-_densities.txt')
open(101, file='diagnostics/123456.txt')
open(102, file='diagnostics/check_rh_terms_o2+.txt')
open(103, file='diagnostics/check_rh_terms_o+.txt')
open(104, file='diagnostics/check_rh_terms_co2+.txt')
open(105, file='diagnostics/check_rh_terms_co+.txt')
open(106, file='diagnostics/check_rh_terms_no+.txt')
open(107, file='diagnostics/check_rh_terms_n2+.txt')
open(108, file='diagnostics/check_rh_terms_oh+.txt')
open(109, file='diagnostics/check_rh_terms_h+.txt')
open(110, file='diagnostics/check_rh_terms_hco+.txt')
open(111, file='diagnostics/check_rotation.txt')
open(112, file='diagnostics/check_prod_lowalt.txt')
open(113, file='diagnostics/check_prod_midalt.txt')
open(114, file='diagnostics/check_prod_highalt.txt')
open(121, file='diagnostics/check_initialconditions_w1.txt')
open(122, file='diagnostics/check_initialconditions_w2.txt')

!! define constants/parameters of code
NGrid1 = (NGrid)-1
SMAX = 640522.         ![m]
DELT = .010            ![s]
DELS = SMAX/(NGrid1)      ![m]
Rmars = 3389.5         ![km]
A0 = 1.0               ![m^2?]
DTDS = DELT/DELS       ![s/m]
DDEN = 1.0E6
DVEL = 1.0E6

!!!!Neutral List
!1=N2; 2=CO2; 3=CO; 4=O; 5=O2; 6=N; 7=H; 8=NO; 9=H2; 10=H2O
!!!!Ion List
!1=O2+; 2=O+; 3=CO2+; 4=CO+; 5=NO+; 6=N2+; 7=OH+; 8=H+; 9=HCO+

!!!define masses of ions [kg]
mass(1)=5.352390E-26
mass(2)=2.676195E-26
mass(3)=7.359536E-26
mass(4)=4.6833E-26
mass(5)=5.0179E-26
mass(6)=4.6833E-26
mass(7)=2.8435E-26
mass(8)=1.6726E-27
mass(9)=4.8506E-26

!!!read in s file, also gives x location and altitude (y location) at each point
read(25,*)
do i=1,n
	read(25,*) toss(i), orig_s(i), x(i), y(i)
	!convert orig values from km to m, where needed
	s(i) = 1000*orig_s(i)
	alt(i) = y(i)
enddo

!!!read in neutral density file, linear interpolate to fill in between file values
!units are: Altitude: [km]; Neutral Density: [cm^-3]
!read(26,*)
!read(26,*)
!do i=1,m-1
!	read(26,*) or_z(i),(or_NeutDen(i,j),j=1,NNeut)
!enddo
!
!do i=1,m-1
!   do j=1,NGrid
!      do k=1,NNeut
!		if(alt(j)==or_z(i)) then
!			NeutDen(j,k)=or_NeutDen(i,k)
!		else if(alt(j)>or_z(i) .and. alt(j)<or_z(i+1)) then
!			call linearint(or_z(i),or_z(i+1),or_NeutDen(i,k),or_NeutDen(i+1,k),alt(j),NeutDen(j,k))
!		endif
!      enddo
!   enddo
!enddo

!write out equations for neutral densities (if applicable)
do j=1,NGrid
	if (alt(j)<200.0) then
		NeutDen(j,1)=  (alt(j)/1155.5)**(-125.0/13.0)			!N2
		NeutDen(j,2)=  (alt(j)/764.40)**(-125.0/9.00)			!CO2
		NeutDen(j,3)=  (alt(j)/504.17)**(-1000./61.0)			!CO
		NeutDen(j,4)=  (alt(j)/5252.5)**(-1000./181.)			!O
		NeutDen(j,5)=  (alt(j)/1797.7)**(-500.0/67.0)			!O2
		NeutDen(j,6)=  (alt(j)/633.28)**(-200.0/19.0)			!N
		NeutDen(j,7)=  (alt(j)/639800.)**(-100./67.0)			!H
		NeutDen(j,8)=  (alt(j)/1259.1)**(-1000./141.)			!NO
		NeutDen(j,9)= -1.1718*(alt(j))**3.+573.88*(alt(j))**2.&		!H2
			-92171.*alt(j)+5.063E6
		NeutDen(j,10)= (alt(j)/1259.1)**(-1000./141.)			!H2O
	elseif (alt(j)>=200.0) then
		NeutDen(j,1)=  exp((alt(j)-571.92)/(-22.06))			!N2
		NeutDen(j,2)=  exp((alt(j)-461.13)/(-14.04))			!CO2
		NeutDen(j,3)=  exp((alt(j)-530.42)/(-22.06))			!CO
		NeutDen(j,4)=  exp((alt(j)-900.00)/(-38.60))			!O
		NeutDen(j,5)=  exp((alt(j)-517.60)/(-19.30))			!O2
		NeutDen(j,6)=  exp((alt(j)-736.25)/(-44.12))			!N
		NeutDen(j,7)=  exp((alt(j)-7626.6)/(-617.7))			!H
		NeutDen(j,8)=  exp((alt(j)-470.17)/(-20.59))			!NO
		NeutDen(j,9)=  exp((alt(j)-2092.3)/(-154.4))			!H2
		NeutDen(j,10)= exp((alt(j)-470.17)/(-20.59))			!H2O
	endif
enddo

!make nden bottom out at 1E-10(if actual value lower, set to 1E-10)
!***revisit this, maybe just make some total neutral density that has to be above 1E-10***
do i=1,NGrid
   do k=1,NNeut
      if(NeutDen(i,k)<1E-10) then
		NeutDen(i,k)=1E-10
      endif
   enddo
enddo

!!!read in production rate file, linear interpolate to fill in between file values
!units are: Altitude: [km]; Production Rate: [cm^-3/s]
read(27,*)
do i=1,m-1
	read(27,*) og_z(i),(og_ProdRate(i,j),j=1,NIon)
enddo

do i=1,m-1
   do j=1,NGrid
      do k=1,NIon
		if (alt(j)==og_z(i)) then
			ProdRate(j,k)=og_ProdRate(i,k)
		elseif(alt(j)>og_z(i) .and. alt(j)<og_z(i+1)) then
			call linearint(og_z(i),og_z(i+1),og_ProdRate(i,k),og_ProdRate(i+1,k),alt(j),ProdRate(j,k))
		endif
      enddo
   enddo
enddo

!make sure prod rate doesn't get too low
do i=1,NGrid
   do j=1,NIon
	if (ProdRate(i,j) .LT. 1E-20) then
		ProdRate(i,j)=1E-20
	else
		continue
	endif
   enddo
enddo

!!!read in electemp file, interpolate electron temperatures on our s grid
!units are: Altitude: [km]; Elec temp: [K]
do i=1,m
	read(28,*) orig_alt(i),orig_te(i)
enddo

do i=1,m
	do j=1,NGrid
		if (alt(j)==orig_alt(i)) then
			t_e(j)=orig_te(i)
		elseif(alt(j)>orig_alt(i) .and. alt(j)<orig_alt(i+1)) then
			call linearint(orig_alt(i), orig_alt(i+1), orig_te(i), orig_te(i+1), alt(j), t_e(j))
		endif
	enddo
enddo

!!!read in iontemp file, interpolate electron temperatures on our s grid
!units are: Altitude: [km]; Ion temp: [K]
read(29,*)
do i=1,mm
	read(29,*) orig_alt2(i),orig_ti(i)
enddo

do i=1,mm
	do j=1,NGrid
		if (alt(j)==orig_alt2(i)) then
			t_i(j)=orig_ti(i)
		elseif(alt(j)>orig_alt2(i) .and. alt(j)<orig_alt2(i+1)) then
			call linearint(orig_alt2(i), orig_alt2(i+1), orig_ti(i), orig_ti(i+1), alt(j), t_i(j))
		endif
	enddo
enddo

!!! run subroutines for R(NGrid), z(NGrid), A(NGrid), gpar(NGrid), atlde(NGrid)
DO i=1,NGrid
	R(i) = Rmars + alt(i)                                     ![km]
	call area(A0, x(i), y(i), A(i))                           ![units of A0]
	call gparallel(Rmars, alt(i), x(i), GPAR(i))              ![N/kg]=[m/s^2]
enddo

call atilde(A(1),A(2), A(1), s(2), s(1), Atld(1)) ![1/m]
do i=2,NGrid1
	call atilde(A(i),A(i+1), A(i-1), s(i+1), s(i-1), Atld(i)) ![1/m]
enddo
call atilde(A(125),A(125), A(124), s(125), s(124), Atld(125)) ![1/m]

!!! Define CSS0 and alpha terms for each ion species
!!! Ion order: 1 O2+; 2 O+; 3 CO2+; 4 CO+; 5 NO+; 6 N2+; 7 OH+; 8 H+; 9 HCO+
! CSS0 = k*T/m [(m/s)^2]
!k=1.38E-23, T=temp(i) = t_e + t_i
DO i=1,NGrid
	DO j=1,NIon
	   CSS0(i,j)= 1.38E-23*(t_i(i)+t_e(i))/mass(j)		![(m/s)^2]
	ENDDO
	alpha(i,1)=2.4E-7*(300/t_e(i))**(0.7)			![cm^3/s], O2+
	alpha(i,2)=4.2E-12*(300/t_e(i))**(0.7)			![cm^3/s], O+
	alpha(i,3)=3.5E-7*(300/t_e(i))**(0.5)			![cm^3/s], CO2+
	alpha(i,4)=2.75E-7*(300/t_e(i))**(0.5)			![cm^3/s], CO+
	alpha(i,5)=4.0E-7*(300/t_e(i))**(0.5)			![cm^3/s], NO+
	alpha(i,6)=2.2E-7*(300/t_e(i))**(0.39)			![cm^3/s], N2+
	alpha(i,7)=3.75E-8*(300/t_e(i))**(0.5)			![cm^3/s], OH+
	alpha(i,8)=5.45E-12*(300/t_e(i))**(0.7)			![cm^3/s], H+
	alpha(i,9)=1.1E-7*(300/t_e(i))**(1.0)			![cm^3/s], HCO+
ENDDO

!DIAGNOSTIC: test out that these subroutines are all being read in correctly
write(35,*)"Checking Subroutines and Inputs"
write(35,*) "i s x altitude radius area a-tilde gpar Te Ti"
write(36,*)"Checking Neutral Density Input [1/cm^3] by species"
write(36,*) "i s altitude N2 CO2 CO O O2 N H NO H2 H2O"
write(37,*)"Checking Photoionization Input [cm^-3/s] by species"
write(37,*) "i s altitude O2+ O+ CO2+ CO+ NO+ N2+ OH+ H+ HCO+"
write(38,*)"Checking alpha values [cm^3/s] by species"
write(38,*)"i s altitude O2+ O+ CO2+ CO+ NO+ N2+ OH+ H+ HCO+"
do i=1,NGrid
	write(35,*) i,s(i),x(i),alt(i),R(i),A(i),Atld(i),gpar(i),t_e(i),t_i(i)
	write(36,*) i,s(i),alt(i),(NeutDen(i,j),j=1,NNeut)
	write(37,*) i,s(i),alt(i),(ProdRate(i,j),j=1,NIon)
	write(38,*) i,s(i),alt(i),(alpha(i,j),j=1,NIon)
enddo

!! set initial conditions, try reading in from init_condtions.txt file (if continuing previous test)
!!***adjust initial conditions (looking at photochemical). ni(NGrid) =/= ne(NGrid)
Do i=2,NGrid1

   w1(i,1)=(A(i)*ProdRate(i,1)*1E6)/(alpha(i,1)*1E-6*NElec(i)/A(i)+&
		4.5E-10*NeutDen(i,8)+1.00E-10*NeutDen(i,6))
   w1(i,2)=(A(i)*ProdRate(i,2)*1E6)/(alpha(i,2)*1E-6*NElec(i)/A(i)+&
		1.10E-9*NeutDen(i,2)+5.70E-10*NeutDen(i,1))
   w1(i,3)=(A(i)*ProdRate(i,3)*1E6)/(alpha(i,3)*1E-6*NElec(i)/A(i)+&
		1.64E-10*NeutDen(i,4)+9.60E-11*NeutDen(i,4))
   w1(i,4)=(A(i)*ProdRate(i,4)*1E6)/(alpha(i,4)*1E-6*NElec(i)/A(i)+&
		1.10E-09*NeutDen(i,2)+1.40E-10*NeutDen(i,4))
   w1(i,5)=(A(i)*ProdRate(i,5)*1E6)/(alpha(i,5)*1E-6*NElec(i)/A(i)+&
		8.00E-10*NeutDen(i,2)+1.30E-10*NeutDen(i,4))
   w1(i,6)=(A(i)*ProdRate(i,6)*1E6)/(alpha(i,6)*1E-6*NElec(i)/A(i)+&
		8.00E-10*NeutDen(i,2)+1.30E-10*NeutDen(i,4))
   w1(i,7)=A(i)*1E5
   w1(i,8)=(A(i)*ProdRate(i,8)*1E6)/(alpha(i,8)*1E-6*NElec(i)/A(i)+&
		3.80E-09*NeutDen(i,2))
   w1(i,9)=A(i)*1E6

!   w1(i,1)=A(i)*1E5
!   w1(i,2)=A(i)*1E5
!   w1(i,3)=A(i)*1E5
!   w1(i,4)=A(i)*1E5
!   w1(i,5)=A(i)*1E5
!   w1(i,6)=A(i)*1E5
!   w1(i,7)=A(i)*1E5
!   w1(i,8)=A(i)*1E5
!   w1(i,9)=A(i)*1E5

   NElec(i)=w1(i,1)+w1(i,2)+w1(i,3)+w1(i,4)+w1(i,5)+w1(i,6)+w1(i,7)+w1(i,8)+w1(i,9)
   Do j=1,NIon
	w2(i,j)=0.0
	F1(i,j)=w2(i,j)
	F2(i,j)=w2(i,j)*w2(i,j)/w1(i,j)+CSS0(i,j)*w1(i,j)
   enddo
enddo

!!!set boundary conditions
!set ion densities to photochemical at boundaries
w1(1,1)=A(1)*6.96E10			!O2+ density [m2]*[1/m3]=[1/m]
w1(1,2)=A(1)*9.028E7			!O+ density [m2]*[1/m3]=[1/m]
w1(1,3)=A(1)*9.221E9			!CO2+ density [m2]*[1/m3]=[1/m]
w1(1,4)=A(1)*5.000E7			!CO+ density [m2]*[1/m3]=[1/m]
w1(1,5)=A(1)*5.195E9			!NO+ density [m2]*[1/m3]=[1/m]
w1(1,6)=A(1)*7.898E6			!N2+ density [m2]*[1/m3]=[1/m]
w1(1,7)=A(1)*5.000E3			!OH+ density [m2]*[1/m3]=[1/m]
w1(1,8)=A(1)*1.000E4			!H+ density [m2]*[1/m3]=[1/m]
w1(1,9)=A(1)*2.392E8			!HCO+ density [m2]*[1/m3]=[1/m]

do j=1,NIon
	!set si=NGrid endpoint to be same as si=1 endpoint
	w1(NGrid,j)=w1(1,j)
	!set velocities at boundary to 0 [m/s]
	w2(1,j)=0.0
	w2(NGrid,j)=0.0
enddo

NElec(1)=w1(1,1)+w1(1,2)+w1(1,3)+w1(1,4)+w1(1,5)+w1(1,6)+w1(1,7)+w1(1,8)+w1(1,9)
NElec(NGrid)=NElec(1)

!do i=1,NGrid
!	write(*,*)i,s(i),NElec(i)
!enddo

!if we have other conditions to start from, uncomment this section to read in files 31,32 
do i=1,NGrid
	read(31,*,end=199)toss(i),w1(i,1),w1(i,2),w1(i,3),w1(i,4),w1(i,5),w1(i,6),&
		w1(i,7),w1(i,8),w1(i,9)
	read(32,*,end=199)toss(i),w2(i,1),w2(i,2),w2(i,3),w2(i,4),w2(i,5),w2(i,6),&
		w2(i,7),w2(i,8),w2(i,9)
	!make sure err=# is same # as the line# of enddo statement
	NElec(i)=w1(i,1)+w1(i,2)+w1(i,3)+w1(i,4)+w1(i,5)+w1(i,6)+w1(i,7)+w1(i,8)+w1(i,9)
	do j=1,NIon
		F1(i,j)=w2(i,j)
		F2(i,j)=w2(i,j)*w2(i,j)/w1(i,j)+CSS0(i,j)*w1(i,j)
		!Fhp1(i)=w2p(i)
		!Fhp2(i)=w2p(i)*w2p(i)/w1p(i)+CSS0(i)*w1p(i)
		!Fhm1(i)=w2m(i)
		!Fhm2(i)=w2m(i)*w2m(i)/w1m(i)+CSS0(i)*w1m(i)
	enddo
enddo
199 continue

!!!DIAGNOSTIC: check that initial values are correct
write(121,*)"Check Initial Conditions - w1 values (Area Weighted Densities)"
write(121,*)"i  s  O2+  O+  CO2+  CO+  NO+  N2+  OH+  H+  HCO+"
write(122,*)"Check Initial Conditions - w2 values (Area Weighted Velocities)"
write(122,*)"i  s  O2+  O+  CO2+  CO+  NO+  N2+  OH+  H+  HCO+"
do i=1,NGrid
	write(121,*) i, s(i), (w1(i,j),j=1,NIon)
	write(122,*) i, s(i), (w2(i,j),j=1,NIon)
enddo

!save initial conditions *as time point 0*
DO i=1,NGrid
	do j=1,NIon
		den_i(i,j)=w1(i,j)/A(i)
		uu_i(i,j)=w2(i,j)/A(i)/den_i(i,j)
	enddo
	elecden_i(i)=NElec(i)/A(i)
enddo

!!!DIAGNOSTIC: setup diagnostic files
write(102,*)"Checking Right Hand Terms at s=",s(30)/1000.0,"and alt=",alt(30)
write(102,*)"Time  Production  Loss  Gravity  Collision  ElecPressure  MassLoading  f1*A  f2*A"
write(103,*)"Checking Right Hand Terms at s=",s(30)/1000.0,"and alt=",alt(30)
write(103,*)"Time  Production  Loss  Gravity  Collision  ElecPressure  MassLoading  f1*A  f2*A"
write(104,*)"Checking Right Hand Terms at s=",s(30)/1000.0,"and alt=",alt(30)
write(104,*)"Time  Production  Loss  Gravity  Collision  ElecPressure  MassLoading  f1*A  f2*A"
write(105,*)"Checking Right Hand Terms at s=",s(30)/1000.0,"and alt=",alt(30)
write(105,*)"Time  Production  Loss  Gravity  Collision  ElecPressure  MassLoading  f1*A  f2*A"
write(106,*)"Checking Right Hand Terms at s=",s(30)/1000.0,"and alt=",alt(30)
write(106,*)"Time  Production  Loss  Gravity  Collision  ElecPressure  MassLoading  f1*A  f2*A"
write(107,*)"Checking Right Hand Terms at s=",s(30)/1000.0,"and alt=",alt(30)
write(107,*)"Time  Production  Loss  Gravity  Collision  ElecPressure  MassLoading  f1*A  f2*A"
write(108,*)"Checking Right Hand Terms at s=",s(30)/1000.0,"and alt=",alt(30)
write(108,*)"Time  Production  Loss  Gravity  Collision  ElecPressure  MassLoading  f1*A  f2*A"
write(109,*)"Checking Right Hand Terms at s=",s(30)/1000.0,"and alt=",alt(30)
write(109,*)"Time  Production  Loss  Gravity  Collision  ElecPressure  MassLoading  f1*A  f2*A"
write(110,*)"Checking Right Hand Terms at s=",s(30)/1000.0,"and alt=",alt(30)
write(110,*)"Time  Production  Loss  Gravity  Collision  ElecPressure  MassLoading  f1*A  f2*A"

!!!!!!!!!!!!!!!!!!!!!
!!!! time loop over K  advance from K to K+1 via 4 steps
!!!!  STEP 1 is a half step from K to K+1/2 to find w1 and w2 at the half time-step
!!!!  K+1/2 from the K values.
!!!!!!!!!!!!!!!!!!!!!

!!Create time array
read (30,*,end=145) newtime  !!Need to read in newtime from init_conditions file, if we have one
145 continue
close(30)

write(*,*) 'Starting time is:', newtime

!!check s grid is correct
!do i=1,NGrid
!	write(*,*) i,s(i)
!enddo

!!!Loop through time
DO K=1,NSteps
   time(k)=k*DELT+newtime 			!Sets the total elapsed time of code from initial conditions
   write(*,*)'time step is:',K,'time is:',time(k)

!!!Loop through ions
   DO J=1,NIon
!	write(*,*) 'Looping through Ion:',j

!Step 0: Update F1,F2,R1,R2 values
      DO i=1,NGrid
	F1(i,j)=w2(i,j)
	F2(i,j)=w2(i,j)*w2(i,j)/w1(i,j)+CSS0(i,j)*w1(i,j)

	!!!!!define prod and loss terms by species
	!prod(i,j)=A(i)*ProdRate(i,j)*1E6
	!loss(i,j)=alpha(i,j)*1E-6*w1(i,j)*NElec(i)/A(i)
	!O2+
	prod(i,1)=A(i)*ProdRate(i,1)*1E6+&				!Photoionization
		1.64E-10*w1(i,3)*NeutDen(i,4)+&				!CO2+ + O -> O2+        
		1.10E-09*w1(i,2)*NeutDen(i,2)+&				!O+ + CO2 -> O2+		
		7.10E-10*w1(i,7)*NeutDen(i,4)+&				!OH+ + O -> O2+		
		1.17E-09*w1(i,8)*NeutDen(i,5)+&				!H+ + O2 -> O2+
		1.50E-10*(300/t_i(i))**(1.1)*w1(i,4)*NeutDen(i,5)+&	!CO+ + O2 -> O2+	
		5.10E-11*(300/t_i(i))**(1.16)*w1(i,6)*NeutDen(i,5)+&	!N2+ + O2 -> O2+	
		1.60E-10*(300/t_i(i))**(0.52)*w1(i,2)*NeutDen(i,5)	!O+ + O2 -> O2+		
	loss(i,1)=alpha(i,1)*1E-6*w1(i,1)*NElec(i)/A(i)+&		!Dissociative recombination
		4.60E-10*w1(i,1)*NeutDen(i,8)+&				!O2+ + NO -> NO+
		1.00E-10*w1(i,1)*NeutDen(i,6)+&				!O2+ + N -> NO+
		1.00E-15*w1(i,1)*NeutDen(i,1)				!O2+ + N2 -> NO+
	!O+
	prod(i,2)=A(i)*ProdRate(i,2)*1E6+&				!Photoionization
		9.60E-11*w1(i,3)*NeutDen(i,4)+&				!CO2+ + O -> O+
		3.75E-10*w1(i,8)*NeutDen(i,4)+&				!H+ + O -> O+
		7.00E-12*(300/t_i(i))**(0.23)*w1(i,6)*NeutDen(i,4)+&	!N2+ + O -> O+
		1.40E-10*w1(i,4)*NeutDen(i,4)				!CO+ + O -> O+
	loss(i,2)=alpha(i,2)*1E-6*w1(i,2)*NElec(i)/A(i)+&		!Radiative recombination
		1.10E-09*w1(i,2)*NeutDen(i,2)+&				!O+ + CO2 -> O2+
		1.60E-11*(300/t_i(i))**(0.52)*w1(i,2)*NeutDen(i,5)+&	!O+ + O2 -> O2+
		8.00E-13*w1(i,2)*NeutDen(i,8)+&				!O+ + NO -> NO+
		1.20E-12*(300/t_i(i))**(0.45)*w1(i,2)*NeutDen(i,1)+&	!O+ + N2 -> NO+
		1.67E-09*w1(i,2)*NeutDen(i,9)+&				!O+ + H2 -> OH+
		6.40E-10*w1(i,2)*NeutDen(i,7)				!O+ + H -> H+
	!CO2+
	prod(i,3)=A(i)*ProdRate(i,3)*1E6+&				!Photoionization
		9.00E-10*(300/t_i(i))**(0.23)*w1(i,6)*NeutDen(i,2)+&	!N2+ + CO2 -> CO2+
		1.10E-09*w1(i,4)*NeutDen(i,2)				!CO+ + CO2 -> CO2+
	loss(i,3)=alpha(i,3)*1E-6*w1(i,3)*NElec(i)/A(i)+&		!Dissociative recombination
		1.64E-10*w1(i,3)*NeutDen(i,4)+&				!CO2+ + O -> O2+
		9.60E-11*w1(i,3)*NeutDen(i,4)+&				!CO2+ + O -> O+
		3.40E-10*w1(i,3)*NeutDen(i,6)+&				!CO2+ + N -> CO+
		1.23E-10*w1(i,3)*NeutDen(i,8)+&				!CO2+ + NO -> NO+
		2.35E-11*w1(i,3)*NeutDen(i,7)+&				!CO2+ + H -> H+
		4.70E-10*w1(i,3)*NeutDen(i,7)				!CO2+ + H -> HCO+
	!CO+
	prod(i,4)=A(i)*ProdRate(i,4)*1E6+&				!Photoionization
		3.40E-10*w1(i,3)*NeutDen(i,6)+&				!CO2+ + N -> CO+
		3.55E-10*w1(i,7)*NeutDen(i,3)+&				!OH+ + CO -> CO+
		7.60E-11*w1(i,6)*NeutDen(i,3)				!N2+ + CO -> CO+
	loss(i,4)=alpha(i,4)*1E-6*w1(i,4)*NElec(i)/A(i)+&		!Dissociative recombination
		1.10E-09*w1(i,4)*NeutDen(i,2)+&				!CO+ + CO2 -> CO2+
		1.50E-10*(300/t_i(i))**(1.1)*w1(i,4)*NeutDen(i,5)+&	!CO+ + O2 -> O2+
		1.40E-10*w1(i,4)*NeutDen(i,4)+&				!CO+ + O -> O+
		4.20E-10*w1(i,4)*NeutDen(i,8)+&				!CO+ + NO -> NO+
		8.20E-11*w1(i,4)*NeutDen(i,6)+&				!CO+ + N -> NO+
		4.00E-10*w1(i,4)*NeutDen(i,7)+&				!CO+ + H -> H+
		7.50E-10*w1(i,4)*NeutDen(i,9)				!CO+ + H2 -> HCO+
	!NO+
	prod(i,5)=A(i)*ProdRate(i,5)*1E6+&				!Photoionization
		1.33E-10*(300/t_i(i))**(0.44)*w1(i,6)*NeutDen(i,4)+&	!N2+ + O -> NO+
		1.23E-10*w1(i,3)*NeutDen(i,8)+&				!CO2+ + NO -> NO+
		4.20E-10*w1(i,4)*NeutDen(i,8)+&				!CO+ + NO -> NO+
		4.60E-10*w1(i,1)*NeutDen(i,8)+&				!O2+ + NO -> NO+
		7.00E-13*(t_i(i)/300)**(0.87)*w1(i,2)*NeutDen(i,8)+&	!O+ + NO -> NO+
		1.90E-09*w1(i,8)*NeutDen(i,8)+&				!H+ + NO -> NO+
		8.20E-11*w1(i,4)*NeutDen(i,6)+&				!CO+ + N -> NO+
		1.00E-10*w1(i,1)*NeutDen(i,6)+&				!O2+ + N -> NO+
		1.20E-12*(300/t_i(i))**(0.45)*w1(i,2)*NeutDen(i,1)+&	!O+ + N2 -> NO+
		1.00E-15*w1(i,1)*NeutDen(i,1)				!O2+ + N2 -> NO+
	loss(i,5)=alpha(i,5)*1E-6*w1(i,5)*NElec(i)/A(i)			!Dissociative recombination
	!N2+
	prod(i,6)=A(i)*ProdRate(i,6)*1E6				!Photoionization
	loss(i,6)=alpha(i,6)*1E-6*w1(i,6)*NElec(i)/A(i)+&		!Dissociative recombination
		5.10E-11*(300/t_i(i))**(1.16)*w1(i,6)*NeutDen(i,5)+&	!N2+ + O2 -> O2+
		7.00E-12*(300/t_i(i))**(0.23)*w1(i,6)*NeutDen(i,4)+&	!N2+ + O -> O+
		9.00E-10*(300/t_i(i))**(0.23)*w1(i,6)*NeutDen(i,2)+&	!N2+ + CO2 -> CO2+
		7.60E-11*w1(i,6)*NeutDen(i,3)+&				!N2+ + CO -> CO+
		1.33E-10*(300/t_i(i))**(0.44)*w1(i,6)*NeutDen(i,4)	!N2+ + O -> NO+
	!OH+
	prod(i,7)=A(i)*ProdRate(i,7)*1E6+&				!Photoionization
		1.67E-09*w1(i,2)*NeutDen(i,9)				!O+ + H2 -> OH+
	loss(i,7)=alpha(i,7)*1E-6*w1(i,7)*NElec(i)/A(i)+&		!Dissociative recombination
		7.10E-09*w1(i,7)*NeutDen(i,4)+&				!OH+ + O -> O2+
		3.55E-10*w1(i,7)*NeutDen(i,3)+&				!OH+ + CO -> CO+
		3.55E-10*w1(i,7)*NeutDen(i,3)+&				!OH+ + CO -> HCO+
		1.10E-09*w1(i,7)*NeutDen(i,2)+&				!OH+ + CO2 -> HCO2+
		2.40E-10*w1(i,7)*NeutDen(i,1)				!OH+ + N2 -> N2H+
	!H+
	prod(i,8)=A(i)*ProdRate(i,8)*1E6+&				!Photoionization
		2.35E-11*w1(i,3)*NeutDen(i,7)+&				!CO2+ + H -> H+
		4.00E-10*w1(i,4)*NeutDen(i,7)+&				!CO+ + H -> H+
		6.40E-10*w1(i,2)*NeutDen(i,7)				!O+ + H -> H+
	loss(i,8)=alpha(i,8)*1E-6*w1(i,8)*NElec(i)/A(i)+&		!Radiative recombination
		3.75E-10*w1(i,8)*NeutDen(i,4)+&				!H+ + O -> O+
		1.90E-09*w1(i,8)*NeutDen(i,8)+&				!H+ + NO -> NO+
		3.80E-09*w1(i,8)*NeutDen(i,2)				!H+ + CO2 -> HCO+
	!HCO+
	prod(i,9)=A(i)*ProdRate(i,9)*1E6+&				!Photoionization
		7.50E-10*w1(i,4)*NeutDen(i,9)+&				!CO+ + H2 -> HCO+
		3.80E-09*w1(i,8)*NeutDen(i,2)+&				!H+ + CO2 -> HCO+
		4.70E-10*w1(i,3)*NeutDen(i,7)+&				!CO2+ + H -> HCO+
		3.55E-10*w1(i,7)*NeutDen(i,3)				!OH+ + CO -> HCO+
	loss(i,9)=alpha(i,9)*1E-6*w1(i,9)*NElec(i)/A(i)			!Dissociative recombination

	!done with production/loss terms

	!!!!!!!define collision terms by species
	!O2+
	coll(i,1)=w1(i,1)*(0.13*w1(i,2)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,1)/w1(i,1)-w2(i,2)/w1(i,2))+&
		0.17*w1(i,3)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,1)/w1(i,1)-w2(i,3)/w1(i,3))+&
		0.15*w1(i,4)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,1)/w1(i,1)-w2(i,4)/w1(i,4))+&
		0.16*w1(i,5)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,1)/w1(i,1)-w2(i,5)/w1(i,5))+&
		0.15*w1(i,6)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,1)/w1(i,1)-w2(i,6)/w1(i,6))+&
		0.13*w1(i,7)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,1)/w1(i,1)-w2(i,7)/w1(i,7))+&
		0.039*w1(i,8)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,1)/w1(i,1)-w2(i,8)/w1(i,8))+&
		0.15*w1(i,9)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,1)/w1(i,1)-w2(i,9)/w1(i,9))+&
		(4.13*NeutDen(i,1)+5.63*NeutDen(i,2)+4.37*NeutDen(i,3)+2.31*NeutDen(i,4)+&
		1.89*NeutDen(i,5)+2.64*NeutDen(i,6)+0.65*NeutDen(i,7))*1E-10*w2(i,1)/w1(i,1))
	!O+
	coll(i,2)=w1(i,2)*(0.26*w1(i,1)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,2)/w1(i,2)-w2(i,1)/w1(i,1))+&
		0.27*w1(i,3)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,2)/w1(i,2)-w2(i,3)/w1(i,3))+&
		0.25*w1(i,4)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,2)/w1(i,2)-w2(i,4)/w1(i,4))+&
		0.26*w1(i,5)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,2)/w1(i,2)-w2(i,5)/w1(i,5))+&
		0.25*w1(i,6)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,2)/w1(i,2)-w2(i,6)/w1(i,6))+&
		0.23*w1(i,7)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,2)/w1(i,2)-w2(i,7)/w1(i,7))+&
		0.077*w1(i,8)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,2)/w1(i,2)-w2(i,8)/w1(i,8))+&
		0.25*w1(i,9)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,2)/w1(i,2)-w2(i,9)/w1(i,9))+&
		(6.82*NeutDen(i,1)+8.95*NeutDen(i,2)+7.22*NeutDen(i,3)+3.34*NeutDen(i,4)+&
		6.64*NeutDen(i,5)+4.62*NeutDen(i,6)+8.59*NeutDen(i,7))*1E-10*w2(i,2)/w1(i,2))
	!CO2+
	coll(i,3)=w1(i,3)*(0.12*w1(i,1)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,3)/w1(i,3)-w2(i,1)/w1(i,1))+&
		0.10*w1(i,2)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,3)/w1(i,3)-w2(i,2)/w1(i,2))+&
		0.12*w1(i,4)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,3)/w1(i,3)-w2(i,4)/w1(i,4))+&
		0.12*w1(i,5)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,3)/w1(i,3)-w2(i,5)/w1(i,5))+&
		0.12*w1(i,6)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,3)/w1(i,3)-w2(i,6)/w1(i,6))+&
		0.10*w1(i,7)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,3)/w1(i,3)-w2(i,7)/w1(i,7))+&
		0.029*w1(i,8)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,3)/w1(i,3)-w2(i,8)/w1(i,8))+&
		0.12*w1(i,9)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,3)/w1(i,3)-w2(i,9)/w1(i,9))+&
		(3.22*NeutDen(i,1)+1.58*NeutDen(i,2)+3.40*NeutDen(i,3)+1.76*NeutDen(i,4)+&
		3.18*NeutDen(i,5)+2.00*NeutDen(i,6)+0.47*NeutDen(i,7))*1E-10*w2(i,3)/w1(i,3))
	!CO+
	coll(i,4)=w1(i,4)*(0.18*w1(i,1)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,4)/w1(i,4)-w2(i,1)/w1(i,1))+&
		0.15*w1(i,2)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,4)/w1(i,4)-w2(i,2)/w1(i,2))+&
		0.19*w1(i,3)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,4)/w1(i,4)-w2(i,3)/w1(i,3))+&
		0.17*w1(i,5)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,4)/w1(i,4)-w2(i,5)/w1(i,5))+&
		0.17*w1(i,6)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,4)/w1(i,4)-w2(i,6)/w1(i,6))+&
		0.15*w1(i,7)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,4)/w1(i,4)-w2(i,7)/w1(i,7))+&
		0.045*w1(i,8)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,4)/w1(i,4)-w2(i,8)/w1(i,8))+&
		0.17*w1(i,9)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,4)/w1(i,4)-w2(i,9)/w1(i,9))+&
		(4.24*NeutDen(i,1)+6.18*NeutDen(i,2)+1.79*NeutDen(i,3)+2.58*NeutDen(i,4)+&
		4.49*NeutDen(i,5)+2.95*NeutDen(i,6)+0.74*NeutDen(i,7))*1E-10*w2(i,4)/w1(i,4))
	!NO+
	coll(i,5)=w1(i,5)*(0.17*w1(i,1)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,5)/w1(i,5)-w2(i,1)/w1(i,1))+&
		0.14*w1(i,2)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,5)/w1(i,5)-w2(i,2)/w1(i,2))+&
		0.18*w1(i,3)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,5)/w1(i,5)-w2(i,3)/w1(i,3))+&
		0.16*w1(i,4)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,5)/w1(i,5)-w2(i,4)/w1(i,4))+&
		0.16*w1(i,6)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,5)/w1(i,5)-w2(i,6)/w1(i,6))+&
		0.14*w1(i,7)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,5)/w1(i,5)-w2(i,7)/w1(i,7))+&
		0.042*w1(i,8)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,5)/w1(i,5)-w2(i,8)/w1(i,8))+&
		0.16*w1(i,9)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,5)/w1(i,5)-w2(i,9)/w1(i,9))+&
		(4.34*NeutDen(i,1)+5.89*NeutDen(i,2)+4.59*NeutDen(i,3)+2.44*NeutDen(i,4)+&
		4.27*NeutDen(i,5)+2.79*NeutDen(i,6)+0.69*NeutDen(i,7))*1E-10*w2(i,5)/w1(i,5))
	!N2+
	coll(i,6)=w1(i,6)*(0.18*w1(i,1)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,6)/w1(i,6)-w2(i,1)/w1(i,1))+&
		0.15*w1(i,2)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,6)/w1(i,6)-w2(i,2)/w1(i,2))+&
		0.19*w1(i,3)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,6)/w1(i,6)-w2(i,3)/w1(i,3))+&
		0.17*w1(i,4)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,6)/w1(i,6)-w2(i,4)/w1(i,4))+&
		0.17*w1(i,5)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,6)/w1(i,6)-w2(i,5)/w1(i,5))+&
		0.15*w1(i,7)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,6)/w1(i,6)-w2(i,7)/w1(i,7))+&
		0.045*w1(i,8)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,6)/w1(i,6)-w2(i,8)/w1(i,8))+&
		0.17*w1(i,9)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,6)/w1(i,6)-w2(i,9)/w1(i,9))+&
		(4.15*NeutDen(i,1)+6.18*NeutDen(i,2)+4.84*NeutDen(i,3)+2.58*NeutDen(i,4)+&
		4.49*NeutDen(i,5)+2.95*NeutDen(i,6)+0.74*NeutDen(i,7))*1E-10*w2(i,6)/w1(i,6))
	!OH+
	coll(i,7)=w1(i,7)*(0.25*w1(i,1)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,7)/w1(i,7)-w2(i,1)/w1(i,1))+&
		0.21*w1(i,2)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,7)/w1(i,7)-w2(i,2)/w1(i,2))+&
		0.26*w1(i,3)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,7)/w1(i,7)-w2(i,3)/w1(i,3))+&
		0.24*w1(i,4)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,7)/w1(i,7)-w2(i,4)/w1(i,4))+&
		0.25*w1(i,5)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,7)/w1(i,7)-w2(i,5)/w1(i,5))+&
		0.24*w1(i,6)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,7)/w1(i,7)-w2(i,6)/w1(i,6))+&
		0.073*w1(i,8)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,7)/w1(i,7)-w2(i,8)/w1(i,8))+&
		0.24*w1(i,9)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,7)/w1(i,7)-w2(i,9)/w1(i,9))+&
		(10.0*NeutDen(i,1)+12.0*NeutDen(i,2)+10.0*NeutDen(i,3)+5.00*NeutDen(i,4)+&
		10.0*NeutDen(i,5)+5.00*NeutDen(i,6)+1.00*NeutDen(i,7))*1E-10*w2(i,7)/w1(i,7))
	!H+
	coll(i,8)=w1(i,8)*(1.25*w1(i,1)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,8)/w1(i,8)-w2(i,1)/w1(i,1))+&
		1.23*w1(i,2)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,8)/w1(i,8)-w2(i,2)/w1(i,2))+&
		1.26*w1(i,3)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,8)/w1(i,8)-w2(i,3)/w1(i,3))+&
		1.25*w1(i,4)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,8)/w1(i,8)-w2(i,4)/w1(i,4))+&
		1.25*w1(i,5)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,8)/w1(i,8)-w2(i,5)/w1(i,5))+&
		1.25*w1(i,6)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,8)/w1(i,8)-w2(i,6)/w1(i,6))+&
		1.23*w1(i,7)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,8)/w1(i,8)-w2(i,7)/w1(i,7))+&
		1.25*w1(i,9)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,8)/w1(i,8)-w2(i,9)/w1(i,9))+&
		(33.6*NeutDen(i,1)+41.4*NeutDen(i,2)+35.6*NeutDen(i,3)+8.59*NeutDen(i,4)+&
		32.0*NeutDen(i,5)+26.1*NeutDen(i,6)+14.7*NeutDen(i,7))*1E-10*w2(i,8)/w1(i,8))
	!HCO+
	coll(i,9)=w1(i,9)*(0.17*w1(i,1)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,9)/w1(i,9)-w2(i,1)/w1(i,1))+&
		0.14*w1(i,2)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,9)/w1(i,9)-w2(i,2)/w1(i,2))+&
		0.18*w1(i,3)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,9)/w1(i,9)-w2(i,3)/w1(i,3))+&
		0.17*w1(i,4)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,9)/w1(i,9)-w2(i,4)/w1(i,4))+&
		0.17*w1(i,5)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,9)/w1(i,9)-w2(i,5)/w1(i,5))+&
		0.17*w1(i,6)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,9)/w1(i,9)-w2(i,6)/w1(i,6))+&
		0.14*w1(i,7)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,9)/w1(i,9)-w2(i,7)/w1(i,7))+&
		0.043*w1(i,8)/A(i)*1E-6/(t_i(i))**(3.0/2.0)*(w2(i,9)/w1(i,9)-w2(i,8)/w1(i,8))+&
		(4.30*NeutDen(i,1)+6.00*NeutDen(i,2)+4.70*NeutDen(i,3)+2.50*NeutDen(i,4)+&
		4.38*NeutDen(i,5)+2.87*NeutDen(i,6)+0.72*NeutDen(i,7))*1E-10*w2(i,9)/w1(i,9))

	grav(i,j)=w1(i,j)*gpar(i)

	if (i==1) then
		epres(i,j)=1.38E-23/mass(j)*((w1(i+1,j)*t_e(i+1)-w1(i,j)*t_e(i))/(DELS)-&
			w1(i,j)/NElec(i)*(NElec(i+1)*t_e(i+1)-NElec(i)*t_e(i))/(DELS))
	elseif (i==NGrid) then
		epres(i,j)=1.38E-23/mass(j)*((w1(i,j)*t_e(i)-w1(i-1,j)*t_e(i-1))/(DELS)-&
			w1(i,j)/NElec(i)*(NElec(i)*t_e(i)-NElec(i-1)*t_e(i-1))/(DELS))
	else
		epres(i,j)=1.38E-23/mass(j)*((w1(i+1,j)*t_e(i+1)-w1(i-1,j)*t_e(i-1))/(2*DELS)-&
			w1(i,j)/NElec(i)*(NElec(i+1)*t_e(i+1)-NElec(i-1)*t_e(i-1))/(2*DELS))
	endif

	mload(i,j)=-w2(i,j)/w1(i,j)*loss(i,j)

	R1(i,j)=prod(i,j)-loss(i,j)+f1(i,j)*Atld(i)
	R2(i,j)=grav(i,j)-coll(i,j)+epres(i,j)+mload(i,j)+f2(i,j)*Atld(i)
      enddo

!!DIAGNOSTICS: Check RHS values at 3 different altitudes
!write(*,*) "Low Altitude - RHS1 for ion",j,"at altitude",alt(12),"is:",R1(12,j)
!write(*,*) "Low Altitude - RHS2 for ion",j,"at location",alt(12),"is:",R2(12,j)
!write(*,*) "Mid Altitude - RHS1 for ion",j,"at altitude",alt(33),"is:",R1(33,j)
!write(*,*) "Mid Altitude - RHS2 for ion",j,"at location",alt(33),"is:",R2(33,j)
!write(*,*) "High Altitude - RHS1 for ion",j,"at altitude",alt(63),"is:",R1(63,j)
!write(*,*) "High Altitude - RHS2 for ion",j,"at location",alt(63),"is:",R2(63,j)
!!DIAGNOSTICS: Check Right Hand side terms
if (j==1) then
write(102,*)time(k),prod(30,1),loss(30,1),grav(30,1),coll(30,1),epres(30,1),mload(30,1),(f1(30,1)*Atld(30)),(f2(30,1)*Atld(30))
elseif (j==2) then
write(103,*)time(k),prod(30,2),loss(30,2),grav(30,2),coll(30,2),epres(30,2),mload(30,2),(f1(30,2)*Atld(30)),(f2(30,2)*Atld(30))
elseif (j==3) then
write(104,*)time(k),prod(30,3),loss(30,3),grav(30,3),coll(30,3),epres(30,3),mload(30,3),(f1(30,3)*Atld(30)),(f2(30,3)*Atld(30))
elseif (j==4) then
write(105,*)time(k),prod(30,4),loss(30,4),grav(30,4),coll(30,4),epres(30,4),mload(30,4),(f1(30,4)*Atld(30)),(f2(30,4)*Atld(30))
elseif (j==5) then
write(106,*)time(k),prod(30,5),loss(30,5),grav(30,5),coll(30,5),epres(30,5),mload(30,5),(f1(30,5)*Atld(30)),(f2(30,5)*Atld(30))
elseif (j==6) then
write(107,*)time(k),prod(30,6),loss(30,6),grav(30,6),coll(30,6),epres(30,6),mload(30,6),(f1(30,6)*Atld(30)),(f2(30,6)*Atld(30))
elseif (j==7) then
write(108,*)time(k),prod(30,7),loss(30,7),grav(30,7),coll(30,7),epres(30,7),mload(30,7),(f1(30,7)*Atld(30)),(f2(30,7)*Atld(30))
elseif (j==8) then
write(109,*)time(k),prod(30,8),loss(30,8),grav(30,8),coll(30,8),epres(30,8),mload(30,8),(f1(30,8)*Atld(30)),(f2(30,8)*Atld(30))
elseif (j==9) then
write(110,*)time(k),prod(30,9),loss(30,9),grav(30,9),coll(30,9),epres(30,9),mload(30,9),(f1(30,9)*Atld(30)),(f2(30,9)*Atld(30))
endif

!!!!!!Step 1a: Find w1m and w2m (half grid step I-1/2, half time step k+1/2)
	DO i=2,NGrid
		w1m(i,j)=0.5*(w1(i-1,j)+w1(i,j))-0.5*DTDS*(F1(i,j)-F1(i-1,j))+&
			0.5*(R1(i,j)+R1(i-1,j))*DELT
		w2m(i,j)=0.5*(w2(i-1,j)+w2(i,j))-0.5*DTDS*(F2(i,j)-F2(i-1,j))+&
			0.5*(R2(i,j)+R2(i-1,j))*DELT
		fhm1(i,j)=w2m(i,j)
		fhm2(i,j)=w2m(i,j)*w2m(i,j)/w1m(i,j)+CSS0(i,j)*w1m(i,j)
	enddo

!!!!!!Step 1b: Find w1p and w2p (half grid step I+1/2, half time step k+1/2)
	DO i=1,NGrid1
 		w1p(i,j)=0.5*(w1(i+1,j)+w1(i,j))-0.5*DTDS*(F1(i+1,j)-F1(i,j))+&
			0.5*(R1(i,j)+R1(i+1,j))*DELT
  		w2p(i,j)=0.5*(w2(i+1,j)+w2(i,j))-0.5*DTDS*(F2(i+1,j)-F2(i,j))+&
			0.5*(R2(i,j)+R2(i+1,j))*DELT
		fhp1(i,j)=w2p(i,j)
		fhp2(i,j)=w2p(i,j)*w2p(i,j)/w1p(i,j)+CSS0(i,j)*w1p(i,j)
	enddo

!!!!!!!STEP 2, Find w1new and w2new (values at I, full time step k to k+1)
	DO i=2,NGrid1
		w1new(i,j)=w1(i,j)-DTDS*(fhp1(i,j)-fhm1(i,j))+R1(i,j)*DELT+&
			DDEN*DTDS/DELS*(w1(i+1,j)+w1(i-1,j)-2*w1(i,j))
		w2new(i,j)=w2(i,j)-DTDS*(fhp2(i,j)-fhm2(i,j))+R2(i,j)*DELT+&
			DVEL*DTDS/DELS*(w2(i+1,j)+w2(i-1,j)-2*w2(i,j))
	enddo

!!!!!!! done with ion loop     
!      write(*,*) 'Finished looping through Ion:',j
   enddo

!update all w1/w2 values
   do j=1,NIon
      do i=2,NGrid1
		w1(i,j)=w1new(i,j)
		w2(i,j)=w2new(i,j)
      enddo
   enddo
!write(*,*) "w1 & w2 values updated"

!Update electron densities
   DO i=1,NGrid
      NElec(i)=0.0
      do j=1,NIon
		NElec(i)=NElec(i)+w1(i,j)
      enddo
   enddo
!write(*,*) "Electron density updated"

!!!!!!! write within time step
   DO i=1,NGrid
      do j=1,NIon
		den=w1(i,j)/A(i)
		uu=w2(i,j)/A(i)/den
		DNN(K,i,j)=den
		UUU(K,i,j)=uu
		elecden(K,i)=NElec(i)/A(i)
      enddo
   enddo

!!!Writing out w terms to init_conditions
   if ((k==(1+(NSteps-1)/4)).or.(k==(1+(NSteps-1)/2)).or.(k==(1+3*(NSteps-1)/4)).or.(k==NSteps)) then
	!check if any values are NaN/infinity
	!make sure GOTO(#) is same as line# for enddo after write statement
	do i=1,NGrid
	   do j=1,NIon
		infinity = huge(0.00)
		if (w1(i,j)/=w1(i,j)) then
			GOTO 328
		elseif (w1(i,j)>infinity) then
			GOTO 328
		elseif (w2(i,j)/=w2(i,j)) then
			GOTO 328
		elseif (w2(i,j)>infinity) then
			GOTO 328
		endif
	   enddo
	enddo

   !write out w terms to new file
	open(30, file='init_conditions_time.txt')
	   write(30,*)time(k)
	close(30)
	open(31, file='init_conditions_w1.txt')
	do i=1,NGrid
	   write(31,*)i,w1(i,1),w1(i,2),w1(i,3),w1(i,4),w1(i,5),w1(i,6),w1(i,7),w1(i,8),w1(i,9)
	enddo
	close(31)
	open(32, file='init_conditions_w2.txt')
	do i=1,NGrid
	   write(32,*)i,w2(i,1),w2(i,2),w2(i,3),w2(i,4),w2(i,5),w2(i,6),w2(i,7),w2(i,8),w2(i,9)
	enddo
	close(32)
	328 continue
   endif

!!!!!!! done with time-loops
write(*,*) "Done with time step",k
enddo

!!check s grid is correct
!do i=1,NGrid
!	write(*,*) i,s(i)
!enddo

!!!!!!Writing out new ion densities and velocities
write(40,*) 'O2+ Densities: [1/m^3]'
write(40,*) '          i   s[m]             alt[km]          t=0.0         ',(k*DELT,k=100000,NSteps,100000)
write(41,*) 'O+ Densities: [1/m^3]'
write(41,*) '          i   s[m]             alt[km]          t=0.0         ',(k*DELT,k=100000,NSteps,100000)
write(42,*) 'CO2+ Densities: [1/m^3]'
write(42,*) '          i   s[m]             alt[km]          t=0.0         ',(k*DELT,k=100000,NSteps,100000)
write(43,*) 'CO+ Densities: [1/m^3]'
write(43,*) '          i   s[m]             alt[km]          t=0.0         ',(k*DELT,k=100000,NSteps,100000)
write(44,*) 'NO+ Densities: [1/m^3]'
write(44,*) '          i   s[m]             alt[km]          t=0.0         ',(k*DELT,k=100000,NSteps,100000)
write(45,*) 'N2+ Densities: [1/m^3]'
write(45,*) '          i   s[m]             alt[km]          t=0.0         ',(k*DELT,k=100000,NSteps,100000)
write(46,*) 'OH+ Densities: [1/m^3]'
write(46,*) '          i   s[m]             alt[km]          t=0.0         ',(k*DELT,k=100000,NSteps,100000)
write(47,*) 'H+ Densities: [1/m^3]'
write(47,*) '          i   s[m]             alt[km]          t=0.0         ',(k*DELT,k=100000,NSteps,100000)
write(48,*) 'HCO+ Densities: [1/m^3]'
write(48,*) '          i   s[m]             alt[km]          t=0.0         ',(k*DELT,k=100000,NSteps,100000)
write(50,*) 'O2+ Velocities: [m/s]'
write(50,*) '          i   s[m]             alt[km]          t=0.0         ',(k*DELT,k=100000,NSteps,100000)
write(51,*) 'O+ Velocities: [m/s]'
write(51,*) '          i   s[m]             alt[km]          t=0.0         ',(k*DELT,k=100000,NSteps,100000)
write(52,*) 'CO2+ Velocities: [m/s]'
write(52,*) '          i   s[m]             alt[km]          t=0.0         ',(k*DELT,k=100000,NSteps,100000)
write(53,*) 'CO+ Velocities: [m/s]'
write(53,*) '          i   s[m]             alt[km]          t=0.0         ',(k*DELT,k=100000,NSteps,100000)
write(54,*) 'NO+ Velocities: [m/s]'
write(54,*) '          i   s[m]             alt[km]          t=0.0         ',(k*DELT,k=100000,NSteps,100000)
write(55,*) 'N2+ Velocities: [m/s]'
write(55,*) '          i   s[m]             alt[km]          t=0.0         ',(k*DELT,k=100000,NSteps,100000)
write(56,*) 'OH+ Velocities: [m/s]'
write(56,*) '          i   s[m]             alt[km]          t=0.0         ',(k*DELT,k=100000,NSteps,100000)
write(57,*) 'H+ Velocities: [m/s]'
write(57,*) '          i   s[m]             alt[km]          t=0.0         ',(k*DELT,k=100000,NSteps,100000)
write(58,*) 'HCO+ Velocities: [m/s]'
write(58,*) '          i   s[m]             alt[km]          t=0.0         ',(k*DELT,k=100000,NSteps,100000)
write(60,*) 'e- Densities: [1/m^3]'
write(60,*) '          i   s[m]             alt[km]          t=0.0         ',(k*DELT,k=100000,NSteps,100000)
DO i=1,NGrid
write(40,*)i,s(i),alt(i),den_i(i,1),(dnn(k,i,1),k=100000,NSteps,100000)
write(41,*)i,s(i),alt(i),den_i(i,2),(dnn(k,i,2),k=100000,NSteps,100000)
write(42,*)i,s(i),alt(i),den_i(i,3),(dnn(k,i,3),k=100000,NSteps,100000)
write(43,*)i,s(i),alt(i),den_i(i,4),(dnn(k,i,4),k=100000,NSteps,100000)
write(44,*)i,s(i),alt(i),den_i(i,5),(dnn(k,i,5),k=100000,NSteps,100000)
write(45,*)i,s(i),alt(i),den_i(i,6),(dnn(k,i,6),k=100000,NSteps,100000)
write(46,*)i,s(i),alt(i),den_i(i,7),(dnn(k,i,7),k=100000,NSteps,100000)
write(47,*)i,s(i),alt(i),den_i(i,8),(dnn(k,i,8),k=100000,NSteps,100000)
write(48,*)i,s(i),alt(i),den_i(i,9),(dnn(k,i,9),k=100000,NSteps,100000)
write(50,*)i,s(i),alt(i),uu_i(i,1),(uuu(k,i,1),k=100000,NSteps,100000)
write(51,*)i,s(i),alt(i),uu_i(i,2),(uuu(k,i,2),k=100000,NSteps,100000)
write(52,*)i,s(i),alt(i),uu_i(i,3),(uuu(k,i,3),k=100000,NSteps,100000)
write(53,*)i,s(i),alt(i),uu_i(i,4),(uuu(k,i,4),k=100000,NSteps,100000)
write(54,*)i,s(i),alt(i),uu_i(i,5),(uuu(k,i,5),k=100000,NSteps,100000)
write(55,*)i,s(i),alt(i),uu_i(i,6),(uuu(k,i,6),k=100000,NSteps,100000)
write(56,*)i,s(i),alt(i),uu_i(i,7),(uuu(k,i,7),k=100000,NSteps,100000)
write(57,*)i,s(i),alt(i),uu_i(i,8),(uuu(k,i,8),k=100000,NSteps,100000)
write(58,*)i,s(i),alt(i),uu_i(i,9),(uuu(k,i,9),k=100000,NSteps,100000)
write(60,*)i,s(i),alt(i),elecden_i(i),(elecden(k,i),k=100000,NSteps,100000)
enddo

end program


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! write functions for R1/R2 here, which are RHS equations

!! R1=Production-Loss
!! =A*(ionization rate)*(nden)-alpha*(w1**2)/A
!! alpha is recombination rate for e- and ions
real function RHS1(j, nden, w1, f1, A, Atilde, PR, alpha)
	integer si
	parameter(si=202)
  integer i, j
	real irate, w1, f1, A, nden, Atilde, PR

	!PR = irate      !1/(cm^3*s), function of altitude
	!alpha  		 !cm^3/s, function of temperature

	!!test with R1=0
	R1=A*PR*1E6-alpha*1E-6*(w1**2)/A
return
end

!! R2=w1*gpar-w2*(bulkflowrate)+w1*dCSSO/ds+f2
real function RHS2(j, w1, gpar, w2, f2, Atilde, nden)
	parameter(si=863)
	integer j
	real w1, gpar, bfr, w2, f2, Atilde, nden

	bfr=1E-9   !cm^3/s
  !!test with R2=0
	R2=-w2*bfr*nden+w1*gpar
return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Define subroutines here

!subroutine Atilde(area, a0, L, R_d, lambda, atlde)
!! get A-tilde = (1/A)*(dA/ds)
!! function of s, but convert to a function of lambda using (dA/ds)=(dA/d[lambda])*(d[lambda]/ds)
!implicit none
!real, intent(in)  ::area, a0, L, R_d, lambda
!real, intent(out) ::atlde
!atlde = (6*a0*L**2)/(area*R_d)*(cos(lambda))**4*sin(lambda)*(1+3*(sin(lambda))**2)**(-1/2)
!return
!end

subroutine Atilde(a0, a1, a2, s1, s2, atlde)
!! get A-tilde = (1/A)*(dA/ds)
!! function of s, but convert to a function of lambda using (dA/ds)=(dA/d[lambda])*(d[lambda]/ds)

implicit none
real, intent(in)  ::a0, a1, a2, s1, s2
real, intent(out) ::atlde

atlde = (1/a0)*((a1-a2)/(s1-s2))

return
end

subroutine radius(lambda, L, R_d, radtotal)
!! get radius from dipole to field line, based on distance along field line
!! function of lambda
!! units of km
implicit none
real, intent(in)   ::lambda, L, R_d
real, intent(out)  ::radtotal

radtotal = L*R_d*(cos(lambda))**2

return
end


subroutine height(radtotal, R_d, z)
!! altitude
!! radtotal is a function of s
!! function of radtotal, thus overall a function of s
!! radtotal is distance from dipole to field line
!! R_d is depth of dipole under the surface
!! units of km
implicit none
real, intent(in)  ::radtotal, R_d
real, intent(out) ::z
z = radtotal - R_d
return
end


subroutine ion_temp(temp_i)
! ion temperature
! units of K
!!!!!!!!!!!!!!!!!!!make sure to come back and review this in the future
implicit none
real, intent(out)  ::temp_i

temp_i = 500

return
end

subroutine elec_temp(temp_e)
! electron temperature
! units of K
!!!!!!!!!!!!!!!!!make sure to come back and review this in the future
implicit none
real, intent(out) :: temp_e

temp_e = 1000

return
end

subroutine area(A0, x, y, A)
!! area of the field lines
!! function of s
!! related to B-field strength
!! units of km^2
implicit none
real, intent(in)   ::A0, x, y
real, intent(out)  ::A
real Rd, d

Rd = 260.5     ![km] Case 1 (max alt 400km)
!Rd = 254.9     ![km] Case 2 (max alt 600km)
!Rd = 252.8     ![km] Case 3 (max alt 800km)

if (x .LT. 250) then
	d = sqrt(x**2+(y+100)**2)
else
	d = sqrt((x-500)**2+(y+100)**2)
endif

A = A0*(d/Rd)**3

return
end


!subroutine gparallel(altitude, lambda, R_m, GPAR)
!! gravitational attraction parallel to the field lines
!! function of r, lambda
!! units of N/kg
!implicit none
!real, intent(in)   ::altitude, lambda, R_m
!real, intent(out)  ::GPAR
!real :: gsat, w, radtotal
!! constants
!gsat = 3.711            ![N/kg]
!w = 1.1324E-08          !rotation speed of Saturn   [rad/s]
!radtotal = altitude + R_m
!GPAR = w**2*radtotal*10**3*cos(lambda)*sin(lambda)-gsat*(R_m/radtotal)**2*2*sin(lambda)/(sqrt(1+3*(sin(lambda))**2))
!return
!end


subroutine gparallel(Rmars, alt, x, GPAR)
!! gravitational attraction parallel to the field lines
!! function of x/y location along field line
!! units of N/kg
implicit none
real, intent(in)   ::Rmars, alt, x
real, intent(out)  ::GPAR
real :: gmars, w, radtotal, slope, angle, sf

!! constants
gmars = 3.73356            ![N/kg]=[m/s^2]
w = 0.241                  !rotation speed of Mars  [rad/s]
radtotal = alt + Rmars     !distance from center of Mars [km]

!!calculate scaling factor (sf)
slope = -0.008*2*(x-250)
angle = atan(slope)
sf = sin(angle)

!GPAR = w**2*radtotal*10**3*cos(lambda)*sin(lambda)-gmars*(Rmars/radtotal)**2*sf
GPAR = (w**2.0/radtotal*10**3.0-gmars*(Rmars/radtotal)**2)*sf
return
end


subroutine linearint(s0, s1, l0, l1, sn, ln)
!! linear interpolater to get lambda(s) from s(lambda)
!! lambdas have units of radians, s have unis of km
implicit none
real, intent(in)  :: s0, s1, l0, l1, sn
real, intent(out) :: ln

ln = l0+(sn-s0)*(l1-l0)/(s1-s0)
return
end


subroutine linearint2(s0, s1, l0, l1, sn, ln)
!! linear interpolater to get lambda(s) from s(lambda)
!! lambdas have units of radians, s have unis of km
implicit none
real, intent(in)  :: s0, s1, l0, l1
real, intent(in) ::  sn
real :: logl0, logl1, logln
real, intent(out) :: ln

logl0 = log(l0)
logl1 = log(l1)

logln = logl0+(sn-s0)*(logl1-logl0)/(s1-s0)

ln = exp(logln)

return
end
