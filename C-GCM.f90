!****************************************************************************!
!                                                                            !
!                            GOECP 1.00 			             !
!                            						     !
!									     !
!                    Mathieu Reverdy / LMD / ESA                             !
!                      last update 11/19/2013                                !
!                                                                            !
!                                                                            !
!  1) Purpose : compute instantaneous and daily mean profiles of Scattering  !
!               Ration (SR), Color Ratio (CR) & Depolarization Ratio (DR)    !
!               over a model grid.                                           !
!                                                                            !
!               Compute various daily/monthly cloudiness over a model grid : !
!                  - Map of Low Mid High cloud cover                         !
!                  - 3D Cloud Fraction                                       !
!                  - 3D Cloud Phase                                          !
!                  - SR Histograms                                           !
!                                                                            !
!  2) Input   : SDS & META variables from CALIPSO hdf level1 Data files.     !
!                                                                            !
!  3) Output  : - 3D_CloudFraction :                                         !
!                    clcalipso(lon,lat,alt,time)                             !
!                    clrcalipso(lon,lat,alt,time)                            !
!                    uncalipso(lon,lat,alt,time)                             !
!                                                                            !
!               - Map of Low Mid High :                                      !
!                    cllcalipso(lon,lat,time)                                !
!                    clmcalipso(lon,lat,time)                                !
!                    clhcalipso(lon,lat,time)                                !
!                    cltcalipso(lon,lat,time)                                !
!                    clccalipso(lon,lat,time)                                !
!                                                                            !
!               - SR_Histo :                                                 !
!                    cfad_lidarsr532_Occ(lon,lat,alt,box,time)               !
!                    cfad_lidarsr532_Occ2(lon,lat,alt,box2,time)             !
!                                                                            !
!               - 3D_CloudPhase :                                            !
!                    ice_cloud(lon,lat,alt,time)                             !
!                    water_cloud(lon,lat,alt,time)                           !
!                                                                            !
!               - instant_SR_CR_DR                                           !
!                    longitude(it)                                           !
!                    latitude(it)                                            !
!                    altitude(it)                                            !
!                    time(it)                                                !
!                    SE(it)                                                  !
!                    instant_SR(it,alt)                                      !
!                    instant_CR(it,alt)                                      !
!                    instant_DR(it,alt)                                      !
!                                                                            !
!  4) Grid :                                                                 !
!               -CFMIP2   : 2° x 2° x 40levels from 0 to 19.2km              !
!                           (180,90,40)                                      !
!               -CFMIP1   : 1° x 1° x 40levels from 0 to 19.2km              !
!                           (360,180,40)                                     !
!               -CFMIP2.5 : 2.5° x 2.5° x 40levels from 0 to 19.2km          !
!                           (144,72,40)                                      !
!               -CFMIP    : 3.75° x 2.5° x 40levels from 0 to 19.2km         !
!                           (96,72,40)                                       !
!               -LMDZ     : 3.75° x 2.5° x 19levels from 0 to 40.4823km      !
!                           (96,72,19)                                       !
!               -LMDZ40   : 3.75° x 2.5° x 40levels from 0 to 40.4823km      !
!                           (96,72,40)                                       !
!               -NASA     : 5° x 5° x 40levels each 480m from 0 to 19.2km    !
!                           (73,37,41)                                       !
!                                                                            !
!  5) Compilation : use the makefile "makefiles.sh"                          !
!      makefiles.sh $1                                                       !
!      ifort $1.f90 -I/usr/include/hdf -L/usr/lib64/hdf -lmfhdf -ldf -ljpeg  !
!                   -lz -I/opt/netcdf-3.6.0-p1-ifort-64/include              !
!                   -L/opt/netcdf-3.6.0-p1-ifort-64/lib  -lnetcdf -o $1.e    !
!                                                                            !
!  6) Last updates & bug fix :                                               !
!   - 15/11/13 Draft 		                                             !
!									     !
!                                                                            !
!****************************************************************************!

!************************* SUBROUTINE SCHEME ***************************!
!                                                                       !
! 									!
!                                                                       !
!-----------------------------------------------------------------------!

program CGCM

	use netcdf
	implicit none
  
	integer i
	character*100 filename, rep, filename2,filename3,filename4,filename5,filename6,filename7
	real, dimension(:,:), allocatable :: CPRReflectivityValues, DoplerVelocityValues,CPRRefl,DoplerVelo,CPRRefl20, DoplerVelo20,CPRRefldbz, CPRReflmmh
	real, dimension(:), allocatable:: altValues, lidAltValues,lonValues
	integer numLon,numAlt,numzscene 
	real, dimension(:,:,:), allocatable :: pressureValues, temperatureValues
	real, dimension(:), allocatable :: zsceneValues,zzsceneValues,z480Values
	real, dimension(:), allocatable :: pressureValues2,temperatureValues2

rep="//bdd/CFMIP/workspace/mreverdy/ESA/C-GCM/Data/"
filename="radar_output.nc"
filename2="env_data_cosmo_II.nc"
filename3="testCPRdbz.nc"
filename4="testDopler.nc"
filename5="testCPR.nc"
filename6="testCPRmmh.nc"
filename7="testCPRnative.nc"

	call radarread(rep,filename,CPRReflectivityValues,DoplerVelocityValues,altValues,numLon,numAlt,lonValues)

	call environmentread(rep,filename2,pressureValues,temperatureValues,numzscene,zsceneValues)

	call verticaldimension(zzscenevalues,z480Values)

	call interp2D(CPRReflectivityValues,CPRRefl20,altValues,zzsceneValues)
	call interp2D(CPRRefl20,CPRRefl,zzsceneValues,z480Values)
		
	call interp2D(DoplerVelocityValues,DoplerVelo20,altValues,zzsceneValues)
	call interp2D(DoplerVelo20,DoplerVelo,zzsceneValues,z480Values)

	call refltodbz(CPRRefl,CPRRefldbz)
	call refltommh(CPRRefldbz,CPRReflmmh)
	
	call radarwrite2D(rep,filename7,CPRReflectivityValues,lonValues,altValues)

	call radarwrite2D(rep,filename5,CPRRefl,lonValues,z480Values)
	call radarwrite2D(rep,filename3,CPRRefldbz,lonValues,z480Values)
	call radarwrite2D(rep,filename6,CPRReflmmh,lonValues,z480Values)
	call radarwrite2D(rep,filename4,DoplerVelo,lonValues,z480Values)
contains

!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*                            END OF PROGRAM                                *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!



!****************************************************************************!
!****************************************************************************!
!**************************    SUBROUTINE LIST    ***************************!
!****************************************************************************!
!****************************************************************************!


!----------------------------------------------------------------------------!
subroutine radarread(rep,filename,CPRReflectivityValues,DoplerVelocityValues,altValues,numLon,numAlt,lonValues)

	use netcdf
	implicit none
	
	integer, dimension(nf90_max_var_dims) :: dimIDs
	character*100 filename, repdata
	character*100 rep
	integer lon, lat,alt, lidAlt, status, ncid
	integer numLat, numLidAlt
	integer numLon, numAlt
	real, dimension(:), allocatable :: lonValues, latValues
	real, dimension(:), allocatable:: altValues, lidAltValues, clidValues    	
	integer DoplerVelocity,CPRReflectivity
	real, dimension(:,:), allocatable :: DoplerVelocityValues, CPRReflectivityValues

repdata=trim(rep)//filename

print *,'Read radar netcdf files'
print *,trim(repdata)

status=nf90_open(repdata, nf90_nowrite, ncid)
status=nf90_inq_varid(ncid, "x_scene", lon)
status=nf90_inquire_variable(ncid, lon, dimids = dimIDs)
status=nf90_inquire_dimension(ncid, dimIDs(1), len = numLon)
allocate(lonValues(numLon))
status=nf90_get_var(ncid, lon, lonValues)

status=nf90_inq_varid(ncid, "y_scene", lat)
status=nf90_inquire_variable(ncid, lat, dimids = dimIDs)
status=nf90_inquire_dimension(ncid, dimIDs(1), len = numLat)
allocate(latValues(numLat))
status=nf90_get_var(ncid, lat, latValues)

status=nf90_inq_varid(ncid, "height", alt)
status=nf90_inquire_variable(ncid, alt, dimids = dimIDs)
status=nf90_inquire_dimension(ncid, dimIDs(1), len = numAlt)
allocate(altValues(numAlt))
status=nf90_get_var(ncid, alt, altValues)

status=nf90_inq_varid(ncid, "CPR_Reflectivity", CPRReflectivity)
status=nf90_inquire_variable(ncid, CPRReflectivity, dimids = dimIDs)
status=nf90_inquire_dimension(ncid, dimIDs(1), len = numLon)
status=nf90_inquire_dimension(ncid, dimIDs(2), len = numAlt)
allocate(CPRReflectivityValues(numLon,numAlt))
status=nf90_get_var(ncid, CPRReflectivity, CPRReflectivityValues)

status=nf90_inq_varid(ncid, "Dopler_Velocity", DoplerVelocity)
status=nf90_inquire_variable(ncid, DoplerVelocity, dimids = dimIDs)
status=nf90_inquire_dimension(ncid, dimIDs(1), len = numLon)
status=nf90_inquire_dimension(ncid, dimIDs(2), len = numAlt)
allocate(DoplerVelocityValues(numLon,numAlt))
status=nf90_get_var(ncid, DoplerVelocity, DoplerVelocityValues)

end subroutine radarread
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!
subroutine environmentread(rep,filename2,pressureValues,temperatureValues,numzscene,zsceneValues)

	use netcdf
	implicit none
	
 
	integer, dimension(nf90_max_var_dims) :: dimIDs
   	character*100 filename2, repdata, rep
	integer xscene, yscene,zscene, status, ncid2,pressure,temperature
	integer numyscene, numxscene, numzscene
	real, dimension(:), allocatable :: xsceneValues, ysceneValues, zsceneValues
	real, dimension(:,:,:), allocatable :: pressureValues, temperatureValues

repdata=trim(rep)//filename2

print *,'Read environment netcdf files'
print *,trim(repdata)

status=nf90_open(repdata, nf90_nowrite, ncid2)
status=nf90_inq_varid(ncid2, "x", xscene)
status=nf90_inquire_variable(ncid2, xscene, dimids = dimIDs)
status=nf90_inquire_dimension(ncid2, dimIDs(1), len = numxscene)
allocate(xsceneValues(numxscene))
status=nf90_get_var(ncid2, xscene, xsceneValues)

status=nf90_inq_varid(ncid2, "y", yscene)
status=nf90_inquire_variable(ncid2, yscene, dimids = dimIDs)
status=nf90_inquire_dimension(ncid2, dimIDs(1), len = numyscene)
allocate(ysceneValues(numyscene))
status=nf90_get_var(ncid2, yscene, ysceneValues)

status=nf90_inq_varid(ncid2, "z", zscene)
status=nf90_inquire_variable(ncid2, zscene, dimids = dimIDs)
status=nf90_inquire_dimension(ncid2, dimIDs(1), len = numzscene)
allocate(zsceneValues(numzscene))
status=nf90_get_var(ncid2, zscene, zsceneValues)

status=nf90_inq_varid(ncid2, "Pressure", pressure)
status=nf90_inquire_variable(ncid2, pressure, dimids = dimIDs)
status=nf90_inquire_dimension(ncid2, dimIDs(1), len = numxscene)
status=nf90_inquire_dimension(ncid2, dimIDs(2), len = numyscene)
status=nf90_inquire_dimension(ncid2, dimIDs(3), len = numzscene)
allocate(pressureValues(numxscene,numyscene,numzscene))
status=nf90_get_var(ncid2, pressure, pressureValues)

status=nf90_inq_varid(ncid2, "Temperature", temperature)
status=nf90_inquire_variable(ncid2, temperature, dimids = dimIDs)
status=nf90_inquire_dimension(ncid2, dimIDs(1), len = numxscene)
status=nf90_inquire_dimension(ncid2, dimIDs(2), len = numyscene)
status=nf90_inquire_dimension(ncid2, dimIDs(3), len = numzscene)
allocate(temperatureValues(numxscene,numyscene,numzscene))
status=nf90_get_var(ncid2, temperature, temperatureValues)

end subroutine environmentread
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!
subroutine verticaldimension(zzsceneValues,z480Values)

	integer i,pas
	real, dimension(:), allocatable :: zzsceneValues,z480Values
	real vertTemp

	120 format (A,E)
 	write(*, '(A)', advance="no") 'vertical resolution (km)? ='
	read(*,*) vertTemp


allocate(zzsceneValues(721)) 
	zzsceneValues(1)=0
	do i=2,721
	zzsceneValues(i)=zzsceneValues(i-1)+0.02
	end do

	pas=(12.96/vertTemp)+1

print *, pas
	allocate(z480Values(pas)) 
	z480Values(1)=0
	do i=2,pas
	z480Values(i)=z480Values(i-1)+vertTemp
	end do
end subroutine verticaldimension
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!
subroutine interp1D(var,var2,fvar,fvar2)
  implicit none

	real, dimension(:), allocatable :: var
	real, dimension(:), allocatable :: fvar, fvar2, var2
	integer i,j
	real a,b

print *,"data interpolation"
allocate(var2(size(fvar2)))

do i=2,size(fvar)
      a=(var(i)-var(i-1))/(fvar(i)-fvar(i-1))
      b=var(i)-a*fvar(i)

do j=1,size(fvar2)

if ((fvar2(j).ge.fvar(i-1)).and.(fvar2(j).lt.fvar(i)))then
      var2(j)=a*fvar2(j)+b

end if
end do
end do

end subroutine interp1D
!----------------------------------------------------------------------------!
!----------------------------------------------------------------------------!
subroutine interp2D(var,var2,fvar,fvar2)
  implicit none

	real, dimension(:,:), allocatable :: var, var2
	real, dimension(:), allocatable :: fvar, fvar2
	integer i,j,k
	real a,b

print *,"data interpolation"
allocate(var2(size(var,1),size(fvar2)))

do k=1,size(var,1)
do i=2,size(fvar)
      a=(var(k,i)-var(k,i-1))/(fvar(i)-fvar(i-1))
      b=var(k,i)-a*fvar(i)

do j=1,size(fvar2)

if ((fvar2(j).ge.fvar(i-1)).and.(fvar2(j).lt.fvar(i)))then
      var2(k,j)=a*fvar2(j)+b

end if
end do
end do
end do
end subroutine interp2D
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!
subroutine refltodbz(varin,varout)
  implicit none

	real, dimension(:,:), allocatable :: varin, varout
	
print *,"compute Reflectivity to dBZ"
allocate(varout(size(varin,1),size(varin,2)))

varout=10*log10(varin)

end subroutine refltodbz
!----------------------------------------------------------------------------!
!----------------------------------------------------------------------------!
subroutine refltommh(varin,varout)
  implicit none

	integer i,j
	real, dimension(:,:), allocatable :: varin, varout
	
print *,"compute Reflectivity to mmh"
allocate(varout(size(varin,1),size(varin,2)))

do i=1,size(varin,1)
do j=1,size(varin,2)

varout(i,j)=((10**(varin(i,j)/10))/200)**(0.625)

end do
end do
end subroutine refltommh
!----------------------------------------------------------------------------!



!----------------------------------------------------------------------------!
!subroutine check(istatus)
!
!	use netcdf
!	implicit none
!	
! 
!	integer, intent (in) :: istatus
!
!if (istatus /= nf90_noerr) then
!write (*,*) TRIM(ADJUSTL(nf90_strerror(istatus)))
!end if
!
!
!end subroutine check
!----------------------------------------------------------------------------!
!----------------------------------------------------------------------------!
subroutine radarwrite2D(rep,filenamewrite,var,xvar,yvar)

	use netcdf
	implicit none

	real, dimension(:), allocatable :: xvar,yvar
 	real, dimension(:,:), allocatable :: var
	integer, dimension(2) :: numatbid
   	character*100 filenamewrite, repdata, rep
	integer status, ncid,numx,numy
	integer numxid,numyid,varid,varid2,varid3

repdata=trim(rep)//filenamewrite

print *,'Write radar netcdf files'
print *,trim(repdata)

numx=size(xvar)
numy=size(yvar)

status=nf90_create(repdata, nf90_clobber, ncid)
status=nf90_def_dim(ncid,"x",numx,numxid)
status=nf90_def_dim(ncid,"y",numy,numyid)

numatbid=(/ numxid, numyid /)

status=nf90_def_var(ncid, "lon", nf90_float,numxid, varid)
status=nf90_def_var(ncid, "alt", nf90_float,numyid, varid2)
status=nf90_def_var(ncid,"CPR",nf90_float,numatbid,varid3)

status=nf90_enddef(ncid)

status=nf90_put_var(ncid,varid,xvar)
status=nf90_put_var(ncid,varid2,yvar)
status=nf90_put_var(ncid,varid3,var)

status=nf90_close(ncid)

end subroutine radarwrite2D
!----------------------------------------------------------------------------!
!----------------------------------------------------------------------------!
subroutine radarwrite1D(rep,filenamewrite,var,xvar)

	use netcdf
	implicit none

	real, dimension(:),allocatable :: xvar,var
   	character*100 filenamewrite, repdata, rep
	integer status, ncid,numx
	integer numxid,varid,varid2

repdata=trim(rep)//filenamewrite	


print *,'Write radar netcdf files'
print *,trim(repdata)
status=nf90_create(repdata, nf90_clobber, ncid)

numx=size(xvar)
status=nf90_def_dim(ncid,"x",numx,numxid)


status=nf90_def_var(ncid, "lon", NF90_FLOAT,numxid, varid)
status=nf90_def_var(ncid,"pmol",NF90_FLOAT,numxid,varid2)

status=nf90_enddef(ncid)


status=nf90_put_var(ncid,varid,xvar)
status=nf90_put_var(ncid,varid2,var)

status=nf90_close(ncid)


end subroutine radarwrite1D
!----------------------------------------------------------------------------!


!****************************************************************************!



end program CGCM

!****************************************************************************!
!                                                                            !
!                               END PROGRAM                                  !
!                                                                            !
!****************************************************************************!


