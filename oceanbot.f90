"""
Copyright 2015 ARC Centre of Excellence for Climate Systems Science

author: Aidan Heerdegen <aidan.heerdegen@anu.edu.au>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

program oceanbot

  ! Read in a MOM 3D ocean file and output corresponding 2D bottom velocities data

  use ncio, only: ncvar, nc_read, nc_open, nc_close, nc_create, nc_write_dim, nc_write, nc_get_att, &
       nc_size, nc_print_attr, nc_v_init
  use cmdline_arguments, only: get_options, bad_options, have_args, option_exists, &
      has_value, get_value, assignment(=), next_arg, num_args
  use iso_varying_string
  use file_functions, only: stderr
  use string_functions, only: join

  implicit none

  character(len=1000) :: datafile, gridfile, outfile

  real, allocatable :: kmu(:,:), u(:,:,:,:), v(:,:,:,:), ubot(:,:,:), vbot(:,:,:), hu(:,:), lon(:), lat(:)
  integer, allocatable :: time(:)

  type(ncvar) :: ncv

  integer :: i, j, k, id_grid, id_data, depth, nlon, nlat, nlvl, nmonth, error

  logical :: initialised = .false.

  type (varying_string) :: myoptions(2)

  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)
  myoptions(1) = 'help'
  myoptions(2) = 'outfile'

  ! This call parses the command line arguments for command line options
  call get_options(myoptions, error)

  ! Check we weren't passed duff options -- spit the dummy if we were
  if (error > 0) then
     write(stderr,*) 'ERROR! Unknown options: ',join(bad_options()," ")
     call usage
     STOP
  end if

  ! Check if we just want to print the usage
  if (option_exists('help')) then
     call usage
     STOP
  end if

  if (num_args() < 2) then
     write(stderr,*) 'ERROR! Must supply ocean grid file and at least one data file as command-line arguments'
     call usage
     STOP
  end if

  if (option_exists('outfile')) then
     ! We have specified variable name for vertical grid variable
     if (.NOT. has_value('outfile')) then
        write(stderr,*) 'Option outfile must specify a value!'
        call usage
        stop
     end if
     outfile = ""
     outfile = get_value('outfile')
  else
     outfile = "botvel.nc"
  end if

  ! Read in ocean grid file
  gridfile = next_arg()
  call nc_open(trim(gridfile), id_grid, writable=.false.)

  ! Initialize the netcdf variable info and load attributes
  call nc_v_init(ncv,trim("kmu"))
  call nc_get_att(id_grid,ncv,readmeta=.TRUE.)
  call nc_print_attr(ncv)
  
  nlon = ncv%dlen(1)
  nlat = ncv%dlen(2)
  
  print *,nlon,'x',nlat
  
  nmonth = num_args() - 1
     
  do k = 1, nmonth

     ! Read in ocean data
     datafile = next_arg()

     print *,'Reading ',k,trim(datafile)

     call nc_open(trim(datafile), id_data, writable=.false.)

     if (.not. initialised) then

        ! Initialize the netcdf variable info and load attributes
        call nc_v_init(ncv,trim("v"))
        call nc_get_att(id_data,ncv,readmeta=.TRUE.)
        call nc_print_attr(ncv)
        
        nlon = ncv%dlen(1)
        nlat = ncv%dlen(2)
        nlvl = ncv%dlen(3)
        
        allocate(kmu(nlon,nlat),u(nlon,nlat,nlvl,1),v(nlon,nlat,nlvl,1))
        allocate(ubot(nlon,nlat,nmonth),vbot(nlon,nlat,nmonth))
        allocate(time(nmonth))

        print *,'read kmu'
        call nc_read(trim(gridfile),"kmu",kmu,ncid=id_grid)

        allocate(lon(nlon),lat(nlat))

        call nc_read(trim(datafile),trim("xu_ocean"),lon(:),ncid=id_data)
        call nc_read(trim(datafile),trim("yu_ocean"),lat(:),ncid=id_data)

        ubot = -1.e20
        vbot = -1.e20

        initialised = .true.

     end if

     u = 0.; v = 0.
        
     ! Read longitude and latitude variables
     print *,'read u'
     call nc_read(trim(datafile),"u",u,ncid=id_data)
     print *,'read v'
     call nc_read(trim(datafile),"v",v,ncid=id_data)
     
     ! Read the time
     call nc_read(trim(datafile),trim("time"),time(k),ncid=id_data)

     do j = 1, nlat
        do i = 1, nlon
           if (kmu(i,j) < 0.) cycle
           depth = nint(kmu(i,j))
           if (depth <= 0) cycle
           if (u(i,j,depth,1) == -1.e20) cycle
           if (v(i,j,depth,1) == -1.e20) cycle
           ubot(i,j,k) = u(i,j,depth,1)
           vbot(i,j,k) = v(i,j,depth,1)
        end do
     end do

     call nc_close(id_data)

  end do
     
  deallocate(kmu,u,v)
  allocate(hu(nlon,nlat))
  print *,'read hu'
  call nc_read(trim(gridfile),"hu",hu,ncid=id_grid)

  call save(trim(outfile),ubot,vbot,hu,time,lon,lat)
     
contains
  
  subroutine save(outfile, ubot, vbot, depth, time, lon, lat)

    character(len=*), intent(in) :: outfile
    real, intent(in)             :: ubot(:,:,:), vbot(:,:,:), depth(:,:), lon(:), lat(:)
    integer, intent(in)          :: time(:)

    integer :: nx, ny

    nx = size(ubot,1)
    ny = size(ubot,2)

    print *,nx,ny

    ! Save result
    call nc_create(outfile,netcdf4=.TRUE.,overwrite=.TRUE.)
    call nc_write_dim(outfile,"xu_ocean",x=lon,long_name="ucell longitude",units="degrees_E")
    call nc_write_dim(outfile,"yu_ocean",x=lat,long_name="ucell latitude",units="degrees_N")
    call nc_write_dim(outfile,"time",x=time,long_name="time",units="days since 0001-01-01 00:00:00")
    call nc_write(outfile,"ubot",ubot,dim1="xu_ocean",dim2="yu_ocean",dim3="time",missing_value=-1.e20,long_name="bottom u velocity",units="m/s")
    call nc_write(outfile,"vbot",vbot,dim1="xu_ocean",dim2="yu_ocean",dim3="time",missing_value=-1.e20,long_name="bottom v velocity",units="m/s")
    call nc_write(outfile,"hu",depth,dim1="xu_ocean",dim2="yu_ocean",missing_value=-1.e20,long_name="ocean depth on u-cells",units="m")
  end subroutine save

 
  subroutine usage    
  
    write(stderr,*)
    write(stderr,*) 'Output 2D bottom velocity field from 3D data'
    write(stderr,*)
    write(stderr,*) 'Usage: oceanbot [--help] gridfile datafiles'
    write(stderr,*)
    write(stderr,*) '  --help      - print this message'
    write(stderr,*) '  --outfile   - name of output file (default: botvel.nc)'
    write(stderr,*)

  end subroutine usage


end program oceanbot

