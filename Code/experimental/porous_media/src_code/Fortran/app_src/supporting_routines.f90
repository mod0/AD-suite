! Author: Alexey Kuznetsov
! Modified: 28/12/2008
! this Fortran90 module contains a collection of subroutines for plotting data,
! including 2D, 3D plots, surfaces, polar coordinates, histograms
! it is a modification of the GNUFOR interface written by John Burkardt:
! http://orion.math.iastate.edu/burkardt/g_src/gnufor/gnufor.html 
!***********************************************************************************
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, version 3 of the License.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!    <http://www.gnu.org/licenses/>.
!***********************************************************************************
	module gnufor2
	implicit none
!***********************************************************************************
! these are default parameters which control linewidth, colors and terminal
!***********************************************************************************
	character(len=3), parameter	:: default_linewidth='1'
	character(len=100), parameter	:: default_color1='blue'
	character(len=100), parameter	:: default_color2='dark-green'
	character(len=100), parameter	:: default_color3='orange-red'
	character(len=100), parameter	:: default_color4='dark-salmon'
	character(len=100), parameter	:: default_terminal='wxt'
	character(len=100), parameter	:: default_palette='CMY'
!***********************************************************************************
	interface plot
		module procedure plot_1
		module procedure plot_2
		module procedure plot_3
		module procedure plot_4
	end interface
!***********************************************************************************
	interface surf
		module procedure surf_1
		module procedure surf_2
		module procedure surf_3
	end interface
!***********************************************************************************
	interface image
		module procedure image_1
		module procedure image_2
		module procedure image_3
		module procedure image_4
	end interface
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	contains
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	function my_date_and_time() result(f_result)
!***********************************************************************************
! this function creates a string with current date and time
! it is a default method to name output files
!***********************************************************************************
	implicit none
        character(len=8)  :: date
        character(len=10) :: time
        character(len=33) :: f_result
!***********************************************************************************
	call date_and_time(date,time)
	f_result= 'date_'//date(7:8)//'-'//date(5:6)//'-'//date(1:4)//'_time_'//time(1:2)//':'//time(3:4)//':'//time(5:10)
!***********************************************************************************
	end function my_date_and_time
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	function output_terminal(terminal) result(f_result)
!***********************************************************************************
	implicit none
	character(len=*),intent(in)	:: terminal
	integer, parameter		:: Nc=35
	character(len=Nc)		:: f_result
!***********************************************************************************
	select case(terminal)
		case('ps')
			f_result='postscript landscape color'
		case default
			f_result=terminal
	end select
!***********************************************************************************
	end function output_terminal
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine image_4(x,y,rgb,pause,terminal,filename,persist,input)
!***********************************************************************************
! this is the most general subroutine for generating 3D plots.
! The data is contained in a 3D array z(:,:,:)
!***********************************************************************************
	implicit none
	real(kind=4), intent(in)	:: x(:), y(:)
	integer,intent(in)		:: rgb(:,:,:)
	real(kind=4), optional		:: pause
	character(len=*),optional	:: terminal, filename, persist, input
	integer				:: nx, ny
	integer 			:: i, j, ierror, ios, file_unit
	character(len=100)		:: data_file_name, command_file_name, my_pause, my_persist,&
		& xrange1,xrange2,yrange1,yrange2
!***********************************************************************************
	nx=size(rgb(1,:,1))
	ny=size(rgb(1,1,:))
	if ((size(x).ne.nx).or.(size(y).ne.ny)) then
		ierror=1
		print *,'image2 ERROR: sizes of x(:),y(:) and gray(:,:) are not compatible'
		stop
	end if
	write (xrange1,'(e15.7)') minval(x)
	write (xrange2,'(e15.7)') maxval(x)
	write (yrange1,'(e15.7)') minval(y)
	write (yrange2,'(e15.7)') maxval(y)
	do j=1,ny
	do i=1,nx
		if ((maxval(rgb(:,i,j))>255).or.(minval(rgb(:,i,j))<0)) then
			print *,'image ERROR: a value of rgb(:,:,:) array is outside [0,255]'
			stop
		end if
	end do
	end do
!***********************************************************************************
	if (present(input)) then
		data_file_name='data_file_'//input//'.txt'
		command_file_name='command_file_'//input//'.txt'		
	else
		data_file_name='data_file.txt'
		command_file_name='command_file.txt'
	end if
!***********************************************************************************
	ierror=0	
	call get_unit(file_unit)	
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=data_file_name, status='replace', iostat=ios)	
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal data file.'
		stop
	end if
!***********************************************************************************	
! here we write the date to the data_file - the gnuplot will read this data later
!***********************************************************************************	
	do j=1,ny
		do i=1,nx
			write (file_unit,'(2E12.4,3I5)') x(i),y(j),rgb(1,i,j),rgb(2,i,j),rgb(3,i,j)
		end do
		write (file_unit,'(a)')
	end do
!***********************************************************************************	
	close (unit=file_unit)
!***********************************************************************************
	ierror = 0
	call get_unit(file_unit)
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=command_file_name, status='replace',	iostat=ios)
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal command file.'
		stop
	end if	
!***********************************************************************************
! here we write the commands to the commands file which gnuplot will execute
!***********************************************************************************
	my_persist='persist '
	if (present(persist).and.(persist=='no')) my_persist=' '
	if (present(terminal)) then
		write ( file_unit, '(a)' ) 'set terminal '// trim(output_terminal(terminal)) 
	if (present(filename)) then
		write ( file_unit, '(a)' ) 'set output "'// trim(filename) //'"' 
	else
		write ( file_unit, '(a)' ) 'set output "'//my_date_and_time()//'"' 
	end if
	else
		write ( file_unit, '(a)' ) 'set terminal ' // trim(default_terminal) // ' ' &
			& //trim(my_persist) // ' title  "Gnuplot"' 
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'set nokey'
	write ( file_unit, '(a)' ) 'set xrange ['// trim(xrange1) // ':'// trim(xrange2) //']'
	write ( file_unit, '(a)' ) 'set yrange ['// trim(yrange1) // ':'// trim(yrange2) //']'
	write ( file_unit, '(a)' ) 'unset colorbox'
!***********************************************************************************	
	write ( file_unit, '(a)' ) 'plot "' // trim ( data_file_name ) // &
	& '" with rgbimage'
!***********************************************************************************
	if (present(pause)) then
		if (pause<0.0) then
			write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
		else 
			write (	my_pause,'(e9.3)') pause
			write ( file_unit, '(a)' ) 'pause ' // trim(my_pause) 
		end if
	else
		write ( file_unit, '(a)' ) 'pause 0'
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'q'
	close ( unit = file_unit )
!***********************************************************************************
	call run_gnuplot (command_file_name)
!***********************************************************************************
	end subroutine image_4
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine image_3(rgb,pause,terminal,filename,persist,input)
!***********************************************************************************
! this is the most general subroutine for generating 3D plots.
! The data is contained in a 3D array z(:,:,:)
!***********************************************************************************
	implicit none
	integer, intent(in)		:: rgb(:,:,:)
	real(kind=4), optional		:: pause
	character(len=*),optional	:: terminal, filename, persist, input
	integer				:: nx, ny
	integer 			:: i, j, ierror, ios, file_unit
	character(len=100)		:: data_file_name, command_file_name, my_pause, my_persist
!***********************************************************************************
	nx=size(rgb(1,:,1))
	ny=size(rgb(1,1,:))
	do j=1,ny
	do i=1,nx
		if ((maxval(rgb(:,i,j))>255).or.(minval(rgb(:,i,j))<0)) then
			print *,'image ERROR: a value of rgb(:,:,:) array is outside [0,255]'
			stop
		end if
	end do
	end do
!***********************************************************************************
	if (present(input)) then
		data_file_name='data_file_'//input//'.txt'
		command_file_name='command_file_'//input//'.txt'		
	else
		data_file_name='data_file.txt'
		command_file_name='command_file.txt'
	end if
!***********************************************************************************
	ierror=0	
	call get_unit(file_unit)	
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=data_file_name, status='replace', iostat=ios)	
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal data file.'
		stop
	end if
!***********************************************************************************	
! here we write the date to the data_file - the gnuplot will read this data later
!***********************************************************************************	
	do j=1,ny
		do i=1,nx
			write (file_unit,'(5I5)') i,j,rgb(1,i,j),rgb(2,i,j),rgb(3,i,j)
		end do
		write (file_unit,'(a)')
	end do
!***********************************************************************************	
	close (unit=file_unit)
!***********************************************************************************
	ierror = 0
	call get_unit(file_unit)
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=command_file_name, status='replace',	iostat=ios)
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal command file.'
		stop
	end if	
!***********************************************************************************
! here we write the commands to the commands file which gnuplot will execute
!***********************************************************************************
	my_persist='persist '
	if (present(persist).and.(persist=='no')) my_persist=' '
	if (present(terminal)) then
		write ( file_unit, '(a)' ) 'set terminal '// trim(output_terminal(terminal)) 
	if (present(filename)) then
		write ( file_unit, '(a)' ) 'set output "'// trim(filename) //'"' 
	else
		write ( file_unit, '(a)' ) 'set output "'//my_date_and_time()//'"' 
	end if
	else
		write ( file_unit, '(a)' ) 'set terminal ' // trim(default_terminal) // ' ' &
			& //trim(my_persist) // ' title  "Gnuplot"' 
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'set nokey'
	write ( file_unit, '(a)' ) 'unset border'
	write ( file_unit, '(a)' ) 'unset xtics'
	write ( file_unit, '(a)' ) 'unset ytics'
	write ( file_unit, '(a)' ) 'unset colorbox'
!***********************************************************************************	
	write ( file_unit, '(a)' ) 'plot "' // trim ( data_file_name ) // &
	& '" with rgbimage'
!***********************************************************************************
	if (present(pause)) then
		if (pause<0.0) then
			write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
		else 
			write (	my_pause,'(e9.3)') pause
			write ( file_unit, '(a)' ) 'pause ' // trim(my_pause) 
		end if
	else
		write ( file_unit, '(a)' ) 'pause 0'
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'q'
	close ( unit = file_unit )
!***********************************************************************************
	call run_gnuplot (command_file_name)
!***********************************************************************************
	end subroutine image_3
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine image_2(x,y,gray,pause,palette,terminal,filename,persist,input)
!***********************************************************************************
! this is the most general subroutine for generating 3D plots.
! The data is contained in a 3D array z(:,:,:)
!***********************************************************************************
	implicit none
	real(kind=4), intent(in)	:: x(:), y(:), gray(:,:)
	real(kind=4), optional		:: pause
	character(len=*),optional	:: palette, terminal, filename, persist, input
	integer				:: nx, ny
	integer 			:: i, j, ierror, ios, file_unit
	character(len=100)		:: data_file_name, command_file_name, my_pause, my_persist,&
		& xrange1,xrange2,yrange1,yrange2
!***********************************************************************************
	nx=size(gray(:,1))
	ny=size(gray(1,:))
	if ((size(x).ne.nx).or.(size(y).ne.ny)) then
		ierror=1
		print *,'image2 ERROR: sizes of x(:),y(:) and gray(:,:) are not compatible'
		stop
	end if
	write (xrange1,'(e15.7)') minval(x)
	write (xrange2,'(e15.7)') maxval(x)
	write (yrange1,'(e15.7)') minval(y)
	write (yrange2,'(e15.7)') maxval(y)
!***********************************************************************************
	if (present(input)) then
		data_file_name='data_file_'//input//'.txt'
		command_file_name='command_file_'//input//'.txt'		
	else
		data_file_name='data_file.txt'
		command_file_name='command_file.txt'
	end if
!***********************************************************************************
	ierror=0	
	call get_unit(file_unit)	
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=data_file_name, status='replace', iostat=ios)	
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal data file.'
		stop
	end if
!***********************************************************************************	
! here we write the date to the data_file - the gnuplot will read this data later
!***********************************************************************************	
	do j=1,ny
		do i=1,nx
			write (file_unit,'(3E12.4)') x(i), y(j), gray(i,j)
		end do
		write (file_unit,'(a)')
	end do
!***********************************************************************************	
	close (unit=file_unit)
!***********************************************************************************
	ierror = 0
	call get_unit(file_unit)
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=command_file_name, status='replace',	iostat=ios)
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal command file.'
		stop
	end if	
!***********************************************************************************
! here we write the commands to the commands file which gnuplot will execute
!***********************************************************************************
	my_persist='persist '
	if (present(persist).and.(persist=='no')) my_persist=' '
	if (present(terminal)) then
		write ( file_unit, '(a)' ) 'set terminal '// trim(output_terminal(terminal)) 
	if (present(filename)) then
		write ( file_unit, '(a)' ) 'set output "'// trim(filename) //'"' 
	else
		write ( file_unit, '(a)' ) 'set output "'//my_date_and_time()//'"' 
	end if
	else
		write ( file_unit, '(a)' ) 'set terminal ' // trim(default_terminal) // ' ' &
			& //trim(my_persist) // ' title  "Gnuplot"' 
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'set nokey'
	write ( file_unit, '(a)' ) 'set xrange ['// trim(xrange1) // ':'// trim(xrange2) //']'
	write ( file_unit, '(a)' ) 'set yrange ['// trim(yrange1) // ':'// trim(yrange2) //']'
	write ( file_unit, '(a)' ) 'unset colorbox'
	if (present(palette)) then
	if ((trim(palette).ne.'RGB').and.(trim(palette).ne.'HSV').and.(trim(palette).ne.'CMY').and.&
		& (trim(palette).ne.'YIQ').and.(trim(palette).ne.'XYZ')) then
		write ( file_unit, '(a)' ) 'set palette '// trim(palette)
	else
		write ( file_unit, '(a)' ) 'set palette model '// trim(palette)
	end if
	else
		write ( file_unit, '(a)' ) 'set palette model '// trim(default_palette)
	end if
!***********************************************************************************	
	write ( file_unit, '(a)' ) 'plot "' // trim ( data_file_name ) // &
	& '" with image'
!***********************************************************************************
	if (present(pause)) then
		if (pause<0.0) then
			write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
		else 
			write (	my_pause,'(e9.3)') pause
			write ( file_unit, '(a)' ) 'pause ' // trim(my_pause) 
		end if
	else
		write ( file_unit, '(a)' ) 'pause 0'
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'q'
	close ( unit = file_unit )
!***********************************************************************************
	call run_gnuplot (command_file_name)
!***********************************************************************************
	end subroutine image_2
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine image_1(gray,pause,palette,terminal,filename,persist,input)
!***********************************************************************************
! this is the most general subroutine for generating 3D plots.
! The data is contained in a 3D array z(:,:,:)
!***********************************************************************************
	implicit none
	real(kind=4), intent(in)	:: gray(:,:)
	real(kind=4), optional		:: pause
	character(len=*),optional	:: palette, terminal, filename, persist, input
	integer				:: nx, ny
	integer 			:: i, j, ierror, ios, file_unit
	character(len=100)		:: data_file_name, command_file_name, my_pause, my_persist
!***********************************************************************************
	nx=size(gray(:,1))
	ny=size(gray(1,:))
!***********************************************************************************
	if (present(input)) then
		data_file_name='data_file_'//input//'.txt'
		command_file_name='command_file_'//input//'.txt'		
	else
		data_file_name='data_file.txt'
		command_file_name='command_file.txt'
	end if
!***********************************************************************************
	ierror=0	
	call get_unit(file_unit)	
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=data_file_name, status='replace', iostat=ios)	
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal data file.'
		stop
	end if
!***********************************************************************************	
! here we write the date to the data_file - the gnuplot will read this data later
!***********************************************************************************	
	do j=1,ny
		do i=1,nx
			write (file_unit,'(I5,I5,E15.7)') i,j,gray(i,j)
		end do
		write (file_unit,'(a)')
	end do
!***********************************************************************************	
	close (unit=file_unit)
!***********************************************************************************
	ierror = 0
	call get_unit(file_unit)
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=command_file_name, status='replace',	iostat=ios)
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal command file.'
		stop
	end if	
!***********************************************************************************
! here we write the commands to the commands file which gnuplot will execute
!***********************************************************************************
	my_persist='persist '
	if (present(persist).and.(persist=='no')) my_persist=' '
	if (present(terminal)) then
		write ( file_unit, '(a)' ) 'set terminal '// trim(output_terminal(terminal)) 
	if (present(filename)) then
		write ( file_unit, '(a)' ) 'set output "'// trim(filename) //'"' 
	else
		write ( file_unit, '(a)' ) 'set output "'//my_date_and_time()//'"' 
	end if
	else
		write ( file_unit, '(a)' ) 'set terminal ' // trim(default_terminal) // ' ' &
			& //trim(my_persist) // ' title  "Gnuplot"' 
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'set nokey'
	write ( file_unit, '(a)' ) 'unset border'
	write ( file_unit, '(a)' ) 'unset xtics'
	write ( file_unit, '(a)' ) 'unset ytics'
	write ( file_unit, '(a)' ) 'unset colorbox'
	if (present(palette)) then
	if ((trim(palette).ne.'RGB').and.(trim(palette).ne.'HSV').and.(trim(palette).ne.'CMY').and.&
		& (trim(palette).ne.'YIQ').and.(trim(palette).ne.'XYZ')) then
		write ( file_unit, '(a)' ) 'set palette '// trim(palette)
	else
		write ( file_unit, '(a)' ) 'set palette model '// trim(palette)
	end if
	else
		write ( file_unit, '(a)' ) 'set palette model '// trim(default_palette)
	end if
!***********************************************************************************	
	write ( file_unit, '(a)' ) 'plot "' // trim ( data_file_name ) // &
	& '" with image'
!***********************************************************************************
	if (present(pause)) then
		if (pause<0.0) then
			write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
		else 
			write (	my_pause,'(e9.3)') pause
			write ( file_unit, '(a)' ) 'pause ' // trim(my_pause) 
		end if
	else
		write ( file_unit, '(a)' ) 'pause 0'
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'q'
	close ( unit = file_unit )
!***********************************************************************************
	call run_gnuplot (command_file_name)
!***********************************************************************************
	end subroutine image_1
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine plot3d(x,y,z,pause,color,terminal,filename,persist,input,linewidth)
!***********************************************************************************
! this subroutine plots 3D curve, given by three arrays x,y,z
!***********************************************************************************
	implicit none
	real(kind=8), intent(in)	:: x(:),y(:),z(:)
	real(kind=4), optional		:: pause, linewidth
	character(len=*),optional	:: color, terminal, filename, persist, input
	integer 			:: i, j, ierror, ios, file_unit, nx
	character(len=100)		:: data_file_name, command_file_name, my_color, my_pause, my_persist, my_linewidth
!***********************************************************************************
! prepare the data
	nx=size(x)
	if ((size(x).ne.size(y)).or.(size(x).ne.size(z))) then
		print *,'subroutine plot3d ERROR: incompatible sizes of x(:),y(:) and z(:)'
		stop
	end if
!***********************************************************************************
	if (present(input)) then
		data_file_name='data_file_'//input//'.txt'
		command_file_name='command_file_'//input//'.txt'		
	else
		data_file_name='data_file.txt'
		command_file_name='command_file.txt'
	end if
!***********************************************************************************
	ierror=0	
	call get_unit(file_unit)	
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=data_file_name, status='replace', iostat=ios)	
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal data file.'
		stop
	end if
!***********************************************************************************	
! here we write the date to the data_file - the gnuplot will read this data later
!***********************************************************************************	
	do i=1,nx
		write (file_unit,'(3E15.7)') x(i), y(i), z(i)
	end do
!***********************************************************************************	
	close (unit=file_unit)
!***********************************************************************************
	ierror = 0
	call get_unit(file_unit)
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=command_file_name, status='replace',	iostat=ios)
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal command file.'
		stop
	end if	
!***********************************************************************************
! here we write the commands to the commands file which gnuplot will execute
!***********************************************************************************
	my_persist='persist'
	if (present(persist).and.(persist=='no')) my_persist=' '
	if (present(terminal)) then
		write ( file_unit, '(a)' ) 'set terminal '// trim(output_terminal(terminal)) 
	if (present(filename)) then
		write ( file_unit, '(a)' ) 'set output "'// trim(filename) //'"' 
	else
		write ( file_unit, '(a)' ) 'set output "'//my_date_and_time()//'"' 
	end if
	else
		write ( file_unit, '(a)' ) 'set terminal ' // trim(default_terminal) // ' ' &
		&// trim(my_persist) // ' title "Gnuplot"' 
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'set nokey'
!***********************************************************************************
	if (present(linewidth)) then
		write (	my_linewidth,'(e9.3)') linewidth
	else
		my_linewidth=trim(default_linewidth)
	end if	
	if (present(color)) then
		my_color='"'//trim(color)//'"'
	else
		my_color='"'//trim(default_color1)//'"'
	end if	
	write ( file_unit, '(a)' ) 'splot "' // trim ( data_file_name ) // &
	'" using 1:2:3 with lines linecolor rgb' // my_color //' linewidth ' // my_linewidth
!***********************************************************************************
	if (present(pause)) then
		if (pause<0.0) then
			write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
		else 
			write (	my_pause,'(e9.3)') pause
			write ( file_unit, '(a)' ) 'pause ' // trim(my_pause) 
		end if
	else
		write ( file_unit, '(a)' ) 'pause 0'
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'q'
	close ( unit = file_unit )
!***********************************************************************************
	call run_gnuplot (command_file_name)
!***********************************************************************************
	end subroutine plot3d
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine hist(x,n,pause,color,terminal,filename,persist,input)
!***********************************************************************************
! this subroutine plots the histogram of data contained in array x, using n bins
!***********************************************************************************
	implicit none
	real(kind=8), intent(in)	:: x(:) !the data to plot
	integer, intent(in)		:: n !the number of intervals
	real(kind=4), optional		:: pause
	character(len=*),optional	:: color, terminal, filename, persist, input
	integer 			:: i, j, ierror, ios, file_unit, nx
	character(len=100)		:: data_file_name, command_file_name, yrange, xrange1, xrange2, my_color, &
				& xtic_start, dxtic, xtic_end, my_pause, my_persist
	real(kind=8)			:: xmin, xmax, xhist(0:n), yhist(n+1), dx
!***********************************************************************************
! prepare the data
	nx=size(x)
	xmin=minval(x)
	xmax=maxval(x)
	dx=(xmax-xmin)/n
	do i=0,n
		xhist(i)=xmin+i*dx
	end do
	yhist=0.0d0
	do i=1,nx
		j=floor((x(i)-xmin)/dx)+1
		yhist(j)=yhist(j)+1
	end do
!***********************************************************************************
	write (dxtic,'(e15.7)') dx
	write (yrange,'(e15.7)') maxval(yhist)*1.2
	write (xrange1,'(e15.7)') xmin-(n/10.0)*dx
	write (xrange2,'(e15.7)') xmax+(n/10.0)*dx
	xtic_start=xrange1
	xtic_end=xrange2
!***********************************************************************************
	if (present(input)) then
		data_file_name='data_file_'//input//'.txt'
		command_file_name='command_file_'//input//'.txt'		
	else
		data_file_name='data_file.txt'
		command_file_name='command_file.txt'
	end if
!***********************************************************************************
	ierror=0	
	call get_unit(file_unit)	
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=data_file_name, status='replace', iostat=ios)	
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal data file.'
		stop
	end if
!***********************************************************************************	
! here we write the date to the data_file - the gnuplot will read this data later
!***********************************************************************************	
	do i=1,n
		write (file_unit,'(2E15.7)') (xhist(i-1)+0.5*dx), yhist(i)
	end do
!***********************************************************************************	
	close (unit=file_unit)
!***********************************************************************************
	ierror = 0
	call get_unit(file_unit)
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=command_file_name, status='replace',	iostat=ios)
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal command file.'
		stop
	end if	
!***********************************************************************************
! here we write the commands to the commands file which gnuplot will execute
!***********************************************************************************
	my_persist='persist'
	if (present(persist).and.(persist=='no')) my_persist=' '
	if (present(terminal)) then
		write ( file_unit, '(a)' ) 'set terminal '// trim(output_terminal(terminal)) 
	if (present(filename)) then
		write ( file_unit, '(a)' ) 'set output "'// trim(filename) //'"' 
	else
		write ( file_unit, '(a)' ) 'set output "'//my_date_and_time()//'"' 
	end if
	else
		write ( file_unit, '(a)' ) 'set terminal ' // trim(default_terminal) // ' ' &
			& //trim(my_persist) // ' title  "Gnuplot"' 
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'set nokey'
	write ( file_unit, '(a)' ) 'set yrange [0.0:'// trim(yrange) //']'
	write ( file_unit, '(a)' ) 'set xrange ['// trim(xrange1) // ':'// trim(xrange2) //']'
	write ( file_unit, '(a)' ) 'set xtic nomirror rotate by -45 '
	write ( file_unit, '(a)' ) 'set xtics '// trim(xtic_start) // ','// trim(dxtic) // ','// trim(xtic_end)
	write ( file_unit, '(a)' ) 'set style data histograms'
	write ( file_unit, '(a)' ) 'set style fill solid border -1'
!***********************************************************************************
	if (present(color)) then
		my_color='"'//color//'"'
	else
		my_color='"'//trim(default_color1)//'"'
	end if	
	write ( file_unit, '(a)' ) 'plot "' // trim ( data_file_name ) // &
	'" using 1:2 with boxes linecolor rgb' // trim(my_color)
!***********************************************************************************
	if (present(pause)) then
		if (pause<0.0) then
			write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
		else 
			write (	my_pause,'(e9.3)') pause
			write ( file_unit, '(a)' ) 'pause ' // trim(my_pause) 
		end if
	else
		write ( file_unit, '(a)' ) 'pause 0'
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'q'
	close ( unit = file_unit )
!***********************************************************************************
	call run_gnuplot (command_file_name)
!***********************************************************************************
	end subroutine hist
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine surf_3(x,y,z,pause,palette,terminal,filename,pm3d,contour,persist,input)
!***********************************************************************************
! this subroutine plots a surface. x and y are arrays needed to generate the x-y grid
! z(:,:) is a 2D array
!***********************************************************************************
	implicit none
	real(kind=8), intent(in)	:: x(:),y(:),z(:,:)
	real(kind=4), optional		:: pause
	real(kind=8)			:: xyz(3,size(z(:,1)),size(z(1,:)))
	character(len=*),optional	:: palette, terminal, filename, pm3d, contour, persist, input
	integer				:: nx, ny
	integer 			:: i, j
!***********************************************************************************
	nx=size(z(:,1))
	ny=size(z(1,:))
	if ((size(x).ne.nx).or.(size(y).ne.ny)) then
		print *,'subroutine surf_3 ERROR: sizes of x(:),y(:), and z(:,:) are incompatible'
		stop
	end if 
!***********************************************************************************
	do i=1,nx
		do j=1,ny
			xyz(1,i,j)=x(i)
			xyz(2,i,j)=y(j)
			xyz(3,i,j)=z(i,j)
		end do
	end do
	call surf_1(xyz,pause,palette,terminal,filename,pm3d,contour,persist,input)
!***********************************************************************************	
	end subroutine surf_3
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine surf_2(z,pause,palette,terminal,filename,pm3d,contour,persist,input)
!***********************************************************************************
! this subroutine plots a surface. The only input is a 2D array z(:,:), the x-y grid 
! is generated automatically
!***********************************************************************************
	implicit none
	real(kind=8), intent(in)	:: z(:,:)
	real(kind=4), optional		:: pause
	real(kind=8)			:: xyz(3,size(z(:,1)),size(z(1,:)))
	character(len=*),optional	:: palette, terminal, filename, pm3d, contour, persist, input
	integer				:: nx, ny
	integer 			:: i, j
!***********************************************************************************
	nx=size(z(:,1))
	ny=size(z(1,:))
	do i=1,nx
		do j=1,ny
			xyz(1,i,j)=dble(i)
			xyz(2,i,j)=dble(j)
			xyz(3,i,j)=z(i,j)
		end do
	end do
	call surf_1(xyz,pause,palette,terminal,filename,pm3d,contour,persist,input)
!***********************************************************************************	
	end subroutine surf_2
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine surf_1(xyz,pause,palette,terminal,filename,pm3d,contour,persist,input)
!***********************************************************************************
! this is the most general subroutine for generating 3D plots.
! The data is contained in a 3D array z(:,:,:)
!***********************************************************************************
	implicit none
	real(kind=8), intent(in)	:: xyz(:,:,:)
	real(kind=4), optional		:: pause
	character(len=*),optional	:: palette, terminal, filename, pm3d, contour, persist, input
	integer				:: nx, ny, nrow
	integer 			:: i, j, ierror, ios, file_unit
	character(len=100)		:: data_file_name, command_file_name, my_pause, my_persist
!***********************************************************************************
	nx=size(xyz(1,:,1))
	ny=size(xyz(1,1,:))
	nrow=nx*ny
!***********************************************************************************
	if (present(input)) then
		data_file_name='data_file_'//input//'.txt'
		command_file_name='command_file_'//input//'.txt'		
	else
		data_file_name='data_file.txt'
		command_file_name='command_file.txt'
	end if
!***********************************************************************************
	ierror=0	
	call get_unit(file_unit)	
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=data_file_name, status='replace', iostat=ios)	
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal data file.'
		stop
	end if
!***********************************************************************************	
! here we write the date to the data_file - the gnuplot will read this data later
!***********************************************************************************	
	do j=1,ny
		do i=1,nx
			write (file_unit,'(3E15.7)') xyz(1:3,i,j)
		end do
		write (file_unit,'(a)')
	end do
!***********************************************************************************	
	close (unit=file_unit)
!***********************************************************************************
	ierror = 0
	call get_unit(file_unit)
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=command_file_name, status='replace',	iostat=ios)
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal command file.'
		stop
	end if	
!***********************************************************************************
! here we write the commands to the commands file which gnuplot will execute
!***********************************************************************************
	my_persist='persist '
	if (present(persist).and.(persist=='no')) my_persist=' '
	if (present(terminal)) then
		write ( file_unit, '(a)' ) 'set terminal '// trim(output_terminal(terminal)) 
	if (present(filename)) then
		write ( file_unit, '(a)' ) 'set output "'// trim(filename) //'"' 
	else
		write ( file_unit, '(a)' ) 'set output "'//my_date_and_time()//'"' 
	end if
	else
		write ( file_unit, '(a)' ) 'set terminal ' // trim(default_terminal) // ' ' &
			& //trim(my_persist) // ' title  "Gnuplot"' 
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'set nokey'
	if (present(palette)) then
	if ((trim(palette).ne.'RGB').and.(trim(palette).ne.'HSV').and.(trim(palette).ne.'CMY').and.&
		& (trim(palette).ne.'YIQ').and.(trim(palette).ne.'XYZ')) then
		write ( file_unit, '(a)' ) 'set palette '// trim(palette)
	else
		write ( file_unit, '(a)' ) 'set palette model '// trim(palette)
	end if
	else
		write ( file_unit, '(a)' ) 'set palette model '// trim(default_palette)
	end if
!***********************************************************************************
	if (present(pm3d)) then
		write ( file_unit, '(a)' ) 'set '// pm3d	
	else
		write ( file_unit, '(a)' ) 'set surface'
		if (present(contour)) then
			if (contour=='surface') then 
				write ( file_unit, '(a)' ) 'set contour surface'
			elseif (contour=='both') then
				write ( file_unit, '(a)' ) 'set contour both'
			else 
				write ( file_unit, '(a)' ) 'set contour'
			end if
		end if
	end if
	write ( file_unit, '(a)' ) 'set hidden3d'
	write ( file_unit, '(a)' ) 'set parametric'
!***********************************************************************************	
	write ( file_unit, '(a)' ) 'splot "' // trim ( data_file_name ) // &
	& '" using 1:2:3 with lines palette'
!***********************************************************************************
	if (present(pause)) then
		if (pause<0.0) then
			write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
		else 
			write (	my_pause,'(e9.3)') pause
			write ( file_unit, '(a)' ) 'pause ' // trim(my_pause) 
		end if
	else
		write ( file_unit, '(a)' ) 'pause 0'
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'q'
	close ( unit = file_unit )
!***********************************************************************************
	call run_gnuplot (command_file_name)
!***********************************************************************************
	end subroutine surf_1
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine plot_4(x1,y1,x2,y2,x3,y3,x4,y4,style,pause,color1,color2,color3,color4,&
						& terminal,filename,polar,persist,input,linewidth)
!***********************************************************************************
! this subroutine plots 4 two-dimensional graphs in the same coordinate system
!***********************************************************************************
	implicit none
	real(kind=8), intent(in)	:: x1(:), y1(:), x2(:), y2(:), x3(:), y3(:), x4(:), y4(:)
	real(kind=4), optional		:: pause,linewidth
	character(len=*),optional	:: style, color1, color2, color3, color4, terminal, filename, polar,&
						& persist, input
	integer 			:: i, ierror, ios, file_unit, Nx1, Nx2, Nx3, Nx4, Nmax
	character(len=100)		:: data_file_name, command_file_name, my_linewidth
	integer, parameter		:: Nc=20
	character(len=Nc)		:: my_line_type1, my_line_type2, my_line_type3, my_line_type4, &
					& my_color1, my_color2, my_color3,  my_color4, my_range, my_pause, my_persist
!***********************************************************************************
	if (present(input)) then
		data_file_name='data_file_'//input//'.txt'
		command_file_name='command_file_'//input//'.txt'		
	else
		data_file_name='data_file.txt'
		command_file_name='command_file.txt'
	end if
!***********************************************************************************
	Nx1=size(x1)
	Nx2=size(x2)
	Nx3=size(x3)
	Nx4=size(x4)
	if ((size(x1).ne.size(y1)).or.(size(x2).ne.size(y2)).or.(size(x3).ne.size(y3)).or.(size(x4).ne.size(y4))) then
		print *,'subroutine plot ERROR: size(x) is not equal to size(y)'
		stop
	end if
	if (present(style).and.(len(style).ne.12)) then
		print *,'subroutine plot ERROR: argument "style" has wrong number of characters'
		stop
	end if
!***********************************************************************************
	ierror=0	
	call get_unit(file_unit)	
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_data - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=data_file_name, status='replace', iostat=ios)	
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal data file.'
		stop
	end if
!***********************************************************************************	
! here we write the date to the data_file - the gnuplot will read this data later
!***********************************************************************************
	Nmax=max(Nx1,Nx2,Nx3,Nx4)	
	do i=1,Nmax
		write (file_unit,'(8E15.7)') x1(min(i,Nx1)), y1(min(i,Nx1)), x2(min(i,Nx2)), y2(min(i,Nx2)), &
		& x3(min(i,Nx3)), y3(min(i,Nx3)), x4(min(i,Nx4)), y4(min(i,Nx4))
	end do
!***********************************************************************************	
	close (unit=file_unit)
!***********************************************************************************
	ierror = 0
	call get_unit(file_unit)
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_data - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=command_file_name, status='replace',	iostat=ios)
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal command file.'
		stop
	end if	
!***********************************************************************************
! here we write the commands to the commands file which gnuplot will execute
!***********************************************************************************
	my_line_type1='lines'
	if (present(style)) then
	if ((style(3:3)=='-')) then
		my_line_type1='linespoints'
	else
		my_line_type1='points'
	end if
	end if
	my_line_type2='lines'
	if (present(style)) then
	if ((style(6:6)=='-')) then
		my_line_type2='linespoints'
	else
		my_line_type2='points'
	end if
	end if
	my_line_type3='lines'
	if (present(style)) then
	if ((style(9:9)=='-')) then
		my_line_type3='linespoints'
	else
		my_line_type3='points'
	end if
	end if
	my_line_type4='lines'
	if (present(style)) then
	if ((style(12:12)=='-')) then
		my_line_type4='linespoints'
	else
		my_line_type4='points'
	end if
	end if
	if (present(linewidth)) then
		write (	my_linewidth,'(e9.3)') linewidth
	else
		my_linewidth=trim(default_linewidth)
	end if	
	if (present(color1)) then
		my_color1='"'//trim(color1)//'"'
	else
		my_color1='"'//trim(default_color1)//'"'
	end if
	if (present(color2)) then
		my_color2='"'//trim(color2)//'"'
	else
		my_color2='"'//trim(default_color2)//'"'
	end if
	if (present(color3)) then
		my_color3='"'//trim(color3)//'"'
	else
		my_color3='"'//trim(default_color3)//'"'
	end if
	if (present(color4)) then
		my_color4='"'//trim(color4)//'"'
	else
		my_color4='"'//trim(default_color4)//'"'
	end if
!***********************************************************************************
	my_persist='persist '
	if (present(persist).and.(persist=='no')) my_persist=' '
	if (present(terminal)) then
		write ( file_unit, '(a)' ) 'set terminal '// trim(output_terminal(terminal)) 
	if (present(filename)) then
		write ( file_unit, '(a)' ) 'set output "'// trim(filename) //'"' 
	else
		write ( file_unit, '(a)' ) 'set output "'//my_date_and_time()//'"' 
	end if
	else
		write ( file_unit, '(a)' ) 'set terminal ' // trim(default_terminal) // ' ' &
			& //trim(my_persist) // ' title  "Gnuplot"' 
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'unset key'
	if (present(polar).and.(polar=='yes')) then
		write (my_range,'(e15.7)') max(maxval(abs(y1)),maxval(abs(y2)),maxval(abs(y3)),maxval(abs(y4)))
		write ( file_unit, '(a)' ) 'set xrange [-'//trim(my_range)//':'//trim(my_range)//']'
		write ( file_unit, '(a)' ) 'set yrange [-'//trim(my_range)//':'//trim(my_range)//']'
		write ( file_unit, '(a)' ) 'set size square'
		write ( file_unit, '(a)' ) 'set polar'
		write ( file_unit, '(a)' ) 'set grid polar'
	else
		write ( file_unit, '(a)' ) 'set grid'
	end if
!***********************************************************************************
	if (present(style)) then
		write ( file_unit, '(a,i2,a)' ) 'plot "' // trim (data_file_name) &
		&//'" using 1:2 with ' // trim(my_line_type1) // ' pointtype ' // &
		& style(1:2) // ' linecolor rgb ' // trim(my_color1) // ' linewidth '// trim(my_linewidth) // ',\' 
		write ( file_unit, '(a,i2,a)' ) '     "'// trim (data_file_name) &
		&//'" using 3:4 with ' // trim(my_line_type2) // ' pointtype ' &
		&// style(4:5) // ' linecolor rgb ' // trim(my_color2) // ' linewidth '// trim(my_linewidth) //',\' 
		write ( file_unit, '(a,i2,a)' ) '     "'// trim (data_file_name) &
		&//'" using 5:6 with ' // trim(my_line_type3) // ' pointtype ' &
		&// style(7:8) // ' linecolor rgb ' // trim(my_color3) // ' linewidth '// trim(my_linewidth) // ',\' 
		write ( file_unit, '(a,i2,a)' ) '     "'// trim (data_file_name) &
		&//'" using 7:8 with ' // trim(my_line_type4) // ' pointtype ' &
		&// style(10:11) // ' linecolor rgb '// trim(my_color4)// ' linewidth '// trim(my_linewidth) 
	else 
		write ( file_unit, '(a,i2,a)' ) 'plot "' // trim (data_file_name) &
		& //'" using 1:2 with ' // trim(my_line_type1)  // ' linecolor rgb '&
		& // trim(my_color1) // ' linewidth '// trim(my_linewidth) // ',\' 
		write ( file_unit, '(a,i2,a)' ) '     "'// trim (data_file_name) &
		& //'" using 3:4 with ' // trim(my_line_type2)  // ' linecolor rgb '&
		& // trim(my_color2) // ' linewidth '// trim(my_linewidth) // ',\' 
		write ( file_unit, '(a,i2,a)' ) '     "'// trim (data_file_name) &
		& //'" using 5:6 with ' // trim(my_line_type3)  // ' linecolor rgb '&
		& // trim(my_color3) // ' linewidth '// trim(my_linewidth) // ',\' 
		write ( file_unit, '(a,i2,a)' ) '     "'// trim (data_file_name) &
		& //'" using 7:8 with ' // trim(my_line_type4)  // ' linecolor rgb '&
		& // trim(my_color4) // ' linewidth '// trim(my_linewidth) 
	end if
!***********************************************************************************
	if (present(pause)) then
		if (pause<0.0) then
			write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
		else 
			write (	my_pause,'(e9.3)') pause
			write ( file_unit, '(a)' ) 'pause ' // trim(my_pause) 
		end if
	else
		write ( file_unit, '(a)' ) 'pause 0'
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'q'
	close ( unit = file_unit )
!***********************************************************************************
	call run_gnuplot (command_file_name)
!***********************************************************************************
	end subroutine plot_4
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine plot_3(x1,y1,x2,y2,x3,y3,style,pause,color1,color2,color3,terminal,filename,polar,persist,input,linewidth)
!***********************************************************************************
! this subroutine plots 3 two-dimensional graphs in the same coordinate system
!***********************************************************************************
	implicit none
	real(kind=8), intent(in)	:: x1(:), y1(:), x2(:), y2(:), x3(:), y3(:)
	real(kind=4), optional		:: pause,linewidth
	character(len=*),optional	:: style, color1, color2, color3, terminal, filename, polar, persist, input
	integer 			:: i, ierror, ios, file_unit, Nx1, Nx2, Nx3, Nmax
	character(len=100)		:: data_file_name, command_file_name, my_linewidth
	integer, parameter		:: Nc=20
	character(len=Nc)		:: my_line_type1, my_line_type2, my_line_type3, my_color1, my_color2,&
					& my_color3, my_range, my_pause, my_persist
!***********************************************************************************
	if (present(input)) then
		data_file_name='data_file_'//input//'.txt'
		command_file_name='command_file_'//input//'.txt'		
	else
		data_file_name='data_file.txt'
		command_file_name='command_file.txt'
	end if
!***********************************************************************************
	Nx1=size(x1)
	Nx2=size(x2)
	Nx3=size(x3)
	if ((size(x1).ne.size(y1)).or.(size(x2).ne.size(y2)).or.(size(x3).ne.size(y3))) then
		print *,'subroutine plot ERROR: size(x) is not equal to size(y)'
		stop
	end if
	if (present(style).and.(len(style).ne.9)) then
		print *,'subroutine plot ERROR: argument "style" has wrong number of characters'
		stop
	end if
!***********************************************************************************
	ierror=0	
	call get_unit(file_unit)	
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_data - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=data_file_name, status='replace', iostat=ios)	
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal data file.'
		stop
	end if
!***********************************************************************************	
! here we write the date to the data_file - the gnuplot will read this data later
!***********************************************************************************
	Nmax=max(Nx1,Nx2,Nx3)	
	do i=1,Nmax
		write (file_unit,'(6E15.7)') x1(min(i,Nx1)), y1(min(i,Nx1)), x2(min(i,Nx2)), y2(min(i,Nx2)), &
		& x3(min(i,Nx3)), y3(min(i,Nx3))
	end do
!***********************************************************************************	
	close (unit=file_unit)
!***********************************************************************************
	ierror = 0
	call get_unit(file_unit)
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_data - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=command_file_name, status='replace',	iostat=ios)
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal command file.'
		stop
	end if	
!***********************************************************************************
! here we write the commands to the commands file which gnuplot will execute
!***********************************************************************************
	my_line_type1='lines'
	if (present(style)) then
	if ((style(3:3)=='-')) then
		my_line_type1='linespoints'
	else
		my_line_type1='points'
	end if
	end if
	my_line_type2='lines'
	if (present(style)) then
	if ((style(6:6)=='-')) then
		my_line_type2='linespoints'
	else
		my_line_type2='points'
	end if
	end if
	my_line_type3='lines'
	if (present(style)) then
	if ((style(9:9)=='-')) then
		my_line_type3='linespoints'
	else
		my_line_type3='points'
	end if
	end if
	if (present(linewidth)) then
		write (	my_linewidth,'(e9.3)') linewidth
	else
		my_linewidth=trim(default_linewidth)
	end if	
	if (present(color1)) then
		my_color1='"'//trim(color1)//'"'
	else
		my_color1='"'//trim(default_color1)//'"'
	end if
	if (present(color2)) then
		my_color2='"'//trim(color2)//'"'
	else
		my_color2='"'//trim(default_color2)//'"'
	end if
	if (present(color3)) then
		my_color3='"'//trim(color3)//'"'
	else
		my_color3='"'//trim(default_color3)//'"'
	end if
!***********************************************************************************
	my_persist='persist '
	if (present(persist).and.(persist=='no')) my_persist=' '
	if (present(terminal)) then
		write ( file_unit, '(a)' ) 'set terminal '// trim(output_terminal(terminal)) 
	if (present(filename)) then
		write ( file_unit, '(a)' ) 'set output "'// trim(filename) //'"' 
	else
		write ( file_unit, '(a)' ) 'set output "'//my_date_and_time()//'"' 
	end if
	else
		write ( file_unit, '(a)' ) 'set terminal ' // trim(default_terminal) // ' ' & 
			&//trim(my_persist) // ' title  "Gnuplot"' 
	end if

!***********************************************************************************
	write ( file_unit, '(a)' ) 'unset key'
	if (present(polar).and.(polar=='yes')) then
		write (my_range,'(e15.7)') max(maxval(abs(y1)),maxval(abs(y2)),maxval(abs(y3)))
		write ( file_unit, '(a)' ) 'set xrange [-'//trim(my_range)//':'//trim(my_range)//']'
		write ( file_unit, '(a)' ) 'set yrange [-'//trim(my_range)//':'//trim(my_range)//']'
		write ( file_unit, '(a)' ) 'set size square'
		write ( file_unit, '(a)' ) 'set polar'
		write ( file_unit, '(a)' ) 'set grid polar'
	else
		write ( file_unit, '(a)' ) 'set grid'
	end if	
!***********************************************************************************
	if (present(style)) then
		write ( file_unit, '(a,i2,a)' ) 'plot "' // trim (data_file_name) &
		&//'" using 1:2 with ' // trim(my_line_type1) // ' pointtype ' // &
		& style(1:2) // ' linecolor rgb ' // trim(my_color1) // ' linewidth '// trim(my_linewidth) // ',\' 
		write ( file_unit, '(a,i2,a)' ) '     "'// trim (data_file_name) &
		&//'" using 3:4 with ' // trim(my_line_type2) // ' pointtype ' &
		&// style(4:5) // ' linecolor rgb ' // trim(my_color2) // ' linewidth '// trim(my_linewidth) //',\' 
		write ( file_unit, '(a,i2,a)' ) '     "'// trim (data_file_name) &
		&//'" using 5:6 with ' // trim(my_line_type3) // ' pointtype ' &
		&// style(7:8) // ' linecolor rgb ' // trim(my_color3) // ' linewidth '// trim(my_linewidth) 
	else 
		write ( file_unit, '(a,i2,a)' ) 'plot "' // trim (data_file_name) &
		& //'" using 1:2 with ' // trim(my_line_type1)  // ' linecolor rgb '&
		& // trim(my_color1) // ' linewidth '// trim(my_linewidth) // ',\' 
		write ( file_unit, '(a,i2,a)' ) '     "'// trim (data_file_name) &
		& //'" using 3:4 with ' // trim(my_line_type2)  // ' linecolor rgb '&
		& // trim(my_color2) // ' linewidth '// trim(my_linewidth) // ',\' 
		write ( file_unit, '(a,i2,a)' ) '     "'// trim (data_file_name) &
		& //'" using 5:6 with ' // trim(my_line_type3)  // ' linecolor rgb '&
		& // trim(my_color3) // ' linewidth '// trim(my_linewidth) 
	end if
!***********************************************************************************
	if (present(pause)) then
		if (pause<0.0) then
			write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
		else 
			write (	my_pause,'(e9.3)') pause
			write ( file_unit, '(a)' ) 'pause ' // trim(my_pause) 
		end if
	else
		write ( file_unit, '(a)' ) 'pause 0'
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'q'
	close ( unit = file_unit )
!***********************************************************************************
	call run_gnuplot (command_file_name)
!***********************************************************************************
	end subroutine plot_3
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine plot_2(x1,y1,x2,y2,style,pause,color1,color2,terminal,filename,polar,persist,input,linewidth)
!***********************************************************************************
! this subroutine plots 2 two-dimensional graphs in the same coordinate system
!***********************************************************************************
	implicit none
	real(kind=8), intent(in)	:: x1(:), y1(:), x2(:), y2(:)
	real(kind=4), optional		:: pause,linewidth
	character(len=*),optional	:: style, color1, color2, terminal, filename, polar, persist, input
	integer 			:: i, ierror, ios, file_unit, Nx1, Nx2, Nmax
	character(len=100)		:: data_file_name, command_file_name, my_linewidth
	integer, parameter		:: Nc=20
	character(len=Nc)		:: my_line_type1, my_line_type2, my_color1, my_color2, my_range, my_pause, my_persist
!***********************************************************************************
	if (present(input)) then
		data_file_name='data_file_'//input//'.txt'
		command_file_name='command_file_'//input//'.txt'		
	else
		data_file_name='data_file.txt'
		command_file_name='command_file.txt'
	end if
!***********************************************************************************
	Nx1=size(x1)
	Nx2=size(x2)
	if ((size(x1).ne.size(y1)).or.(size(x2).ne.size(y2))) then
		print *,'subroutine plot ERROR: size(x) is not equal to size(y)'
		stop
	end if
	if (present(style).and.(len(style).ne.6)) then
		print *,'subroutine plot ERROR: argument "style" has wrong number of characters'
		stop
	end if
!***********************************************************************************
	ierror=0	
	call get_unit(file_unit)	
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_data - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=data_file_name, status='replace', iostat=ios)	
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal data file.'
		stop
	end if
!***********************************************************************************	
! here we write the date to the data_file - the gnuplot will read this data later
!***********************************************************************************
	Nmax=max(Nx1,Nx2)	
	do i=1,Nmax
		write (file_unit,'(4E15.7)') x1(min(i,Nx1)), y1(min(i,Nx1)), x2(min(i,Nx2)), y2(min(i,Nx2))
	end do
!***********************************************************************************	
	close (unit=file_unit)
!***********************************************************************************
	ierror = 0
	call get_unit(file_unit)
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_data - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=command_file_name, status='replace',	iostat=ios)
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal command file.'
		stop
	end if	
!***********************************************************************************
! here we write the commands to the commands file which gnuplot will execute
!***********************************************************************************
	my_line_type1='lines'
	if (present(style)) then
	if ((style(3:3)=='-')) then
		my_line_type1='linespoints'
	else
		my_line_type1='points'
	end if
	end if
	my_line_type2='lines'
	if (present(style)) then
	if ((style(6:6)=='-')) then
		my_line_type2='linespoints'
	else
		my_line_type2='points'
	end if
	end if
	if (present(linewidth)) then
		write (	my_linewidth,'(e9.3)') linewidth
	else
		my_linewidth=trim(default_linewidth)
	end if	
	if (present(color1)) then
		my_color1='"'//trim(color1)//'"'
	else
		my_color1='"'//trim(default_color1)//'"'
	end if
	if (present(color2)) then
		my_color2='"'//trim(color2)//'"'
	else
		my_color2='"'//trim(default_color2)//'"'
	end if
!***********************************************************************************
	my_persist='persist '
	if (present(persist).and.(persist=='no')) my_persist=' '
	if (present(terminal)) then
		write ( file_unit, '(a)' ) 'set terminal '// trim(output_terminal(terminal)) 
	if (present(filename)) then
		write ( file_unit, '(a)' ) 'set output "'// trim(filename) //'"' 
	else
		write ( file_unit, '(a)' ) 'set output "'//my_date_and_time()//'"' 
	end if
	else
		write ( file_unit, '(a)' ) 'set terminal ' // trim(default_terminal) // ' ' &
			& //trim(my_persist) // ' title  "Gnuplot"' 
	end if

!***********************************************************************************
	write ( file_unit, '(a)' ) 'unset key'
	if (present(polar).and.(polar=='yes')) then
		write (my_range,'(e15.7)') max(maxval(abs(y1)),maxval(abs(y2)))
		write ( file_unit, '(a)' ) 'set xrange [-'//trim(my_range)//':'//trim(my_range)//']'
		write ( file_unit, '(a)' ) 'set yrange [-'//trim(my_range)//':'//trim(my_range)//']'
		write ( file_unit, '(a)' ) 'set size square'
		write ( file_unit, '(a)' ) 'set polar'
		write ( file_unit, '(a)' ) 'set grid polar'
	else
		write ( file_unit, '(a)' ) 'set grid'
	end if		
!***********************************************************************************
	if (present(style)) then
		write ( file_unit, '(a,i2,a)' ) 'plot "' // trim (data_file_name) &
		&//'" using 1:2 with ' // trim(my_line_type1) // ' pointtype ' // &
		& style(1:2) // ' linecolor rgb ' // trim(my_color1) // ' linewidth '// trim(my_linewidth) // ',\' 
		write ( file_unit, '(a,i2,a)' ) '     "'// trim (data_file_name) &
		&//'" using 3:4 with ' // trim(my_line_type2) // ' pointtype ' &
		&// style(4:5) // ' linecolor rgb ' // trim(my_color2) // ' linewidth '// trim(my_linewidth) 
	else 
		write ( file_unit, '(a,i2,a)' ) 'plot "' // trim (data_file_name) &
		& //'" using 1:2 with ' // trim(my_line_type1)  // ' linecolor rgb '&
		& // trim(my_color1) // ' linewidth '// trim(my_linewidth) // ',\' 
		write ( file_unit, '(a,i2,a)' ) '     "'// trim (data_file_name) &
		& //'" using 3:4 with ' // trim(my_line_type2)  // ' linecolor rgb '&
		& // trim(my_color2) // ' linewidth '// trim(my_linewidth) 
	end if
!***********************************************************************************
	if (present(pause)) then
		if (pause<0.0) then
			write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
		else 
			write (	my_pause,'(e9.3)') pause
			write ( file_unit, '(a)' ) 'pause ' // trim(my_pause) 
		end if
	else
		write ( file_unit, '(a)' ) 'pause 0'
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'q'
	close ( unit = file_unit )
!***********************************************************************************
	call run_gnuplot (command_file_name)
!***********************************************************************************
	end subroutine plot_2
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine plot_1(x1,y1,style,pause,color1,terminal,filename,polar,persist,input,linewidth)
!***********************************************************************************
! this subroutine plots a two-dimensional graph
!***********************************************************************************
	implicit none
	real(kind=8), intent(in)	:: x1(:), y1(:)
	real(kind=4), optional		:: pause,linewidth
	character(len=*),optional	:: style, color1, terminal, filename, polar, persist, input
	integer 			:: i, ierror, ios, file_unit, Nx1
	character(len=100)		:: data_file_name, command_file_name, my_linewidth
	integer, parameter		:: Nc=20
	character(len=Nc)		:: my_line_type1, my_color1, my_range, my_pause, my_persist
!***********************************************************************************
	if (present(input)) then
		data_file_name='data_file_'//input//'.txt'
		command_file_name='command_file_'//input//'.txt'		
	else
		data_file_name='data_file.txt'
		command_file_name='command_file.txt'
	end if
!***********************************************************************************
	Nx1=size(x1)
	if ((size(x1).ne.size(y1))) then
		print *,'subroutine plot ERROR: size(x) is not equal to size(y)'
		stop
	end if
	if (present(style).and.(len(style).ne.3)) then
		print *,'subroutine plot ERROR: argument "style" has wrong number of characters'
		stop
	end if
!***********************************************************************************
	ierror=0	
	call get_unit(file_unit)	
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_data - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=data_file_name, status='replace', iostat=ios)	
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal data file.'
		stop
	end if
!***********************************************************************************	
! here we write the date to the data_file - the gnuplot will read this data later
!***********************************************************************************	
	do i=1,Nx1
		write (file_unit,'(2E15.7)') x1(i), y1(i)
	end do
!***********************************************************************************	
	close (unit=file_unit)
!***********************************************************************************
	ierror = 0
	call get_unit(file_unit)
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_data - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=command_file_name, status='replace',	iostat=ios)
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal command file.'
		stop
	end if	
!***********************************************************************************
! here we write the commands to the commands file which gnuplot will execute
!***********************************************************************************
	my_line_type1='lines'
	if (present(style)) then
	if ((style(3:3)=='-')) then
		my_line_type1='linespoints'
	else
		my_line_type1='points'
	end if
	end if
	if (present(linewidth)) then
		write (	my_linewidth,'(e9.3)') linewidth
	else
		my_linewidth=trim(default_linewidth)
	end if	
	if (present(color1)) then
		my_color1='"'//trim(color1)//'"'
	else
		my_color1='"'//trim(default_color1)//'"'
	end if
!***********************************************************************************
	my_persist='persist '
	if (present(persist).and.(persist=='no')) my_persist=' '
	if (present(terminal)) then
		write ( file_unit, '(a)' ) 'set terminal '// trim(output_terminal(terminal)) 
	if (present(filename)) then
		write ( file_unit, '(a)' ) 'set output "'// trim(filename) //'"' 
	else
		write ( file_unit, '(a)' ) 'set output "'//my_date_and_time()//'"' 
	end if
	else
		write ( file_unit, '(a)' ) 'set terminal ' // trim(default_terminal) // ' ' &
			& //trim(my_persist) //' title  "Gnuplot"' 
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'unset key'
	if (present(polar).and.(polar=='yes')) then
		write (my_range,'(e15.7)') maxval(abs(y1))
		write ( file_unit, '(a)' ) 'set xrange [-'//trim(my_range)//':'//trim(my_range)//']'
		write ( file_unit, '(a)' ) 'set yrange [-'//trim(my_range)//':'//trim(my_range)//']'
		write ( file_unit, '(a)' ) 'set size square'
		write ( file_unit, '(a)' ) 'set polar'
		write ( file_unit, '(a)' ) 'set grid polar'
	else
		write ( file_unit, '(a)' ) 'set grid'
	end if	
!***********************************************************************************
	if (present(style)) then
		write ( file_unit, '(a,i2,a)' ) 'plot "' // trim (data_file_name) &
		&//'" using 1:2 with ' // trim(my_line_type1) // ' pointtype ' // &
		& style(1:2) // ' linecolor rgb ' // trim(my_color1) // ' linewidth '// trim(my_linewidth) 
	else 
		write ( file_unit, '(a,i2,a)' ) 'plot "' // trim (data_file_name) &
		& //'" using 1:2 with ' // trim(my_line_type1)  // ' linecolor rgb '&
		& // trim(my_color1) // ' linewidth '// trim(my_linewidth)  
	end if
!***********************************************************************************
	if (present(pause)) then
		if (pause<0.0) then
			write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
		else 
			write (	my_pause,'(e9.3)') pause
			write ( file_unit, '(a)' ) 'pause ' // trim(my_pause) 
		end if
	else
		write ( file_unit, '(a)' ) 'pause 0'
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'q'
	close ( unit = file_unit )
!***********************************************************************************
	call run_gnuplot (command_file_name) 
!***********************************************************************************
	end subroutine plot_1
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine run_gnuplot(command_file_name)
!***********************************************************************************
	implicit none
	character (len = 100) command
	character (len = *) command_file_name
	integer status
	integer system
!***********************************************************************************
!  Issue a command to the system that will startup GNUPLOT, using
!  the file we just wrote as input.
!***********************************************************************************
	write (command, *) 'gnuplot ' // trim (command_file_name)		
	status=system(trim(command))	
	if (status.ne.0) then
		print *,'RUN_GNUPLOT - Fatal error!'
		stop
	end if	
	return
!***********************************************************************************
	end subroutine run_gnuplot
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine get_unit(iunit)
!***********************************************************************************
	implicit none
	integer i
	integer ios
	integer iunit
	logical lopen
!***********************************************************************************	
	iunit=0
	do i=1,99
		if (i/= 5 .and. i/=6) then	
			inquire (unit=i, opened=lopen, iostat=ios)
			if (ios==0) then
				if (.not.lopen) then
					iunit=i
					return
				end if
			end if
		
		end if
	end do	
	return
	end subroutine get_unit
!***********************************************************************************
	end module gnufor2
module utils
  public :: sort
  
  interface sort
    module procedure sort1
    module procedure sort2
  end interface sort
contains
! from: http://www.cs.mtu.edu/~shene/courses/cs201/notes/chap08/sorting.f90
! --------------------------------------------------------------------
! integer function  findminimum():
!    this function returns the location of the minimum in the section
! between start and end.
! --------------------------------------------------------------------

   integer function  findminimum(x, start, end)
      implicit  none
      integer, dimension(1:), intent(in) :: x
      integer, intent(in)                :: start, end
      integer                            :: minimum
      integer                            :: location
      integer                            :: i

      minimum  = x(start)               ! assume the first is the min
      location = start                  ! record its position
      do i = start+1, end               ! start with next elements
         if (x(i) < minimum) then       !   if x(i) less than the min?
            minimum  = x(i)             !      yes, a new minimum found
            location = i                !      record its position
         end if
      end do
      findminimum = location            ! return the position
   end function  findminimum

! --------------------------------------------------------------------
! subroutine  swap():
!    this subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

   subroutine  swap1(a, b)
      implicit  none
      integer, intent(inout) :: a, b
      integer                :: temp

      temp = a
      a    = b
      b    = temp
   end subroutine  swap1

! --------------------------------------------------------------------
! subroutine  sort():
!    this subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

   subroutine  sort1(x, size)
      implicit  none
      integer, dimension(1:), intent(inout) :: x
      integer, intent(in)                   :: size
      integer                               :: i
      integer                               :: location

      do i = 1, size-1                  ! except for the last
         location = findminimum(x, i, size)     ! find min from this to last
         call  swap1(x(i), x(location))  ! swap this and the minimum
      end do
   end subroutine  sort1
   
! --------------------------------------------------------------------
! subroutine  swap(): [Key, Value]'
!    this subroutine swaps the [key, value]' of its two formal arguments.
! --------------------------------------------------------------------

   subroutine  swap2(a, b)
      implicit  none
      integer, dimension(2), intent(inout) :: a, b
      integer, dimension(2)                :: temp

      temp = a
      a    = b
      b    = temp
   end subroutine  swap2

! --------------------------------------------------------------------
! subroutine  sort(): [Key, Value]'
!    this subroutine receives a 2xN array x() and sorts it into ascending
! order based on the key stored in the first row
! This is currently doing insertion sort, need to do quick/merge sort
! --------------------------------------------------------------------

   subroutine  sort2(x, size)
      implicit  none
      integer, dimension(2,size), intent(inout) :: x
      integer, intent(in)                       :: size
      integer                                   :: i
      integer                                   :: location

      do i = 1, size-1                           ! except for the last
         location = findminimum(x(1,:), i, size) ! find min key from this to last
         call  swap2(x(:,i), x(:,location))      ! swap this and the minimum
      end do
   end subroutine  sort2
end module utils
module mgmres
contains
subroutine atx_cr ( n, nz_num, ia, ja, a, x, w )

!*****************************************************************************80
!
!! ATX_CR computes A'*x for a matrix stored in sparse compressed row form.
!
!  Discussion:
!
!    The Sparse Compressed Row storage format is used.
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA[I] through IA[I+1]-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A'.
!
!    Output, real ( kind = 8 ) W(N), the value of A'*X.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  w(1:n) = 0.0D+00

  do i = 1, n
    k1 = ia(i)
    k2 = ia(i+1) - 1
    w(ja(k1:k2)) = w(ja(k1:k2)) + a(k1:k2) * x(i)
  end do

  return
end
subroutine atx_st ( n, nz_num, ia, ja, a, x, w )

!*****************************************************************************80
!
!! ATX_ST computes A'*x for a matrix stored in sparset triplet form.
!
!  Discussion:
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  For instance, the K-th nonzero
!    entry in the matrix is stored by:
!
!      A(K) = value of entry,
!      IA(K) = row of entry,
!      JA(K) = column of entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(NZ_NUM), JA(NZ_NUM), the row and column
!    indices of the matrix values.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A'.
!
!    Output, real ( kind = 8 ) W(N), the value of A'*X.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(nz_num)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  w(1:n) = 0.0D+00

  do k = 1, nz_num
    i = ia(k)
    j = ja(k)
    w(j) = w(j) + a(k) * x(i)
  end do

  return
end
subroutine ax_cr ( n, nz_num, ia, ja, a, x, w )

!*****************************************************************************80
!
!! AX_CR computes A*x for a matrix stored in sparse compressed row form.
!
!  Discussion:
!
!    The Sparse Compressed Row storage format is used.
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA[I] through IA[I+1]-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) W(N), the value of A*X.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  w(1:n) = 0.0D+00

  do i = 1, n
    k1 = ia(i)
    k2 = ia(i+1) - 1
    w(i) = w(i) + dot_product ( a(k1:k2), x(ja(k1:k2)) )
  end do

  return
end
subroutine ax_st ( n, nz_num, ia, ja, a, x, w )

!*****************************************************************************80
!
!! AX_ST computes A*x for a matrix stored in sparset triplet form.
!
!  Discussion:
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  For instance, the K-th nonzero
!    entry in the matrix is stored by:
!
!      A(K) = value of entry,
!      IA(K) = row of entry,
!      JA(K) = column of entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(NZ_NUM), JA(NZ_NUM), the row and column
!    indices of the matrix values.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) W(N), the value of A*X.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(nz_num)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  w(1:n) = 0.0D+00

  do k = 1, nz_num
    i = ia(k)
    j = ja(k)
    w(i) = w(i) + a(k) * x(j)
  end do

  return
end
subroutine diagonal_pointer_cr ( n, nz_num, ia, ja, ua )

!*****************************************************************************80
!
!! DIAGONAL_POINTER_CR finds diagonal entries in a sparse compressed row matrix.
!
!  Discussion:
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA[I] through IA[I+1]-1.
!
!    The array UA can be used to locate the diagonal elements of the matrix.
!
!    It is assumed that every row of the matrix includes a diagonal element,
!    and that the elements of each row have been ascending sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!    On output, the order of the entries of JA may have changed because of
!    the sorting.
!
!    Output, integer ( kind = 4 ) UA(N), the index of the diagonal element
!    of each row.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) ua(n)

  ua(1:n) = -1

  do i = 1, n
    do k = ia(i), ia(i+1) - 1
      if ( ja(k) == i ) then
        ua(i) = k
      end if
    end do
  end do

  return
end
subroutine ilu_cr ( n, nz_num, ia, ja, a, ua, l )

!*****************************************************************************80
!
!! ILU_CR computes the incomplete LU factorization of a matrix.
!
!  Discussion:
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input, integer ( kind = 4 ) UA(N), the index of the diagonal element
!    of each row.
!
!    Output, real ( kind = 8 ) L(NZ_NUM), the ILU factorization of A.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) iw(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jrow
  integer ( kind = 4 ) jw
  integer ( kind = 4 ) k
  real ( kind = 8 ) l(nz_num)
  real ( kind = 8 ) tl
  integer ( kind = 4 ) ua(n)
!
!  Copy A.
!
  l(1:nz_num) = a(1:nz_num)

  do i = 1, n
!
!  IW points to the nonzero entries in row I.
!
    iw(1:n) = -1

    do k = ia(i), ia(i+1) - 1
      iw(ja(k)) = k
    end do

    do j = ia(i), ia(i+1) - 1                  ! ith row entries
      jrow = ja(j)                             ! column corresponding to each entry
      if ( i <= jrow ) then
        exit
      end if
      tl = l(j) * l(ua(jrow))
      l(j) = tl
      do jj = ua(jrow) + 1, ia(jrow+1) - 1
        jw = iw(ja(jj))
        if ( jw /= -1 ) then
          l(jw) = l(jw) - tl * l(jj)
        end if
      end do
    end do

    ua(i) = j

    if ( jrow /= i ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ILU_CR - Fatal error!'
      write ( *, '(a)' ) '  JROW ~= I'
      write ( *, '(a,i8)' ) '  JROW = ', jrow
      write ( *, '(a,i8)' ) '  I    = ', i
      stop
    end if

    if ( l(j) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ILU_CR - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step I = ', i
      write ( *, '(a,i8,a)' ) '  L(', j, ') = 0.0'
      stop
    end if

    l(j) = 1.0D+00 / l(j)

  end do

  l(ua(1:n)) = 1.0D+00 / l(ua(1:n))

  return
end
subroutine lus_cr ( n, nz_num, ia, ja, l, ua, r, z )

!*****************************************************************************80
!
!! LUS_CR applies the incomplete LU preconditioner.
!
!  Discussion:
!
!    The linear system M * Z = R is solved for Z.  M is the incomplete
!    LU preconditioner matrix, and R is a vector supplied by the user.
!    So essentially, we're solving L * U * Z = R.
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = 8 ) L(NZ_NUM), the matrix values.
!
!    Input, integer ( kind = 4 ) UA(N), the index of the diagonal element
!    of each row.
!
!    Input, real ( kind = 8 ) R(N), the right hand side.
!
!    Output, real ( kind = 8 ) Z(N), the solution of the system M * Z = R.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(nz_num)
  real ( kind = 8 ) l(nz_num)
  real ( kind = 8 ) r(n)
  integer ( kind = 4 ) ua(n)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) z(n)
!
!  Copy R in.
!
  w(1:n) = r(1:n)
!
!  Solve L * w = w where L is unit lower triangular.
!
  do i = 2, n
    do j = ia(i), ua(i) - 1
      w(i) = w(i) - l(j) * w(ja(j))
    end do
  end do
!
!  Solve U * w = w, where U is upper triangular.
!
  do i = n, 1, -1
    do j = ua(i) + 1, ia(i+1) - 1
      w(i) = w(i) - l(j) * w(ja(j))
    end do
    w(i) = w(i) / l(ua(i))
  end do
!
!  Copy Z out.
!
  z(1:n) = w(1:n)

  return
end
subroutine mgmres_st ( n, nz_num, ia, ja, a, x, rhs, itr_max, mr, tol_abs, &
  tol_rel, verbose )

!*****************************************************************************80
!
!! MGMRES_ST applies restarted GMRES to a sparse triplet matrix.
!
!  Discussion:
!
!    The linear system A*X=B is solved iteratively.
!
!    The matrix A is assumed to be stored in sparse triplet form.  Only
!    the nonzero entries of A are stored.  For instance, the K-th nonzero
!    entry in the matrix is stored by:
!
!      A(K) = value of entry,
!      IA(K) = row of entry,
!      JA(K) = column of entry.
!
!    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
!    corrections to the code on 31 May 2007.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero matrix values.
!
!    Input, integer ( kind = 4 ) IA(NZ_NUM), JA(NZ_NUM), the row and column
!    indices of the matrix values.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input/output, real ( kind = 8 ) X(N); on input, an approximation to
!    the solution.  On output, an improved approximation.
!
!    Input, real ( kind = 8 ) RHS(N), the right hand side of the linear system.
!
!    Input, integer ( kind = 4 ) ITR_MAX, the maximum number of (outer)
!    iterations to take.
!
!    Input, integer ( kind = 4 ) MR, the maximum number of (inner) iterations
!    to take.  0 < MR <= N.
!
!    Input, real ( kind = 8 ) TOL_ABS, an absolute tolerance applied to the
!    current residual.
!
!    Input, real ( kind = 8 ) TOL_REL, a relative tolerance comparing the
!    current residual to the initial residual.
!
  implicit none

  integer ( kind = 4 ) mr
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) av
  real ( kind = 8 ) c(1:mr)
  real ( kind = 8 ), parameter :: delta = 1.0D-03
  real ( kind = 8 ) g(1:mr+1)
  real ( kind = 8 ) h(1:mr+1,1:mr)
  real ( kind = 8 ) htmp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(nz_num)
  integer ( kind = 4 ) itr
  integer ( kind = 4 ) itr_max
  integer ( kind = 4 ) itr_used
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_copy
  real ( kind = 8 ) mu
  real ( kind = 8 ) r(1:n)
  real ( kind = 8 ) rho
  real ( kind = 8 ) rho_tol
  real ( kind = 8 ) rhs(1:n)
  real ( kind = 8 ) s(1:mr)
  real ( kind = 8 ) tol_abs
  real ( kind = 8 ) tol_rel
  real ( kind = 8 ) v(1:n,1:mr+1)
  logical :: verbose 
  real ( kind = 8 ) x(1:n)
  real ( kind = 8 ) y(1:mr+1)

  itr_used = 0

  if ( n < mr ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MGMRES_ST - Fatal error!'
    write ( *, '(a)' ) '  N < MR.'
    write ( *, '(a,i8)' ) '  N = ', n
    write ( *, '(a,i8)' ) '  MR = ', mr
    stop
  end if

  do itr = 1, itr_max

    call ax_st ( n, nz_num, ia, ja, a, x, r )

    r(1:n) = rhs(1:n) - r(1:n)

    rho = sqrt ( dot_product ( r(1:n), r(1:n) ) )

    if ( verbose ) then
      write ( *, '(a,i8,a,g14.6)' ) '  ITR = ', itr, '  Residual = ', rho
    end if

    if ( itr == 1 ) then
      rho_tol = rho * tol_rel
    end if

    v(1:n,1) = r(1:n) / rho

    g(1) = rho
    g(2:mr+1) = 0.0D+00

    h(1:mr+1,1:mr) = 0.0D+00

    do k = 1, mr

      k_copy = k

      call ax_st ( n, nz_num, ia, ja, a, v(1:n,k), v(1:n,k+1) )

      av = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      do j = 1, k
        h(j,k) = dot_product ( v(1:n,k+1), v(1:n,j) )
        v(1:n,k+1) = v(1:n,k+1) - h(j,k) * v(1:n,j)
      end do

      h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      if ( av + delta * h(k+1,k) == av ) then

        do j = 1, k
          htmp = dot_product ( v(1:n,k+1), v(1:n,j) )
          h(j,k) = h(j,k) + htmp
          v(1:n,k+1) = v(1:n,k+1) - htmp * v(1:n,j)
        end do

        h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      end if

      if ( h(k+1,k) /= 0.0D+00 ) then
        v(1:n,k+1) = v(1:n,k+1) / h(k+1,k)
      end if

      if ( 1 < k ) then

        y(1:k+1) = h(1:k+1,k)

        do j = 1, k - 1
          call mult_givens ( c(j), s(j), j, y(1:k+1) )
        end do

        h(1:k+1,k) = y(1:k+1)

      end if

      mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )
      c(k) = h(k,k) / mu
      s(k) = -h(k+1,k) / mu
      h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
      h(k+1,k) = 0.0D+00
      call mult_givens ( c(k), s(k), k, g(1:k+1) )
      rho = abs ( g(k+1) )

      itr_used = itr_used + 1

      if ( verbose ) then
        write ( *, '(a,i8,a,g14.6)' ) '  K =   ', k, '  Residual = ', rho
      end if

      if ( rho <= rho_tol .and. rho <= tol_abs ) then
        exit
      end if

    end do

    k = k_copy - 1

    y(k+1) = g(k+1) / h(k+1,k+1)

    do i = k, 1, -1
      y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) ) / h(i,i)
    end do

    do i = 1, n
      x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
    end do

    if ( rho <= rho_tol .and. rho <= tol_abs ) then
      exit
    end if

  end do

  if ( verbose ) then
    write ( *, '(a)'       ) ' '
    write ( *, '(a)'       ) 'MGMRES_ST:'
    write ( *, '(a,i8)'    ) '  Iterations = ', itr_used
    write ( *, '(a,g14.6)' ) '  Final residual = ', rho
    write ( *, '(a)'       ) ' '
  end if

  return
end
subroutine mult_givens ( c, s, k, g )

!*****************************************************************************80
!
!! MULT_GIVENS applies a Givens rotation to two successive entries of a vector.
!
!  Discussion:
!
!    In order to make it easier to compare this code with the Original C,
!    the vector indexing is 0-based.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) C, S, the cosine and sine of a Givens
!    rotation.
!
!    Input, integer ( kind = 4 ) K, indicates the location of the first
!    vector entry.
!
!    Input/output, real ( kind = 8 ) G(1:K+1), the vector to be modified.
!    On output, the Givens rotation has been applied to entries G(K) and G(K+1).
!
  implicit none

  integer ( kind = 4 ) k

  real ( kind = 8 ) c
  real ( kind = 8 ) g(1:k+1)
  real ( kind = 8 ) g1
  real ( kind = 8 ) g2
  real ( kind = 8 ) s

  g1 = c * g(k) - s * g(k+1)
  g2 = s * g(k) + c * g(k+1)

  g(k)   = g1
  g(k+1) = g2

  return
end
subroutine pmgmres_ilu_cr ( n, nz_num, ia, ja, a, x, rhs, itr_max, mr, &
  tol_abs, tol_rel, verbose )

!*****************************************************************************80
!
!! PMGMRES_ILU_CR applies the preconditioned restarted GMRES algorithm.
!
!  Discussion:
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!    This routine uses the incomplete LU decomposition for the
!    preconditioning.  This preconditioner requires that the sparse
!    matrix data structure supplies a storage position for each diagonal
!    element of the matrix A, and that each diagonal element of the
!    matrix A is not zero.
!
!    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
!    corrections to the code on 31 May 2007.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 2012
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero matrix values.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column indices
!    of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input/output, real ( kind = 8 ) X(N); on input, an approximation to
!    the solution.  On output, an improved approximation.
!
!    Input, real ( kind = 8 ) RHS(N), the right hand side of the linear system.
!
!    Input, integer ( kind = 4 ) ITR_MAX, the maximum number of (outer)
!    iterations to take.
!
!    Input, integer ( kind = 4 ) MR, the maximum number of (inner) iterations
!    to take.  MR must be less than N.
!
!    Input, real ( kind = 8 ) TOL_ABS, an absolute tolerance applied to the
!    current residual.
!
!    Input, real ( kind = 8 ) TOL_REL, a relative tolerance comparing the
!    current residual to the initial residual.
!
  implicit none

  integer ( kind = 4 ) mr
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) av
  real ( kind = 8 ) c(mr+1)
  real ( kind = 8 ), parameter :: delta = 1.0D-03
  real ( kind = 8 ) g(mr+1)
  real ( kind = 8 ) h(mr+1,mr)
  real ( kind = 8 ) htmp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) itr
  integer ( kind = 4 ) itr_max
  integer ( kind = 4 ) itr_used
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_copy
  real ( kind = 8 ) l(ia(n+1)+1)
  real ( kind = 8 ) mu
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) rho
  real ( kind = 8 ) rho_tol
  real ( kind = 8 ) rhs(n)
  real ( kind = 8 ) s(mr+1)
  real ( kind = 8 ) tol_abs
  real ( kind = 8 ) tol_rel
  integer ( kind = 4 ) ua(n)
  real ( kind = 8 ) v(n,mr+1);
  logical :: verbose
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(mr+1)

  itr_used = 0

  call rearrange_cr ( n, nz_num, ia, ja, a )

  call diagonal_pointer_cr ( n, nz_num, ia, ja, ua )

  call ilu_cr ( n, nz_num, ia, ja, a, ua, l )

  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMGMRES_ILU_CR'
    write ( *, '(a,i4)' ) '  Number of unknowns = ', n
  end if

  do itr = 1, itr_max

    call ax_cr ( n, nz_num, ia, ja, a, x, r )

    r(1:n) = rhs(1:n) - r(1:n)

    call lus_cr ( n, nz_num, ia, ja, l, ua, r, r )

    rho = sqrt ( dot_product ( r, r ) )

    if ( verbose ) then
      write ( *, '(a,i4,a,g14.6)' ) '  ITR = ', itr, '  Residual = ', rho
    end if

    if ( itr == 1 ) then
      rho_tol = rho * tol_rel
    end if

    v(1:n,1) = r(1:n) / rho

    g(1) = rho
    g(2:mr+1) = 0.0D+00

    h(1:mr+1,1:mr) = 0.0D+00

    do k = 1, mr

      k_copy = k

      call ax_cr ( n, nz_num, ia, ja, a, v(1:n,k), v(1:n,k+1) )

      call lus_cr ( n, nz_num, ia, ja, l, ua, v(1:n,k+1), v(1:n,k+1) )

      av = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      do j = 1, k
        h(j,k) = dot_product ( v(1:n,k+1), v(1:n,j) )
        v(1:n,k+1) = v(1:n,k+1) - v(1:n,j) * h(j,k)
      end do

      h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      if ( ( av + delta * h(k+1,k)) == av ) then
        do j = 1, k
          htmp = dot_product ( v(1:n,k+1), v(1:n,j) )
          h(j,k) = h(j,k) + htmp
          v(1:n,k+1) = v(1:n,k+1) - htmp * v(1:n,j)
        end do
        h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )
      end if

      if ( h(k+1,k) /= 0.0D+00 ) then
        v(1:n,k+1) = v(1:n,k+1) / h(k+1,k)
      end if

      if ( 1 < k ) then
        y(1:k+1) = h(1:k+1,k)
        do j = 1, k - 1
          call mult_givens ( c(j), s(j), j, y )
        end do
        h(1:k+1,k) = y(1:k+1)
      end if

      mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )

      c(k) = h(k,k) / mu
      s(k) = -h(k+1,k) / mu
      h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
      h(k+1,k) = 0.0D+00
      call mult_givens ( c(k), s(k), k, g )

      rho = abs ( g(k+1) )

      itr_used = itr_used + 1

      if ( verbose ) then
        write ( *, '(a,i4,a,g14.6)' ) '  K = ', k, '  Residual = ', rho
      end if

      if ( rho <= rho_tol .and. rho <= tol_abs ) then
        exit
      end if

    end do

    k = k_copy - 1

    y(k+1) = g(k+1) / h(k+1,k+1)

    do i = k, 1, -1
      y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) ) / h(i,i)
    end do

    do i = 1, n
      x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
    end do

    if ( rho <= rho_tol .and. rho <= tol_abs ) then
      exit
    end if

  end do

  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMGMRES_ILU_CR:'
    write ( *, '(a,i6)' ) '  Iterations = ', itr_used
    write ( *, '(a,g14.6)' ) '  Final residual = ', rho
  end if

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of real ( kind = 8 ) values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine rearrange_cr ( n, nz_num, ia, ja, a )

!*****************************************************************************80
!
!! REARRANGE_CR sorts a sparse compressed row matrix.
!
!  Discussion:
!
!    This routine guarantees that the entries in the CR matrix
!    are properly sorted.
!
!    After the sorting, the entries of the matrix are rearranged in such
!    a way that the entries of each column are listed in ascending order
!    of their column values.
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), the compressed row indices.
!
!    Input/output, integer ( kind = 4 ) JA(NZ_NUM), the column indices.
!    On output, these may have been rearranged by the sorting.
!
!    Input/output, real ( kind = 8 ) A(NZ_NUM), the matrix values.  On output,
!    the matrix values may have been moved somewhat because of the sorting.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) i4temp
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) r8temp

  do i = 1, n

    do k = ia(i), ia(i+1) - 2
      do l = k + 1, ia(i+1) - 1

        if ( ja(l) < ja(k) ) then
          i4temp = ja(l)
          ja(l)  = ja(k)
          ja(k)  = i4temp

          r8temp = a(l)
          a(l)   = a(k)
          a(k)   = r8temp
        end if

      end do
    end do

  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
end module mgmres
