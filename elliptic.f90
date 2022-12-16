program elliptic
	implicit none

	character (len=1) :: char
	integer :: p
	integer :: a, b
	integer :: numPoints
	integer :: root, i,j,k,l, input, x1,y1,z1,x2,y2,z2,out_x,out_y,out_z
	integer :: findRoot
	integer, dimension (:), allocatable :: points_x, points_y, points_z

	print '("Enter the field (must be a prime not equal to 2 or 3)")',
	read *, p

	allocate(points_x(p+1))
	allocate(points_y(p+1))
	allocate(points_z(p+1))

	if (p==2 .or. p==3) then
		print '("Error: Field cannot be over 2 or 3")',
		stop
	end if

	print '("Field is: ", I0)', p

	print '("Enter the coefficients of the Elliptic curve: y^2 = x^3 + ax + b")',
	print '("Enter a and b:")',
	read *,a,b
	!print '("Now enter b:")'
	!read *,b

	if (a < 0 .or. b < 0 .or. a >= p .or. b >= p) then
		print '("Error: Ensure that the coefficients lie between 0 and p-1")',
		stop
	end if

	! Generate points
	points_x(0) = 0
	points_y(0) = 1
	points_z(0) = 0
	print '("List of points:")',
	numPoints = 1
	print '(I0,": (",I0,",",I0,",",I0,")")', numPoints, 0, 1, 0
	do i=0, p-1
		do j=0, p-1
			if (MOD(i*i,p) == MOD(j*j*j+a*j + b,p)) then
				points_x(numPoints) = j
				points_y(numPoints) = i
				points_z(numPoints) = 1
				numPoints = numPoints+1
				print '(I0,": (",I0,",",I0,",",I0,")")', numPoints, j, i, 1
			end if
		end do
	end do

	print '("Total number of points: ",I0)', numPoints

	do while (.TRUE.)
		print '("")',
		print '("Enter either 1, 2 or 3:")',
		print '("1: Add two points on the elliptic curve")',
		print '("2: Multiply a point by an integer")',
		print '("3: Close")',

		read *, input

		if (input==1) then
			print '("Is the first point at infinity? (y/n)")',
			read *, char
			if(char=='y') then
				x1=0
				y1=1
				z1=0
			else if (char=='n') then
				print '("Enter the coordinates of the first point.")',
				read *, x1,y1
				z1 = 1
			else
				print '("Error: Please enter either y or n.")',
				exit
			end if

			print '("Is the second point at infinity? (y/n)")',
			read *, char
			if(char=='y') then
				x2=0
				y2=1
				z2=0
			else if (char=='n') then
				print '("Enter the coordinates of the second point.")',
				read *, x2,y2
				z2 = 1
			else
				print '("Error: Please enter either y or n.")'
				exit
			end if

			call addPoints(x1,y1,z1,x2,y2,z2,a,b,p,out_x,out_y,out_z)
			print '("The sum of (", I0, ",", I0, ",", I0, ") and (", I0,",",I0,",",I0, ") is (",I0,",",I0,",",I0,").")' &
				, x1, y1, z1, x2,y2,z2, out_x,out_y,out_z
			

		else if (input==2) then
		else if (input == 3) then
			print '("Closing.")',
			stop
		else
			print '("Command not recognized. Input either 1, 2 or 3.")',
		end if
	end do

end program elliptic

! Computes n*P where P=(x,y,z), and stores the result in (out_x,out_z,out_z)
! z-coordinate indicates whether the point is infinity -- following convention from sage
subroutine pointMultiplication(x,y,z,n,a,b,p,out_x,out_y,out_z)
	implicit none

	integer :: x,y,z,n,a,b,p,out_x,out_y,out_z
end subroutine pointMultiplication

! Adds two points P1=(x1,y1,z1) and P2=(x2,y2,z2), and stores the result in (out_x,out_y,out_z)
! z-coordinates indicate whether the point is infinity -- following convention from sage
subroutine addPoints(x1,y1,z1,x2,y2,z2,a,b,p,out_x,out_y,out_z)
	implicit none

	integer :: x1,y1,z1,x2,y2,z2,a,b,p,out_x,out_y,out_z
	integer :: s

	if (z1==0) then  ! if P1 is point at infinity then return P2
		out_x = x2
		out_y = y2
		out_z = z2
	else if (z2==0) then ! if P2 is point at infinity then return P1
		out_x = x1
		out_y = y1
		out_z = z1
	! P1,P2 are not infinity from this point onwards
	else if (x1==x2 .and. y1==-y2) then  ! if slope is vertical then return point at infinity
		out_x = 0   
		out_y = 1
		out_z = 0
	else if (x1 .ne. x2) then  ! addition when points are not equal
		s = (y2-y1)/(x2-x1)

		out_x = MOD((s**2) - x1-x2,p)
		out_y = MOD(s*(x1-out_x)-y1,p)
		out_z = 1
	else if (x1==x2 .and. y1==y2) then  ! addition when points are equal
		s = (3*(x1**2)+a)/(2*y1)

		out_x = MOD((s**2)-2*x1,p)
		out_y = MOD(s*(x1-out_x)-y1,p)
		out_z = 1

	else   ! This should cover all cases, if this else statement is executed, then there is a bug
		print *, "Error: Code went down a branch that should not be possible"
	end if

	! Make numbers positive if the modulos made them negative
	if (out_x < 0) then
		out_x = out_x + p
	end if
	if (out_y < 0) then
		out_y = out_y + p
	end if
end subroutine addPoints

! Checks if x has a square root modulo p
function findRoot (x,p)
	implicit none

	integer :: findRoot
	integer :: p,x,i

	findRoot = -1
	do i=0,p-1
		if (MOD(i*i,p)==MOD(x,p)) then
			findRoot = i
		end if
	end do
end function findRoot