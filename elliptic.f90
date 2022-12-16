program elliptic
	implicit none

	character (len=1) :: isReal
	integer :: p
	integer :: a, b
	integer :: numPoints
	integer :: root, i,j,k,l
	integer :: findRoot
	integer, dimension (:), allocatable :: points_x, points_y, points_z

	print *, "Enter the field (must be a prime not equal to 2 or 3)"
	read *, p

	allocate(points_x(p+1))
	allocate(points_y(p+1))
	allocate(points_z(p+1))

	if (p==2 .or. p==3) then
		print *, "Error: Field cannot be over 2 or 3"
		stop
	end if

	print *, "Field is: ", p

	print *, "Enter the coefficients of the Elliptic curve: y^2 = x^3 + ax + b"
	print *, "First enter a:"
	read *,a
	print *, "Now enter b:"
	read *,b

	if (a < 0 .or. b < 0 .or. a >= p .or. b >= p) then
		print *, "Error: Ensure that the coefficients lie between 0 and p-1"
		stop
	end if

	! Generate points
	points_x(0) = 0
	points_y(0) = 1
	points_z(0) = 0
	print *, "Found point: (0,1,0)"
	numPoints = 1
	do i=0, p-1
		do j=0, p-1
			if (MOD(i*i,p) == MOD(j*j*j+a*j + b,p)) then
				points_x(numPoints) = j
				points_y(numPoints) = i
				points_z(numPoints) = 1
				print *,"Found point: ", j, i
				numPoints = numPoints+1
			end if
		end do
	end do

	print *, "Total number of points: ", numPoints
end program elliptic

! Computes n*P where P=(x,y,z) and stores the result in (out_x,out_z,out_z)
! z-coordinate indicates whether the point is infinity -- following convention from sage
subroutine pointMultiplication(x,y,z,n,a,b,p,out_x,out_y,out_z)
	implicit none

	integer :: x,y,z,a,b,p,out_x,out_y,out_z
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
		out_x = (s**2) - x1-x2
		out_y = s*(x1-out_x)-y1
		out_z = 1
	else if (x1==x2 .and. y1==y2) then  ! addition when points are equal
		s = (3*(x1**2)+a)/(2*y1)
		out_x = (s**2)-2*x1
		out_y = s*(x1-out_x)-y1
		out_z = 1
	else   ! This should cover all cases, if this else statement is executed, then there is a bug
		print *, "Error: Code went down a branch that should not be possible"
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