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
end program elliptic