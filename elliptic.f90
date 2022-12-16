program elliptic
	implicit none

	character (len=1) :: isReal
	integer :: p
	integer :: a, b
	integer :: root, i,j,k,l
	integer :: findRoot
	integer, dimension (:), allocatable :: points_x, points_y

	print *, "Enter the field"
	read *, p

	allocate(points_x(p))
	allocate(points_y(p))

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
	do i=0, p-1
		do j=0, p-1
			if (MOD(i*i,p) == MOD(j*j*j+a*j + b,p)) then
				points_x(i) = j
				points_y(i) = i
				print *,"Found point: ", i, j
			else
				points_x(i) = -1
				points_y(i) = -1
			end if
		end do
	end do

	print *, "Found point: inf"
end program elliptic

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