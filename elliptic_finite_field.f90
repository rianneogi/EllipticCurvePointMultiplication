program elliptic
	implicit none

	character (len=1) :: char
	integer :: p
	integer :: a, b
	integer :: numPoints
	integer :: root, i,j,k,l, input, x1,y1,z1,x2,y2,z2,out_x,out_y,out_z,n
	integer :: findRoot, checkPoint, inv
	integer, dimension (:), allocatable :: points_x, points_y, points_z

	print '("Enter the field (must be a prime not equal to 2 or 3)")',
	read *, p

	allocate(points_x(2*p+1))
	allocate(points_y(2*p+1))
	allocate(points_z(2*p+1))

	if (p==2 .or. p==3) then
		print '("Error: Field cannot be over 2 or 3")',
		stop
	end if

	print '("Field is: ", I0)', p

	print '("Enter the coefficients of the Elliptic curve: y^2 = x^3 + ax + b")',
	print '("Enter a and b:")',
	read *,a,b
	!print '("Enter a:")'
	!read *,a
	!print '("Now enter b:")'
	!read *,b

	if (a < 0 .or. b < 0 .or. a >= p .or. b >= p) then
		print '("Error: Ensure that the coefficients lie between 0 and p-1")',
		stop
	end if

	! Generate points
	points_x(1) = 0
	points_y(1) = 1
	points_z(1) = 0
	print '("List of points:")',
	numPoints = 1
	print '(I0,": (",I0,",",I0,",",I0,")")', numPoints, 0, 1, 0
	do i=0, p-1
		do j=0, p-1
			if (MOD(i*i,p) == MOD(j*j*j+a*j + b,p)) then
				numPoints = numPoints+1
				points_x(numPoints) = j
				points_y(numPoints) = i
				points_z(numPoints) = 1
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
		print '("3: Compute the inverse of an element modulo p")',
		print '("4: Close")',

		read *, input

		if (input==1) then
			print '("Enter the x,y,z-coordinates of the first point. &
				The z-coordinate denotes whether the point is at infinity -- following the convention from sage.")',
			read *,x1,y1,z1
			
			if (checkPoint(x1,y1,z1,a,b,p) == 0) then
				print '("Error: This is not a point on the curve.")'
				exit
			end if

			print '("Enter the x,y,z-coordinates of the second point. &
				The z-coordinate denotes whether the point is at infinity -- following the convention from sage.")',
			read *,x2,y2,z2
			
			if (checkPoint(x2,y2,z2,a,b,p) == 0) then
				print '("Error: This is not a point on the curve.")'
				exit
			end if

			call addPoints(x1,y1,z1,x2,y2,z2,a,b,p,out_x,out_y,out_z)
			print '("The sum of (", I0, ",", I0, ",", I0, ") and (", I0,",",I0,",",I0, ") is (",I0,",",I0,",",I0,").")' &
				, x1, y1, z1, x2,y2,z2, out_x,out_y,out_z
			

		else if (input==2) then
			print '("Enter the x,y,z-coordinates of the point. &
				The z-coordinate denotes whether the point is at infinity -- following the convention from sage.")',
			read *,x1,y1,z1
			
			if (checkPoint(x1,y1,z1,a,b,p) == 0) then
				print '("Error: This is not a point on the curve.")'
				exit
			end if

			print '("Enter the integer n.")'
			read *, n

			call pointMultiplication(x1,y1,z1,n,a,b,p,out_x,out_y,out_z)
			print '("(", I0, ",", I0, ",", I0, ") times ", I0, " is (",I0,",",I0,",",I0,").")' &
				,x1,y1,z1,n,out_x,out_y,out_z

		else if (input == 3) then
			print '("Enter the element x.")'
			read *, n
			print '("The inverse of ", I0, " modulo ", I0, " is ", I0)', n,p,inv(n,p)

		else if (input == 4) then
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

	integer :: x,y,z,n,a,b,p,out_x,out_y,out_z, t_x,t_y,t_z
	integer :: m, len, i
	real :: log2
	integer, dimension (:), allocatable :: binary

	if (n < 0) then
		m = -n
	else
		m = n
	end if

	out_x = x
	out_y = y
	out_z = z
	t_x=x
	t_y=y
	t_z=z

	! Generate the binary representation of m
	allocate(binary(2*INT(log2(REAL(m)))+1))
	len=1
	do while (m >= 2)
		if (MOD(m,2)==0) then
			binary(len)=0
			m = m/2
		else
			binary(len)=1
			m = (m-1)/2
		end if
		len=len+1
	end do

	! Traverse the generated binary representation in reverse
	do i=1,len-1
		if (binary(len-i)==0) then  ! if the i-th bit is 0, add current point to itself
			call addPoints(out_x,out_y,out_z,out_x,out_y,out_z,a,b,p,t_x,t_y,t_z)
			out_x = t_x
			out_y = t_y
			out_z = t_z

			m = m/2
			print '("Doubled.")'
		else  ! if the i-th bit is 1, add current point to itself and then add P to it
			call addPoints(out_x,out_y,out_z,out_x,out_y,out_z,a,b,p,t_x,t_y,t_z)
			out_x = t_x
			out_y = t_y
			out_z = t_z

			call addPoints(out_x,out_y,out_z,x,y,z,a,b,p,t_x,t_y,t_z)
			out_x = t_x
			out_y = t_y
			out_z = t_z

			m = (m-1)/2
			print '("Doubled and added.")'
		end if
	end do

	if ((n < 0) .and. (out_z==1) .and. (out_y .ne. 0)) then  ! if n is negative and output is not the point at infinity, then flip the y-coordinate
		out_y = p-out_y
	end if
end subroutine pointMultiplication

! Adds two points P1=(x1,y1,z1) and P2=(x2,y2,z2), and stores the result in (out_x,out_y,out_z)
! z-coordinates indicate whether the point is infinity -- following convention from sage
subroutine addPoints(x1,y1,z1,x2,y2,z2,a,b,p,out_x,out_y,out_z)
	implicit none

	integer :: x1,y1,z1,x2,y2,z2,a,b,p,out_x,out_y,out_z
	integer :: s
	integer :: inv

	if (z1==0) then  ! if P1 is point at infinity then return P2
		out_x = x2
		out_y = y2
		out_z = z2
	else if (z2==0) then ! if P2 is point at infinity then return P1
		out_x = x1
		out_y = y1
		out_z = z1
	! P1,P2 are not infinity from this point onwards
	else if ((x1==x2) .and. ((y1==p-y2) .or. (y1==y2 .and. y1==0))) then  ! if slope is vertical then return point at infinity (change this for reals)
		out_x = 0   
		out_y = 1
		out_z = 0
	else if (x1 .ne. x2) then  ! addition when points are not equal
		s = MOD((y2-y1)*inv(x2-x1,p),p)

		out_x = MOD((s**2) - x1-x2,p) 
		out_y = MOD(s*(x1-out_x)-y1,p)
		out_z = 1
	else if (x1==x2 .and. y1==y2) then  ! addition when points are equal
		s = MOD((3*(x1**2)+a)*inv(2*y1,p),p)  ! Need to take inverse modulo p here

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

! Computes the inverse of x modulo p by running extended GCD
function inv(x,p)
	implicit none

	integer :: x,p
	integer :: inv
	integer :: r,s,t,old_r,old_s,old_t,q,tmp_r,tmp_t,tmp_s

	old_r = MOD(x,p)
	if (old_r < 0) then
		old_r = old_r+p
	end if
	r = p
	old_s = 1
	s = 0
	old_t = 0
	t = 1

	do while (r .ne. 0)
		q = old_r/r

		tmp_r = r
		r = old_r - q*r
		old_r = tmp_r

		tmp_s = s
		s = old_s - q*s
		old_s = tmp_s

		tmp_t = t
		t = old_t - q*t
		old_t = tmp_t
	end do

	inv = MOD(old_s,p)
	if (inv < 0) then
		inv = inv+p
	end if
end function inv

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

! Checks if the point x,y is a point on the curve
function checkPoint(x,y,z,a,b,p)
	implicit none

	integer :: x,y,z,a,b,p
	integer :: checkPoint

	if (z==0) then  ! if z is 0 then this must be the point at infinity
		if (x==0 .and. y==1) then
			checkPoint = 1
		else
			checkPoint = 0
		end if
	else if (z==1) then ! otherwise check if RHS is equal to LHS of the Elliptic curve
		if (MOD(y*y,p) == MOD(x*x*x + a*x + b,p)) then
			checkPoint = 1
		else 
			checkPoint = 0
		end if
	else ! z!=0 and z!=1
		checkPoint = 0
	end if
end function checkPoint

! Computes log of x base 2
real function log2(x)
  implicit none
  real, intent(in) :: x

  log2 = log(x) / log(2.)
end function