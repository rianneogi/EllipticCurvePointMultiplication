program elliptic
	implicit none

	character (len=1) :: char
	real :: a, b
	real :: x1,y1,z1,x2,y2,z2,out_x,out_y,out_z, root
	integer :: i,j,k,l, input, n
	integer :: checkPoint

	print '("Enter the coefficients of the Elliptic curve: y^2 = x^3 + ax + b")',
	print '("Enter a and b:")',
	read *,a,b
	!print '("Enter a:")'
	!read *,a
	!print '("Now enter b:")'
	!read *,b

	do while (.TRUE.)
		print '("")',
		print '("Enter either 1, 2 or 3:")',
		print '("1: Add two points on the elliptic curve")',
		print '("2: Multiply a point by an integer")',
		print '("3: Close")',

		read *, input

		if (input==1) then
			print '("Enter the x,y,z-coordinates of the first point. &
				The z-coordinate denotes whether the point is at infinity -- following the convention from sage.")',
			read *,x1,y1,z1
			
			if (checkPoint(x1,y1,z1,a,b) == 0) then
				print '("Error: This is not a point on the curve.")'
				exit
			end if

			print '("Enter the x,y,z-coordinates of the second point. &
				The z-coordinate denotes whether the point is at infinity -- following the convention from sage.")',
			read *,x2,y2,z2
			
			if (checkPoint(x2,y2,z2,a,b) == 0) then
				print '("Error: This is not a point on the curve.")'
				exit
			end if

			call addPoints(x1,y1,z1,x2,y2,z2,a,b,out_x,out_y,out_z)
			print '("The sum of (", F0.2, ",", F0.2, ",", F0.2, ") and (", F0.2,",",F0.2,",",F0.2, ") is (",F0.2,",",F0.2,",",F0.2,").")' &
				, x1, y1, z1, x2,y2,z2, out_x,out_y,out_z
			

		else if (input==2) then
			print '("Enter the x,y,z-coordinates of the point. &
				The z-coordinate denotes whether the point is at infinity -- following the convention from sage.")',
			read *,x1,y1,z1
			
			if (checkPoint(x1,y1,z1,a,b) == 0) then
				print '("Error: This is not a point on the curve.")'
				exit
			end if

			print '("Enter the integer n.")'
			read *, n

			call pointMultiplication(x1,y1,z1,n,a,b,out_x,out_y,out_z)
			print '("(", F0.2, ",", F0.2, ",", F0.2, ") times ", I0, " is (",F0.2,",",F0.2,",",F0.2,").")' &
				,x1,y1,z1,n,out_x,out_y,out_z

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
subroutine pointMultiplication(x,y,z,n,a,b,out_x,out_y,out_z)
	implicit none

	real :: x,y,z,a,b,out_x,out_y,out_z, t_x,t_y,t_z
	integer :: n
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
			call addPoints(out_x,out_y,out_z,out_x,out_y,out_z,a,b,t_x,t_y,t_z)
			out_x = t_x
			out_y = t_y
			out_z = t_z

			m = m/2
		else  ! if the i-th bit is 1, add current point to itself and then add P to it
			call addPoints(out_x,out_y,out_z,out_x,out_y,out_z,a,b,t_x,t_y,t_z)
			out_x = t_x
			out_y = t_y
			out_z = t_z

			call addPoints(out_x,out_y,out_z,x,y,z,a,b,t_x,t_y,t_z)
			out_x = t_x
			out_y = t_y
			out_z = t_z

			m = (m-1)/2
		end if
	end do

	if (n < 0 .and. out_z==1) then  ! if n is negative and output is not the point at infinity, then flip the y-coordinate
		out_y = -out_y
	end if
end subroutine pointMultiplication

! Adds two points P1=(x1,y1,z1) and P2=(x2,y2,z2), and stores the result in (out_x,out_y,out_z)
! z-coordinates indicate whether the point is infinity -- following convention from sage
subroutine addPoints(x1,y1,z1,x2,y2,z2,a,b,out_x,out_y,out_z)
	implicit none

	real :: x1,y1,z1,x2,y2,z2,a,b,out_x,out_y,out_z
	real :: s

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

! Checks if the point x,y is a point on the curve
function checkPoint(x,y,z,a,b)
	implicit none

	real :: x,y,z,a,b
	integer :: checkPoint

	if (z==0) then  ! if z is 0 then this must be the point at infinity
		if (x==0 .and. y==1) then
			checkPoint = 1
		else
			checkPoint = 0
		end if
	else if (z==1) then ! otherwise check if RHS is equal to LHS of the Elliptic curve
		if (y*y == x*x*x + a*x + b) then
			checkPoint = 1
		else 
			checkPoint = 0
		end if
	else ! z!=0 and z!=1
		checkPoint = 0
	end if
end function checkPoint

real function log2(x)
  implicit none
  real, intent(in) :: x

  log2 = log(x) / log(2.)
end function