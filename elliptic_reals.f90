program elliptic
	implicit none

	character (len=1) :: isReal
	integer :: p
	real :: a, b

	print *,"Is the Field is over reals? (y/n)"

	read *, isReal

	if (isReal=="y") then
		print *, "The Field is over the reals"

		print *, "Enter the coefficients of the Elliptic curve: y^2 = x^3 + ax + b"
		print *, "First enter a:"
		read *,a
		print *, "Now enter b:"
		read *,b

	else if (isReal=="n") then
		print *, "Enter the field"
		read *, p

		if (p==2 .or. p==3) then
			print *, "Field cannot be over 2 or 3"
			stop
		end if

		print *, "Field is: ", p

		print *, "Enter the coefficients of the Elliptic curve: y^2 = x^3 + ax + b"
		print *, "First enter a:"
		read *,a
		print *, "Now enter b:"
		read *,b
	else
		print *, "Error: Please input either y or n"
		stop 
	end if
end program elliptic

!function overReals