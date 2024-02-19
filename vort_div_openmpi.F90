MODULE WIND_DECOMP
REAL, ALLOCATABLE, DIMENSION(:,:) :: U_CHI, U_PSI, V_PSI, V_CHI
CONTAINS

SUBROUTINE VORT_DIV_WIND(DX, DY, LEN_VORT_X, LEN_VORT_Y, VORT_MAG, DIVERGENCE, LIMIT, start_index, end_index)
! Built by: Michael P. Vossen
! Last Updated: 2/15/2024

! Subroutine to calculate rotational wind and divergent wind
! Comes from Oertel and Schemm (2021)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INPUT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! LEN_VORT_X (INTEGER): size of the x demension of the vorticity and divergence array (Not needed in Python Argument)
! LEN_VORT_Y (INTEGER): size of the y demension of the vorticity and divergence array (Not needed in Python Argument)
! DX (REAL): The grid size in the x demension.  Must be in meters
! DY (REAL): The grid size in the y demension.  Must be in meters 
! VORT_MAG (REAL ARRAY): Vorticity Array
! DIVERGENCE (REAL ARRAY): Diveregence Array
! LIMIT (INTEGER): the size of one side of a box (meters) that limits the area to use for the calculation


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! OUTPUT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!start_index (INTEGER): The starting i index that the current processor uses
!end_index (INTEGER): The ending i index that the current processor uses


USE mpi
IMPLICIT NONE

!f2py intent(in) :: LEN_VORT_X, LEN_VORT_Y, LIMIT
INTEGER, INTENT(IN) :: LEN_VORT_X, LEN_VORT_Y, LIMIT 

!f2py intent(in) :: DX, DY
REAL, INTENT(IN) :: DX, DY

!f2py intent(in) :: VORT_MAG, DIVERGENCE
REAL, INTENT(IN), DIMENSION(1:LEN_VORT_Y, 1:LEN_VORT_X) :: VORT_MAG, DIVERGENCE


!f2py intent(out) :: start_index, end_index
INTEGER, INTENT(OUT) :: start_index, end_index

REAL :: UPSI_SUM, VPSI_SUM, UCHI_SUM, VCHI_SUM, R, PI, DIFF_X, DIIF_Y, R2, DIFF_Y

INTEGER :: I, J, X1, Y1, &
PASS, START_X, START_Y, END_X, &
END_Y, LIM_X, LIM_Y, HALF_X, &
HALF_Y, ierr, &
size, rank, iterations_per_process, namelength,&
grid_remain, i_size, actual_i

character(len=15) :: processorname


!for the limit find how many grid points the distance is
LIM_X = LIMIT/DX
LIM_Y = LIMIT/DY

PI=4.D0*DATAN(1.D0)

!for the limit find how many grid point in each direction we need to go
HALF_X = LIM_X/2
HALF_Y = LIM_Y/2

! Initialize MPI
call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
call MPI_GET_PROCESSOR_NAME(processorName, namelength, ierr)

! Calculate the number of iterations each process will handle
! Assuming LEN_VORT_X and LEN_VORT_Y are divisible by size


iterations_per_process = LEN_VORT_X / size

! Calculate the starting and ending indices for this process
start_index = (rank * iterations_per_process) + 1
end_index = start_index + iterations_per_process - 1

grid_remain = LEN_VORT_X - end_index

if (grid_remain .LT. iterations_per_process) THEN
end_index = end_index+grid_remain
endif

i_size = (end_index - start_index) + 1

allocate(U_PSI(LEN_VORT_Y, i_size))
allocate(V_PSI(LEN_VORT_Y, i_size))
allocate(U_CHI(LEN_VORT_Y, i_size))
allocate(V_CHI(LEN_VORT_Y, i_size))

print *, "RANK =", RANK, " I SIZE=", i_size, "START =", start_index, " END =", end_index

DO I = start_index, end_index
    DO J = 1,LEN_VORT_Y
        UPSI_SUM = 0
        VPSI_SUM = 0
        UCHI_SUM = 0
        VCHI_SUM = 0

        START_X = I - HALF_X
        END_X = I + HALF_X
        START_Y = J - HALF_Y
        END_Y = J + HALF_Y

        !the following is for the limit if it falls outside the domain.
        IF (START_X .LT. 1) THEN
            START_X = 1
        ENDIF
        IF (START_Y .LT. 1) THEN
            START_Y = 1
        ENDIF
        IF (END_X .GT. LEN_VORT_X) THEN
            END_X = LEN_VORT_X
        ENDIF
        IF (END_Y .GT. LEN_VORT_Y) THEN
            END_Y = LEN_VORT_Y
        ENDIF

        
        DO X1 = START_X,END_X
            DO Y1 = START_Y,END_Y
                
                IF (I .eq. X1 .AND. J .eq. Y1) THEN
                    PASS=1
                ELSE
                   DIFF_X = I - X1
                   DIFF_Y = J - Y1
                   R = SQRT((DIFF_X * DIFF_X) + (DIFF_Y * DIFF_Y)) * DX
                   R2 = R * R
                   UPSI_SUM = UPSI_SUM + (VORT_MAG(Y1, X1) * (( -1 * DY * DIFF_Y) / R2) * DX * DY)
                   VPSI_SUM = VPSI_SUM + (VORT_MAG(Y1, X1) * ((DIFF_X *DX)/R2) * DX * DY)
                   UCHI_SUM = UCHI_SUM + (DIVERGENCE(Y1, X1) * ((DIFF_X * DX) / R2) * DX * DY)
                   VCHI_SUM = VCHI_SUM + (DIVERGENCE(Y1, X1) * ((DIFF_Y * DY) / R2) * DX * DY)

                ENDIF
            ENDDO
        ENDDO
	actual_i = (I-start_index) + 1
        U_PSI(J,actual_i) = (1/(2*PI)) * UPSI_SUM
        V_PSI(J,actual_i) = (1/(2*PI)) * VPSI_SUM
        U_CHI(J,actual_i) = (1/(2*PI)) * UCHI_SUM
        V_CHI(J,actual_i) = (1/(2*PI)) * VCHI_SUM


    ENDDO
ENDDO

END SUBROUTINE

END MODULE WIND_DECOMP



