CFLAGS = -std=c++11 -Wall -O3 -DTIMER


test_mpi: CC = mpic++
test_mpi: main_mpi.o LanczosInverseSparse.o LanczosSparse.o CG_mpi.o qr.o Sch_Matrix.o util.o Schrodinger.o
	$(CC) $(CFLAGS) main_mpi.o LanczosInverseSparse.o LanczosSparse.o CG_mpi.o qr.o Sch_Matrix.o util.o Schrodinger.o -o test_mpi

main_mpi.o: main_mpi.cpp
	$(CC) $(CFLAGS) -c main_mpi.cpp

LanczosInverseSparse.o: LanczosInverseSparse.c
	$(CC) $(CFLAGS) -c LanczosInverseSparse.c

LanczosSparse.o: LanczosSparse.c
	$(CC) $(CFLAGS) -c LanczosSparse.c

CG_mpi.o: CG_mpi.cpp
	$(CC) $(CFLAGS) -c CG_mpi.cpp

qr.o: qr.c
	$(CC) $(CFLAGS) -c qr.c

Sch_Matrix.o: Sch_Matrix.c
	$(CC) $(CFLAGS) -c Sch_Matrix.c

util.o: util.c
	$(CC) $(CFLAGS) -c util.c	

Schrodinger.o: Schrodinger.cpp
	$(CC) $(CFLAGS) -c Schrodinger.cpp

clean:
	rm -rf *o test_mpi