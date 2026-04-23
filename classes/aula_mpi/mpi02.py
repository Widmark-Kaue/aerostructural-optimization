from mpi4py import MPI

comm = MPI.COMM_WORLD # grupo de procesos 


myid = comm.Get_rank() # id do processo
numprocs = comm.Get_size() # número de processos

if myid != 3: 
    print('Hello do proc %d de %d procs!' % (myid, numprocs))
else:
    print('Não quero dar oi!')