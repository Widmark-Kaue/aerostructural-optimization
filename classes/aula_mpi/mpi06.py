from mpi4py import MPI

comm = MPI.COMM_WORLD # grupo de procesos 


myid = comm.Get_rank() # id do processo
numprocs = comm.Get_size() # número de processos

data = myid

print('Antes Proc:', myid, 'data:', data)

data = comm.reduce(data, MPI.SUM, root=0) # reduce: reduz os dados de todos os processos do grupo para um único valor, usando a operação especificada (neste caso, soma)

print('Depois Proc:', myid, 'data:', data)
