from mpi4py import MPI

comm = MPI.COMM_WORLD # grupo de procesos 


myid = comm.Get_rank() # id do processo
numprocs = comm.Get_size() # número de processos


if myid == 0:
    data = list(range(numprocs)) # cria uma lista de dados para cada processo
    print('Antes Proc 0 data:', data)
else:
    data = None

data = comm.scatter(data, root=0) # scatter: distribui os dados do processo raiz para todos os processos do grupo

print('Depois Proc:', myid, 'data:', data)
