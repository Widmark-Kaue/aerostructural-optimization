from mpi4py import MPI

comm = MPI.COMM_WORLD # grupo de procesos 


myid = comm.Get_rank() # id do processo
numprocs = comm.Get_size() # número de processos

data = None

if myid == 0 : 
    data = 1234

data = comm.bcast(data, root=0) # broadcast: envia a mensagem para todos os processos do grupo
                                # bloqueador, ou seja, o processo fica travado até receber a mensagem

print('Proc:', myid, 'data:', data)
