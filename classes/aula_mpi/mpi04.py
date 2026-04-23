from mpi4py import MPI

comm = MPI.COMM_WORLD # grupo de procesos 


myid = comm.Get_rank() # id do processo
numprocs = comm.Get_size() # número de processos

data = None

if myid == 0 : 
    data = 1234
    comm.send(data, dest=1, tag=42)

elif myid == 1:
    data = comm.recv(source=0, tag=42) # trava a execução do processo até receber a mensagem

print('Proc:', myid, 'data:', data)