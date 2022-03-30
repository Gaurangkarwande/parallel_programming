/*MPI  Test program*/

#include "mpi.h"
#include  "stdio.h"

int main(int argc, char *argv[])
{
   int rank, numtasks, dest, source, count, tag=1;
   char inmsg, outmsg='x';
   MPI_Status Stat;
   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
   if (rank==0)
   {
      dest=1; source=1;
      for(dest=1, source=1; dest < numtasks && source<numtasks; source++, dest++)
      {
         MPI_Send(&outmsg, 1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
         MPI_Recv(&inmsg, 1, MPI_CHAR, source, tag, MPI_COMM_WORLD, &Stat);
      }
   }
   else
   {
      dest=0; source=0;
      MPI_Recv(&inmsg, 1, MPI_CHAR, source, tag, MPI_COMM_WORLD, &Stat);
      MPI_Send(&outmsg, 1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
   }
   MPI_Get_count(&Stat, MPI_CHAR, &count);
   printf("Task %d: Received %d char(s) from task %d with tag %d \n", 
            rank, count, Stat.MPI_SOURCE, Stat.MPI_TAG);
   MPI_Finalize();
   return 0;
}