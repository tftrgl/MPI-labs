#include <iostream>
#include <cstdlib>
#include <chrono>
#include <thread>
#include "mpi.h"

using namespace std;
void CounterProcess()
{
	MPI_Status Status;
	int ProcNum, ProcRank;
	int recvMsg = 0;
	const int stopMsg = -1;
	int counter = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	srand((unsigned)time(NULL) + ProcRank);
	while (true)
	{
		this_thread::sleep_for(chrono::milliseconds(rand() % 1000 + 10));
		MPI_Send(&counter, 1, MPI_INT, (ProcRank + 1) % ProcNum, 0, MPI_COMM_WORLD);
		if (counter == -1) break;
		MPI_Recv(&recvMsg, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &Status);
		cout << "Process " << ProcRank << " received " << recvMsg << endl;
		counter++;
		if (rand() % 10 == 0 || recvMsg == -1 || recvMsg > 9)
			counter = -1;
	}
}




int main(int argc, char* argv[])
{
	int ProcNum, ProcRank, RecvRank;
	MPI_Status Status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	CounterProcess();
	MPI_Finalize();
	return 0;
}




