#include <iostream>
#include <cstdlib>
#include <chrono>
#include <thread>
#include "mpi.h"

using namespace std;
void CounterProcess()
{
	MPI_Request request;
	MPI_Status status;
	int ProcNum, ProcRank;
	int recvMsg = 0;
	const int stopMsg = -1;
	int counter = 0;
	int result = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Bcast(&result, 1, MPI_INT, 0, MPI_COMM_WORLD);
	srand((unsigned)time(NULL) + ProcRank);

	while (true)
	{
		this_thread::sleep_for(chrono::milliseconds(rand() % 1000 + 10));
		MPI_Isend(&counter, 1, MPI_INT, ProcRank % (ProcNum-1) +1, 0, MPI_COMM_WORLD, &request);
		if (counter == -1) break;
		MPI_Recv(&recvMsg, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		cout << "Process " << ProcRank << " received " << recvMsg << endl;
		counter++;
		if (rand() % 10 == 0 || recvMsg == -1 || recvMsg > 9)
			counter = -1;
		result *= 2;
	}
	int ProcSum = 0;
	MPI_Reduce(&result, &ProcSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
}

int CollectorProccess()
{
	int val0 = 2;
	MPI_Bcast(&val0, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int ProcSum = 0, TotalSum = 0;
	MPI_Reduce(&ProcSum, &TotalSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	return TotalSum;
}




int main(int argc, char* argv[])
{
	int ProcNum, ProcRank, RecvRank;
	int result = -1;
	MPI_Status Status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	if (ProcRank == 0)
	{
		result = CollectorProccess();
		cout << endl << "Result: " << result << endl;
	}
	else
		CounterProcess();

	MPI_Finalize();

	return 0;
}


