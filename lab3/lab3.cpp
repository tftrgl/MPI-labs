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
	int result = 1;
	int factor = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Scatter(nullptr, 1, MPI_INT,
		&factor, 1, MPI_INT, 0, MPI_COMM_WORLD);
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
		result *= factor;
	}
	int ProcSum = 0;

	MPI_Gather(&result, 1, MPI_INT,
		&ProcSum, 1, MPI_INT, 0, MPI_COMM_WORLD);

}

int* CollectorProccess()
{
	int val0 = 2;
	int ProcNum, recvMsg;
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	int* sendBuf = new int[ProcNum];
	for (int i = 0; i < ProcNum; i++)
		sendBuf[i] = i;
	MPI_Scatter(sendBuf, 1, MPI_INT,
		&recvMsg, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int ProcSum = 0;
	int* recvBuf = new int[ProcNum];

	MPI_Gather(&ProcSum, 1, MPI_INT,
		recvBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);

	delete[] sendBuf;
	return recvBuf;
}




int main(int argc, char* argv[])
{
	int ProcNum, ProcRank, RecvRank;
	int* result;
	MPI_Status Status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	if (ProcRank == 0)
	{
		result = CollectorProccess();
		for (int i = 1; i < ProcNum; i++)
			cout << "Result(" << i <<"): " << result[i] << endl;
		delete[] result;
	}
	else
		CounterProcess();

	MPI_Finalize();

	return 0;
}


