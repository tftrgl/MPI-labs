#include <iostream>
#include <cstdlib>
#include <chrono>
#include <thread>
#include <vector>
#include "mpi.h"

using namespace std;

class Polynomial
{
	vector<MPI_Datatype> types;

public:
	~Polynomial() 
	{
		for (int i = 0; i < types.size(); i++)
		{
			MPI_Type_free(&types[i]);
		}
	}

	MPI_Datatype& get(int count)
	{
		MPI_Datatype type;
		MPI_Type_contiguous(count, MPI_DOUBLE, &type);
		types.push_back(type);
		MPI_Type_commit(&type);
		return type;
	}
};


void mySumFunction(void* inputBuffer, void* outputBuffer, int* len, MPI_Datatype* datatype)
{
	double* input = (double*)inputBuffer;
	double* output = (double*)outputBuffer;
	MPI_Status status;
	int size = 0;
	MPI_Type_size(*datatype, &size);

	int count = size / sizeof(double);
	for (int i = 0; i < count; i++)
	{
		output[i] += input[i];
	}
}


int main(int argc, char* argv[])
{
	int ProcNum, ProcRank;
	double* result;

	int N, M;
	double *poly1, *poly2;

	double *minArray, *maxArray;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	Polynomial mpiPoly;
	if (ProcRank == 0)
	{
		cout << "Polynomial ranks (N,M)" << endl;
		cin >> N >> M;
	}
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
	poly1 = new double[N];
	poly2 = new double[M];
	MPI_Datatype MyPolyType = mpiPoly.get(N + M - 1);
	MPI_Datatype PolyN = mpiPoly.get(N);
	MPI_Datatype PolyM = mpiPoly.get(M);
	if (ProcRank == 0)
	{
		cout << "Polynomial 1:" << endl;
		for (int i = N - 1; i >= 0; i--)
			cin >> poly1[i];
		cout << "Polynomial 2:" << endl;
		for (int i = M - 1; i >= 0; i--)
			cin >> poly2[i];
	}
	MPI_Bcast(poly1, 1, PolyN, 0, MPI_COMM_WORLD);
	MPI_Bcast(poly2, 1, PolyM, 0, MPI_COMM_WORLD);

	int minSize = min(N, M);
	int maxSize = max(N, M);
	if (minSize == N)
	{
		minArray = poly1;
		maxArray = poly2;
	}
	else
	{
		minArray = poly2;
		maxArray = poly1;
	}

	result = new double[N + M -1];
	for (int i = 0; i < N + M -1; i++)
		result[i] = 0;
	for (int i = ProcRank / (ProcNum-1) * minSize; i < (ProcRank + 1) / (ProcNum-1) * minSize; i++)
	{
		for (int j = 0; j < maxSize; j++)
		{
			result[i + j] += minArray[i] * maxArray[j];
		}
	}

	MPI_Op MySum;
	MPI_Op_create(&mySumFunction, 1, &MySum);
	double* finalResult = new double[N + M -1];
	for (int i = 0; i < N + M -1; i++)
		finalResult[i] = 0;
	MPI_Reduce(result, finalResult, 1, MyPolyType, MySum, 0, MPI_COMM_WORLD);

	
	if (ProcRank == 0)
	{
		cout << "Final result:\n";
		for (int i = N + M - 2; i >= 0; i--)
		{
			cout << finalResult[i] << " ";
		}
		cout << endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Op_free(&MySum);
	mpiPoly.~Polynomial();
	MPI_Finalize();

	delete[] result, finalResult, poly1, poly2;
	return 0;
}


