#include <iostream>
#include <cstdlib>
#include <chrono>
#include <thread>
#include <vector>
#include "mpi.h"
#include <complex>
 
using namespace std;
 
const double PI = 3.141592653589793238460;
 
 
class PolynomialComplex
{
    vector<MPI_Datatype> types;
 
public:
    ~PolynomialComplex()
    {
        for (int i = 0; i < types.size(); i++)
        {
            MPI_Type_free(&types[i]);
        }
    }
 
    MPI_Datatype get(int count)
    {
        MPI_Datatype type;
        MPI_Type_contiguous(count, MPI_DOUBLE_COMPLEX, &type);
        types.push_back(type);
        MPI_Type_commit(&type);
        return type;
    }
};
/// ////////////////////////////////////////////////////////
 
 
typedef complex<double> Complex;
typedef vector<complex<double>> vectorComplex;
 
vectorComplex fft(const vectorComplex& as)
{
    int n = as.size();
 
    if (n == 1) return vectorComplex(1, as[0]);
 
    vectorComplex w(n);
    for (int i = 0; i < n; i++)
    {
        double alpha = 2 * PI * i / n;
        w[i] = Complex(cos(alpha), sin(alpha));
    }
 
    vectorComplex A(n / 2), B(n / 2);
    for (int i = 0; i < n / 2; i++)
    {
        A[i] = as[i * 2];
        B[i] = as[i * 2 + 1];
    }
    vectorComplex Av = fft(A);
    vectorComplex Bv = fft(B);
    vectorComplex res(n);
    for (int i = 0; i < n; i++)
        res[i] = Av[i % (n / 2)] + w[i] * Bv[i % (n / 2)];
    return res;
}
 
vectorComplex fftPartial(const vectorComplex& A, const vectorComplex& B)
{
    int n = A.size() + B.size();
 
    vectorComplex w(n);
    for (int i = 0; i < n; i++) 
    {
        double alpha = 2 * PI * i / n;
        w[i] = Complex(cos(alpha), sin(alpha));
    }
 
    vectorComplex Av = A;
    vectorComplex Bv = B;
    vectorComplex res(n);
    for (int i = 0; i < n; i++)
        res[i] = Av[i % (n / 2)] + w[i] * Bv[i % (n / 2)];
    return res;
}
 
vectorComplex ifft(const vectorComplex& as)
{
    vectorComplex res = fft(as);
    for (int i = 0; i < (int)res.size(); i++) 
        res[i] /= as.size();
    reverse(res.begin() + 1, res.end());
    return res;
}
 
vectorComplex ifftPartial(const vectorComplex& A, const vectorComplex& B)
{
    vectorComplex res = fftPartial(A, B);
    for (int i = 0; i < (int)res.size(); i++)
        res[i] /= res.size();
    reverse(res.begin() + 1, res.end());
    return res;
}
 
vector<int> carry(const vectorComplex& poly)
{
    vector<int> res;
    for (int i = poly.size() - 1; i >= 0; i--)
        res.push_back(round(poly[i].real()));
    for (int i = 0; i < res.size() - 1; i++)
    {
        res[i + 1] += res[i] / 10;
        res[i] %= 10;
    }
    while (res[res.size() - 1] > 9)
    {
        res.push_back(res[res.size() - 1] / 10);
        res[res.size() - 2] %= 10;
    }
    for (int i = res.size() - 1; i > 0; i--)
    {
        if (res[i] == 0)
            res.erase(res.begin() + i);
        else
            break;
    }
    return res;
}
 
int main(int argc, char* argv[])
{
    int ProcNum, ProcRank;
    double* result;
 
    string N, M;
    vectorComplex poly, polyHalf, polyHalfRight;
    vector<int> poly1, poly2;
    Complex *fftPoly1, *fftPoly2;
 
    MPI_Init(&argc, &argv);
    int index[] = { 2,4,5,6 };
    int edges[] = { 2,1, 3,0, 0, 1 };
    MPI_Comm SSComm;
    MPI_Status Status;
    MPI_Graph_create(MPI_COMM_WORLD, 4, index, edges, 1, &SSComm);
 
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    PolynomialComplex mpiPolyComplex;
    int NSize, MSize;
 
    if (ProcRank == 0)
    {
        cout << "Enter numbers to multiply" << endl;
        cin >> N >> M;
        NSize = N.size();
        MSize = M.size();
    }
    char* buffN, * buffM;
    MPI_Bcast(&NSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&MSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (ProcRank == 0)
    {
        buffN = (char*)N.c_str();
        buffM = (char*)M.c_str();
    }
    else
    {
        buffN = new char[NSize];
        buffM = new char[MSize];
    }
    MPI_Bcast(buffN, NSize, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(buffM, MSize, MPI_CHAR, 0, MPI_COMM_WORLD);
    
    if (ProcRank != 0)
    {
        for (int i = 0; i < NSize; i++)
            N.push_back(buffN[i]);
        for (int i = 0; i < MSize; i++)
            M.push_back(buffM[i]);
    }
 
    size_t prodSize = N.size() + M.size();
    size_t twoSize = 4;
    while (twoSize < prodSize) twoSize *= 2;
    prodSize = twoSize;
    string zeroes;
    while (zeroes.size() < prodSize - N.size()) zeroes.push_back('0');
    N = zeroes + N;
    zeroes.clear();
    while (zeroes.size() < prodSize - M.size()) zeroes.push_back('0');
    M = zeroes + M;
 
    if (ProcRank == 0)
    {
        for (int i = 0; i < N.size() / 2; i++)
        {
            polyHalf.push_back(Complex((int)(N[2 * i] - '0')));
        }
        polyHalf = fft(polyHalf);
        polyHalfRight.resize(polyHalf.size());
        MPI_Recv(polyHalfRight.data(), 1, mpiPolyComplex.get(polyHalf.size()), 2, 0, SSComm, &Status);
        poly = fftPartial(polyHalf, polyHalfRight);
        vectorComplex polyRight(poly.size());
        MPI_Recv(polyRight.data(), 1, mpiPolyComplex.get(polyRight.size()), 1, 0, SSComm, &Status);
        for (int i = 0; i < polyRight.size(); i++)
            poly[i] *= polyRight[i];
        for (int i = 0; i < poly.size() / 2; i++)
        {
            polyHalf[i] = poly[2 * i];
            polyHalfRight[i] = poly[2 * i + 1];
        }
        MPI_Send(polyHalfRight.data(), 1, mpiPolyComplex.get(polyHalfRight.size()), 1, 0, SSComm);
 
        vectorComplex polyHalfHalf(polyHalf.size() / 2);
        vectorComplex polyHalfHalfRight(polyHalf.size() / 2);
        for (int i = 0; i < polyHalf.size() / 2; i++)
        {
            polyHalfHalf[i] = polyHalf[2 * i];
            polyHalfHalfRight[i] = polyHalf[2 * i + 1];
        }
        MPI_Send(polyHalfHalfRight.data(), 1, mpiPolyComplex.get(polyHalfHalfRight.size()), 2, 0, SSComm);
 
        polyHalfHalf = fft(polyHalfHalf);
        MPI_Recv(polyHalfHalfRight.data(), 1, mpiPolyComplex.get(polyHalfHalf.size()), 2, 0, SSComm, &Status);
        polyHalf = fftPartial(polyHalfHalf, polyHalfHalfRight);
        MPI_Recv(polyHalfRight.data(), 1, mpiPolyComplex.get(polyHalf.size()), 1, 0, SSComm, &Status);
        poly = ifftPartial(polyHalf, polyHalfRight);
 
        poly.erase(poly.end()-1);
        vector<int> result = carry(poly);
        for (int i = result.size() - 1; i >= 0; i--)
            cout << result[i];
        cout << endl;
    }
    else if (ProcRank == 2)
    {
        for (int i = 0; i < N.size() / 2; i++)
        {
            polyHalf.push_back(Complex((int)(N[2 * i + 1] - '0')));
        }
        polyHalf = fft(polyHalf);
        MPI_Send(polyHalf.data(), 1, mpiPolyComplex.get(polyHalf.size()), 0, 0, SSComm);
 
        vectorComplex polyHalfHalf(polyHalf.size() / 2);
        MPI_Recv(polyHalfHalf.data(), 1, mpiPolyComplex.get(polyHalfHalf.size()), 0, 0, SSComm, &Status);
        polyHalfHalf = fft(polyHalfHalf);
        MPI_Send(polyHalfHalf.data(), 1, mpiPolyComplex.get(polyHalfHalf.size()), 0, 0, SSComm);
    }
    else if (ProcRank == 1)
    {
        for (int i = 0; i < M.size() / 2; i++)
        {
            polyHalf.push_back(Complex((int)(M[2 * i] - '0')));
        }
        polyHalf = fft(polyHalf);
        polyHalfRight.resize(polyHalf.size());
        MPI_Recv(polyHalfRight.data(), 1, mpiPolyComplex.get(polyHalf.size()), 3, 0, SSComm, &Status);
        poly = fftPartial(polyHalf, polyHalfRight);
        MPI_Send(poly.data(), 1, mpiPolyComplex.get(poly.size()), 0, 0, SSComm);
 
        MPI_Recv(polyHalf.data(), 1, mpiPolyComplex.get(polyHalf.size()), 0, 0, SSComm, &Status);
        vectorComplex polyHalfHalf(polyHalf.size() / 2);
        vectorComplex polyHalfHalfRight(polyHalf.size() / 2);
        for (int i = 0; i < polyHalf.size() / 2; i++)
        {
            polyHalfHalf[i] = polyHalf[2 * i];
            polyHalfHalfRight[i] = polyHalf[2 * i + 1];
        }
        MPI_Send(polyHalfHalfRight.data(), 1, mpiPolyComplex.get(polyHalfHalfRight.size()), 3, 0, SSComm);
        polyHalfHalf = fft(polyHalfHalf);
        MPI_Recv(polyHalfHalfRight.data(), 1, mpiPolyComplex.get(polyHalfHalfRight.size()), 3, 0, SSComm, &Status);
        polyHalf = fftPartial(polyHalfHalf, polyHalfHalfRight);
        MPI_Send(polyHalf.data(), 1, mpiPolyComplex.get(polyHalf.size()), 0, 0, SSComm);
    }
    else if (ProcRank == 3) 
    {
        for (int i = 0; i < M.size() / 2 ; i++)
        {
            polyHalf.push_back(Complex((int)(M[2 * i + 1] - '0')));
        }
        polyHalf = fft(polyHalf);
        MPI_Send(polyHalf.data(), 1, mpiPolyComplex.get(polyHalf.size()), 1, 0, SSComm);
 
        vectorComplex polyHalfHalf(polyHalf.size() / 2);
        MPI_Recv(polyHalfHalf.data(), 1, mpiPolyComplex.get(polyHalfHalf.size()), 1, 0, SSComm, &Status);
        polyHalfHalf = fft(polyHalfHalf);
        MPI_Send(polyHalfHalf.data(), 1, mpiPolyComplex.get(polyHalfHalf.size()), 1, 0, SSComm);
    }
 
    MPI_Barrier(MPI_COMM_WORLD);
    mpiPolyComplex.~PolynomialComplex();
    MPI_Finalize();
 
    return 0;
}