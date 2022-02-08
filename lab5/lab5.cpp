#include <iostream>
#include <cstdlib>
#include <chrono>
#include <thread>
#include <vector>
#include <string>
#include "mpi.h"

using namespace std;


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

int getRank16(unsigned long long num, int scale)
{
    int base = 1;
    while (num /= scale) base++;
    return base;
}

int power(int x, int y)
{
    int result = 1;
    for (int i = 0; i < y; i++)
        result *= x;
    return result;
}


class EndlessInt
{
    vector<MPI_Datatype> types;
    vector<uint32_t> digits;
    bool negative = false;
    uint64_t limit = UINT32_MAX;


    void deleteExcessDigits()
    {
        for (int i = digits.size() - 1; i >= 0; i--)
        {
            if (digits[i] == 0)
                digits.erase(digits.end() - 1);
            else
                break;
        }
    }

    EndlessInt& plus(const uint32_t& r, int firstDigit)
    {
        while (digits.size() < firstDigit+1) digits.push_back(0);

        int overflow = 0;
        if (digits[firstDigit] > (limit - r))   // overflow
        {
            digits[firstDigit] -= (limit - r);
            overflow = 1;
        }
        else                                // no overflow
        {
            digits[firstDigit] += r;
            overflow = 0;
        }
        int i = firstDigit + 1;
        while (overflow && i < digits.size())
        {
            if (digits[i] == limit)
            {
                digits[i] = 0;
            }
            else
            {
                digits[i]++;
                overflow = 0;
            }
            i++;
        }
        if (overflow) digits.push_back(1);
        return *this;
    }

    EndlessInt& minus(const uint32_t& r)
    {
        if (size() == 1 && r > digits[0])
        {
            uint32_t copy = digits[0];
            digits[0] = r;
            minus(copy);
            negate();
            return *this;
        }

        int underflow = 0;
        if (digits[0] < r)  // underflow
        {
            digits[0] += (limit - r);
            underflow = 1;
        }
        else                // no underflow
        {
            digits[0] -= r;
        }
        int i = 1;
        while (underflow && i < digits.size())
        {
            if (digits[i] == 0)
            {
                digits[i] = limit - underflow;
            }
            else
            {
                digits[i]--;
                underflow = 0;
            }
            i++;
        }
        deleteExcessDigits();
        return *this;
    }

    EndlessInt& plus(const EndlessInt& r)
    {
        int overflow = 0;
        while (size() < r.size()) digits.push_back(0);
        int minSize = size() < r.size() ? size() : r.size();
        for (int i = 0; i < minSize; i++)
        {
            if (digits[i] > (limit - r.digits[i]) || (digits[i] > (limit - r.digits[i] - overflow))) // overflow
            {
                digits[i] -= (limit - r.digits[i]) - overflow;
                overflow = 1;
            }
            else // no overflow
            {
                digits[i] += r.digits[i] + overflow;
                overflow = 0;
            }
        }
        if (overflow) digits.push_back(1);
        return *this;
    }

    EndlessInt& minus(const EndlessInt& r)
    {
        if (r > * this)
        {
            EndlessInt copy(0);
            copy = *this;
            *this = r;
            this->negative = !negative;
            int underflow = 0;
            int minSize = size() < copy.size() ? size() : copy.size();
            for (int i = 0; i < minSize; i++)
            {
                if (digits[i] < copy.digits[i] + underflow) // underflow
                {
                    digits[i] += (limit - (copy.digits[i] + underflow));
                    underflow = 1;
                }
                else // no underflow
                {
                    digits[i] -= copy.digits[i] + underflow;
                    underflow = 0;
                }
            }
        }
        else
        {
            int underflow = 0;
            for (int i = 0; i < size(); i++)
            {
                if (digits[i] < r.digits[i] + underflow) // underflow
                {
                    digits[i] += (limit - (r.digits[i] + underflow));
                    underflow = 1;
                }
                else // no underflow
                {
                    digits[i] -= r.digits[i] + underflow;
                    underflow = 0;
                }
            }
        }
        deleteExcessDigits();
        return *this;
    }

    EndlessInt& operator-=(const uint32_t& r)
    {
        return minus(r);
    }

    EndlessInt& operator+=(const uint32_t& r)
    {
        return plus(r, 0);
    }

public:
    EndlessInt(uint32_t* nums, int size)
    {
        for (int i = 0; i < size; i++)
        {
            digits.push_back(nums[i]);
        }
    }
    EndlessInt(uint32_t num)
    {
        digits.push_back(num);
    }
    EndlessInt()
    {
        digits.push_back(0);
    }

    bool isNegative() const
    {
        return negative;
    }

    int size() const
    {
        return digits.size();
    }

    uint32_t* ptr()
    {
        return &digits[0];
    }

    string toString() const
    {
        string result = "";
        if (negative) result += "-";
        for (int i = digits.size()-1; i >= 0; i--)
        {
            result += "|" + to_string(digits[i]);
        }
        return result;
    }

    ~EndlessInt()
    {
        for (int i = 0; i < types.size(); i++)
        {
            MPI_Type_free(&types[i]);
        }
    }

    void resize(int newSize)
    {
        digits.resize(newSize);
    }

    MPI_Datatype getType()
    {
        MPI_Datatype type;
        MPI_Type_contiguous(size(), MPI_UINT32_T, &type);
        types.push_back(type);
        MPI_Type_commit(&type);
        return type;
    }

    bool operator>(const EndlessInt& r) const
    {
        if (size() == r.size())
        {
            for (int i = size()-1; i >= 0; i--)
            {
                if (digits[i] > r.digits[i])
                    return true;
            }
            return false;
        }
        else
            return size() > r.size();
    }

    EndlessInt operator+(const EndlessInt& r) const
    {
        EndlessInt copy(0);
        copy += r;
        return copy;
    }

    EndlessInt shift(int n) const
    {
        EndlessInt copy(0);

        int shiftedDigits = n / 32;
        int shift = n % 32;
        while (copy.size() < size() - shiftedDigits) 
            copy.digits.push_back(0);

        uint32_t shifted = 0;
        copy.digits[0] = (digits[shiftedDigits] >> shift);
        for (int i = 1; i < size() - shiftedDigits; i++)
        {
            shifted = digits[i + shiftedDigits] << shift;
            copy.digits[i] = (digits[i + shiftedDigits] >> shift);
            copy.digits[i-1] |= shifted;
        }
        copy.deleteExcessDigits();
        return copy;
    }

    EndlessInt& negate()
    {
        negative = !negative;
        return *this;
    }

    EndlessInt& negate(bool sign)
    {
        negative = sign;
        return *this;
    }

    EndlessInt& divideBy2()
    {
        int shifted = 0;
        int preshifted = 0;
        for (int i = size()-1; i >= 0; i--)
        {
            preshifted = digits[i] % 2;
            digits[i] = (digits[i] >> 1);
            if (shifted)
                digits[i] |= 0xF0000000;
            shifted = preshifted;
        }
        deleteExcessDigits();
        return *this;
    }

    EndlessInt& divideBy3()
    {
        int prevSize = size();
        EndlessInt roughMultiplier(0x55555556);
        for (int i = 1; i < size(); i++)
            roughMultiplier.digits.push_back(0x55555556);
        *this *= roughMultiplier;
        *this = shift(32 * prevSize);

        deleteExcessDigits();
        return *this;
    }

    EndlessInt& operator=(const EndlessInt& r)
    {
        digits = r.digits;
        negative = r.negative;
        return *this;
    }

    EndlessInt& operator=(const uint32_t& r)
    {
        digits.clear();
        digits.push_back(r);
        return *this;
    }

    EndlessInt& operator+=(const EndlessInt& r)
    {
        if (negative != r.negative)
            minus(r);
        else 
            plus(r);
        return *this;
    }

    EndlessInt& operator-=(const EndlessInt& r)
    {
        if (negative != r.negative)
            plus(r);
        else
            minus(r);
        return *this;
    }

    EndlessInt& operator*=(const EndlessInt& r)
    {
        negative ^= r.negative;
        uint64_t multiplied = 1;
        EndlessInt result(0);
        for (int i = 0; i < r.size(); i++)
        {
            for (int j = 0; j < size(); j++)
            {
                multiplied = digits[j];
                multiplied *= r.digits[i];
                result.plus(multiplied % (limit + 1), j);
                if (multiplied / (limit + 1) > 0) 
                    result.plus(multiplied / (limit + 1), j+1);
            }
        }
        digits = result.digits;
        return *this;
    }

    friend ostream& operator<<(ostream& os, const EndlessInt& dt)
    {
        return (os << dt.toString());
    }
};

int main(int argc, char* argv[])
{
    int ProcNum, ProcRank;
    const int k = 3;
    uint64_t M = 831275469, N = 897512436;
    uint32_t m[k];
    uint32_t n[k];


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    if (ProcNum != 4) return 1;
    //if (ProcRank == 0)
    //{
    //    cout << "Enter integer numbers to multiply" << endl;
    //    cin >> M >> N;
    //}
    //MPI_Bcast(&M, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    //MPI_Bcast(&N, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);

    long long int maxNum = M > N ? M : N;

    int scale = 10;
    int rank = getRank16(maxNum, scale);
    int base = rank / 3;
    if (rank > base * 3) base++;

    m[0] = M / power(scale, 2 * base);
    n[0] = N / power(scale, 2 * base);
    m[1] = M / power(scale, base) % power(scale, base);
    n[1] = N / power(scale, base) % power(scale, base);
    m[2] = M % power(scale, base);
    n[2] = N % power(scale, base);


    MPI_Group defaultGroup, group13, group12, group23, group02, group023;
    MPI_Comm  comm13, comm12, comm23, comm02, comm023;
    MPI_Comm_group(MPI_COMM_WORLD, &defaultGroup);
    int groupExludes13[] = { 0,2 };
    int groupExludes03[] = { 0,3 };
    int groupExludes23[] = { 0,1 };
    int groupExludes02[] = { 1,3 };
    int groupExludes023[] = { 1 };
    MPI_Group_excl(defaultGroup, 2, groupExludes13, &group13);
    MPI_Group_excl(defaultGroup, 2, groupExludes03, &group12);
    MPI_Group_excl(defaultGroup, 2, groupExludes23, &group23);
    MPI_Group_excl(defaultGroup, 2, groupExludes02, &group02);
    MPI_Group_excl(defaultGroup, 1, groupExludes023, &group023);
    MPI_Comm_create(MPI_COMM_WORLD, group13, &comm13);
    MPI_Comm_create(MPI_COMM_WORLD, group12, &comm12);
    MPI_Comm_create(MPI_COMM_WORLD, group23, &comm23);
    MPI_Comm_create(MPI_COMM_WORLD, group02, &comm02);
    MPI_Comm_create(MPI_COMM_WORLD, group023, &comm023);
    MPI_Group_free(&defaultGroup);
    MPI_Group_free(&group13);
    MPI_Group_free(&group12);
    MPI_Group_free(&group23);
    MPI_Group_free(&group02);
    MPI_Group_free(&group023);


    EndlessInt p[5], q[5];
    EndlessInt result;
    int size = 0;
    int negative = false;

    switch (ProcRank)
    {
    case 0:
        p[0] = m[0];
        q[0] = n[0];
        p[0] *= q[0]; //r0
        p[4] = m[2];
        q[4] = n[2];
        p[4] *= q[4]; //r4

        size = p[0].size();
        MPI_Bcast(&size, 1, MPI_INT, 0, comm02);
        negative = p[0].isNegative();
        MPI_Bcast(&negative, 1, MPI_INT, 0, comm02);
        MPI_Bcast(p[0].ptr(), 1, p[0].getType(), 0, comm02); //send

        //cout << p[0] << " p0 send" << endl;
        //cout << p[4] << " p4 send" << endl;

        size = p[4].size();
        MPI_Bcast(&size, 1, MPI_INT, 0, comm023);
        negative = p[4].isNegative();
        MPI_Bcast(&negative, 1, MPI_INT, 0, comm023);
        MPI_Bcast(p[4].ptr(), 1, p[4].getType(), 0, comm023); //send

        //MPI_Reduce(p[0].ptr(), &result, 1, p[0].getType(), MPI_SUM, 0, MPI_COMM_WORLD);
        //cout << result << endl;
        MPI_Barrier(MPI_COMM_WORLD);
        cout << p[0] << endl;
        cout << p[4] << endl;
        break;
    case 1:
        p[1] = m[0] + m[1] + m[2];
        q[1] = n[0] + n[1] + m[2];
        p[1] *= q[1]; //r1

        size = p[1].size();
        MPI_Bcast(&size, 1, MPI_INT, 0, comm13);
        negative = p[1].isNegative();
        MPI_Bcast(&negative, 1, MPI_INT, 0, comm13);
        MPI_Bcast(p[1].ptr(), 1, p[1].getType(), 0, comm13); //send

        //cout << p[1] << " p1 send" << endl;

        MPI_Bcast(&size, 1, MPI_INT, 1, comm12);
        p[2].resize(size);
        MPI_Bcast(&negative, 1, MPI_INT, 0, comm12);
        p[2].negate(negative);
        MPI_Bcast(p[2].ptr(), 1, p[2].getType(), 1, comm12); //recv

        //cout << p[2] << " p2 recv" << endl;
        cout << p[1] << " " << p[2] << endl;
        p[1] -= p[2];
        cout << "p1 - p2: " << p[1] << endl;
        p[1].divideBy2();
        cout << "(p1 - p2)/2: " << p[1] << endl;

        size = p[1].size();
        MPI_Bcast(&size, 1, MPI_INT, 0, comm12);
        negative = p[1].isNegative();
        MPI_Bcast(&negative, 1, MPI_INT, 0, comm12);
        MPI_Bcast(p[1].ptr(), 1, p[1].getType(), 0, comm12); //send

        //cout << p[1] << " p1.1 send" << endl;

        MPI_Bcast(&size, 1, MPI_INT, 1, comm13);
        p[3].resize(size);
        MPI_Bcast(&negative, 1, MPI_INT, 1, comm13);
        p[3].negate(negative);
        MPI_Bcast(p[3].ptr(), 1, p[3].getType(), 1, comm13); //recv

        //cout << p[3] << " p3 recv" << endl;
        cout << p[1] << " " << p[3] << endl;
        p[1] -= p[3];
        cout << "p1 - p3: " << p[1] << endl;
        //MPI_Reduce(p[1].ptr(), 0, 1, p[1].getType(), MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        cout << p[1] << endl;
        break;
    case 2:
        p[2] = m[0] - m[1] + m[2];
        q[2] = n[0] - n[1] + n[2];
        p[2] *= q[2]; //r2

        size = p[2].size();
        MPI_Bcast(&size, 1, MPI_INT, 1, comm12);
        negative = p[2].isNegative();
        MPI_Bcast(&negative, 1, MPI_INT, 1, comm12);
        MPI_Bcast(p[2].ptr(), 1, p[2].getType(), 1, comm12); //send

        cout << p[2] << " p2 send" << endl;

        MPI_Bcast(&size, 1, MPI_INT, 0, comm02);
        p[0].resize(size);
        MPI_Bcast(&negative, 1, MPI_INT, 0, comm02);
        p[0].negate(negative);
        MPI_Bcast(p[0].ptr(), 1, p[0].getType(), 0, comm02); //recv

        cout << p[0] << " p0 recv" << endl;
        cout << p[2] << " " << p[0] << endl;
        p[2] -= p[0];
        cout << "p2 - p0: " << p[2] << endl;

        size = p[2].size();
        MPI_Bcast(&size, 1, MPI_INT, 0, comm23);
        negative = p[2].isNegative();
        MPI_Bcast(&negative, 1, MPI_INT, 0, comm23);
        MPI_Bcast(p[2].ptr(), 1, p[2].getType(), 0, comm23); //send

        cout << p[2] << " p2 send" << endl;

        MPI_Bcast(&size, 1, MPI_INT, 0, comm023);
        p[4].resize(size);
        MPI_Bcast(&negative, 1, MPI_INT, 0, comm023);
        p[4].negate(negative);
        MPI_Bcast(p[4].ptr(), 1, p[4].getType(), 0, comm023); //recv

        cout << p[4] << " p4 recv" << endl;
        cout << "aaaaaa" << endl;
        MPI_Bcast(&size, 1, MPI_INT, 0, comm12);
        cout << size << endl;
        p[1].resize(size);
        cout << "aaaaaa" << endl;
        MPI_Bcast(&negative, 1, MPI_INT, 0, comm12);
        cout << "aaaaaa" << endl;
        p[1].negate(negative);
        cout << "aaaaaa" << endl;
        MPI_Bcast(p[1].ptr(), 1, p[1].getType(), 0, comm12); //recv

        cout << p[1] << " p1.1 recv" << endl;
        cout << p[2] << " " << p[1] << " " << p[4] << endl;
        p[2] += p[1];
        cout << "p2 + p1: " << p[2] << endl;
        p[2] -= p[4];
        cout << "p2 - p4: " << p[2] << endl;
        //MPI_Reduce(p[2].ptr(), 0, 1, p[2].getType(), MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        cout << p[2] << endl;
        break;
    case 3:
        p[3] = m[0] - 2 * m[1] + 4 * m[2];
        q[3] = n[0] - 2 * n[1] + 4 * n[2];
        p[3] *= q[3]; //r3

        MPI_Bcast(&size, 1, MPI_INT, 0, comm13);
        p[1].resize(size);
        MPI_Bcast(&negative, 1, MPI_INT, 0, comm13);
        p[1].negate(negative);
        MPI_Bcast(p[1].ptr(), 1, p[1].getType(), 0, comm13); //recv

        //cout << p[1] << " p1 recv" << endl;
        cout << p[3] << " " << p[1] << endl;
        p[3] -= p[1];
        cout << "p3 - p1: " << p[3] << endl;
        p[3].divideBy3();
        cout << "(p3 - p1)/3: " << p[3] << endl;

        MPI_Bcast(&size, 1, MPI_INT, 0, comm023);
        p[4].resize(size);
        MPI_Bcast(&negative, 1, MPI_INT, 0, comm023);
        p[4].negate(negative);
        MPI_Bcast(p[4].ptr(), 1, p[4].getType(), 0, comm023); //recv

        //cout << p[4] << " p4 recv" << endl;

        MPI_Bcast(&size, 1, MPI_INT, 0, comm23);
        p[2].resize(size);
        MPI_Bcast(&negative, 1, MPI_INT, 0, comm23);
        p[2].negate(negative);
        MPI_Bcast(p[2].ptr(), 1, p[2].getType(), 0, comm23); //recv

        //cout << p[2] << " p2 recv" << endl;
        cout << p[3] << " " << p[2] << " " << p[4] << endl;
        p[3] -= p[2];
        cout << "p3 - p2: " << p[3] << endl;
        p[3].negate();
        cout << "-(p3 - p2): " << p[3] << endl;
        p[3].divideBy2();
        cout << "-(p3 - p2)/2: " << p[3] << endl;
        p[4] *= EndlessInt(2);
        cout << "p4*2: " << p[4] << endl;
        p[3] += p[4];
        cout << "-(p3 - p2)/2 + 2*p4: " << p[3] << endl;

        size = p[3].size();
        MPI_Bcast(&size, 1, MPI_INT, 1, comm13);
        negative = p[3].isNegative();
        MPI_Bcast(&negative, 1, MPI_INT, 1, comm13);
        MPI_Bcast(p[3].ptr(), 1, p[3].getType(), 1, comm13); //send

        //cout << p[3] << " p3 send" << endl;
        //MPI_Reduce(p[3].ptr(), 0, 1, p[3].getType(), MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        cout << p[3] << endl;
        break;
    }

    //MPI_Comm_free(&comm13);
    //MPI_Comm_free(&comm12);
    //MPI_Comm_free(&comm23);
    //MPI_Comm_free(&comm02);
    //MPI_Comm_free(&comm023);
    for (int i = 0; i < 5; i++)
    {
        p[i].~EndlessInt();
        q[i].~EndlessInt();
    }
    result.~EndlessInt();
    MPI_Finalize();

    return 0;
}



