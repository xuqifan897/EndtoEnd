#ifndef utils_h
#define utils_h

#include <array>
#include <iostream>

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& arr)
{
    os << "[";
    int N = arr.size();
    for (int i=0; i<N; i++)
    {
        os << arr[i];
        if (i < N-1)
            os << ", ";
    }
    os << "]";
    return os;
}

#endif