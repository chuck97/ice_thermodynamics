#pragma once
#include <iostream>

#define ERR(message) {std::cerr << "Error: " << message  << std::endl; exit(1);}

template <typename T>
std::ostream& operator<<(std::ostream& stream, const std::vector<T>& vec)
{
    for (const auto& item: vec)
    {
        stream << item << " "; 
    }
    return stream;
};

template <typename T>
std::vector<T> operator-(const std::vector<T>& vec1, const std::vector<T>& vec2)
{
    if (vec1.size() != vec2.size())
    {
        ERR("can't subtract vectors with different sizes");
    }

    std::vector<T> res(vec1.size());

    for (unsigned int i = 0; i < vec1.size(); ++i)
    {
        res[i] = vec1[i] - vec2[i];
    }

    return res;
};

template <typename T>
std::vector<T> operator+(const std::vector<T>& vec1, const std::vector<T>& vec2)
{
    if (vec1.size() != vec2.size())
    {
        ERR("can't add vectors with different sizes");
    }

    std::vector<T> res(vec1.size());

    for (unsigned int i = 0; i < vec1.size(); ++i)
    {
        res[i] = vec1[i] + vec2[i];
    }

    return res;
};