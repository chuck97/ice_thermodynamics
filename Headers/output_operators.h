#pragma once
#include <iostream>

template <typename T>
std::ostream& operator<<(std::ostream& stream, const std::vector<T>& vec)
{
    for (const auto& item: vec)
    {
        stream << item << " "; 
    }
    return stream;
};