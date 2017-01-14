#pragma once
#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <vector>
#include <stdexcept>

template<class T>
class Vector : public std::vector<T>
{
	public:
		using std::vector<T>::vector;
#ifndef NDEBUG
	T& operator[](size_t n) /*override*/ { 
		DCHECK_LT(n, this->size());
		return std::vector<T>::at(n); 
	}
	const T& operator[](size_t n) const /*override*/ { 
		DCHECK_LT(n, this->size());
		return std::vector<T>::at(n);
	}
#endif
};

#endif /* VECTOR_HPP */
