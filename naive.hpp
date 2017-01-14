#pragma once
#ifndef NAIVE_HPP
#define NAIVE_HPP

#include "main.hpp"

len_t lcs_naive(const std::string& text, const len_t a, const len_t b) {
	const len_t m = std::min(a,b);
	len_t i=0;
	while(m >= i && text[a-i] == text[b-i]) { ++i;}
	return i;
}

len_t lcp_naive(const std::string& text, const len_t a, const len_t b) {
	len_t i=0;
	while(text[a+i] == text[b+i]) { ++i;}
	return i;
}

struct Square {
	len_t start;
	len_t period;
	Square(len_t _start, len_t _period) : start(_start), period(_period) {}
} __attribute__((__packed__));

/**
 * Finds the leftmost occurrences of all squares naively in O(n^3 lg n) time.
 */
Vector<Square> find_naively(const std::string& text) {
	const len_t n = text.length()+1;
	std::set<std::string> check;
	Vector<Square> squares;

	auto report = [&] (const len_t pos, const len_t period) {
		const std::string sub = text.substr(pos,period*2);
		if(check.find(sub) != check.end()) return;
		check.insert(sub);

		DVLOG(1) << "T[" << pos << "," << (pos+period*2-1) << "] = " << text.substr(pos,period) << "," << text.substr(pos+period,period) << " | ";
		for(len_t i = 0; i < period; ++i) { 
			DCHECK_EQ(text[pos+i],text[pos+period+i]);
		}
		squares.emplace_back(pos, period);
	};


	for(len_t i = 0; i < n; ++i) {
		for(len_t p = 1; p <= (n-i)/2; ++p) {
			if(lcp_naive(text,i,i+p) >= p) {
				report(i,p);
			}
		}
	}
	return squares;
}


void check_square_algo(std::string text) {

	Vector<Square> naive_squares { find_naively(text) };
	Vector<Square> squares;
	compute_distinct_squares(text, [&squares] (len_t pos, len_t period) { squares.emplace_back(pos,period); });

	std::sort(squares.begin(), squares.end(), [] (const Square& a, const Square& b) 
			{
			if(a.start == b.start) return a.period < b.period;
			return a.start < b.start; 
			});

	for(size_t s = 0; s < std::min(naive_squares.size(),squares.size()); ++s) {
		if(squares[s].start != naive_squares[s].start || squares[s].period != naive_squares[s].period) {
			std::cout << "Naive : " << naive_squares[s].start << "," << naive_squares[s].period << std::endl;
			std::cout << "Own   : " << squares[s].start << "," << squares[s].period << std::endl;
		}
		// DCHECK_EQ(squares[s].start, naive_squares[s].start);
		// DCHECK_EQ(squares[s].period, naive_squares[s].period);
	}
	DCHECK_EQ(squares.size(), naive_squares.size());

}

#endif /* NAIVE_HPP */
