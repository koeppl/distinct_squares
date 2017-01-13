#include <sdsl/suffix_trees.hpp>
#include "naive.hpp"


using namespace sdsl;

/**
 * Naively checks whether the lpf array is correct
 */
template<class lpf_t>
void check_lpf(const std::string& text, const lpf_t& lpf) {
	for(size_t i = 0; i < lpf.size(); ++i) {
		{
			const std::string pattern = text.substr(i, lpf[i]);
			const size_t firstMatch = text.find(pattern,0);
			DCHECK_NE(firstMatch, std::string::npos);
			DCHECK_LE(firstMatch, i);
		}
		if(i+1+lpf[i] < text.length())
		{
			const std::string pattern = text.substr(i, 1+lpf[i]);
			const size_t firstMatch = text.find(pattern,0);
			DCHECK_NE(firstMatch, std::string::npos);
			DCHECK_EQ(firstMatch, i);
		}
	}
}


void check_lpf(const std::string& text) {
	const len_t n = text.length()+1;
	Vector<len_t> sa(n);
	suffix_sort((const uint8_t*)text.c_str(), sa.data(), n);

	Vector<len_t> plcp { construct_phi_array<Vector<len_t>,decltype(sa)>(sa) };
	phi_algorithm(plcp,text);
	const Vector<len_t> lcp { construct_lcp_array<Vector<len_t>,decltype(plcp),decltype(sa)>(plcp,sa) };
	Vector<len_t>().swap(plcp); // delete plcp
	const sdsl::rmq_succinct_sct<> rmqlcp { &lcp };

	Vector<len_t> isa(n);
	for(size_t i = 0; i < n; ++i) {
		DCHECK_LT(sa[i], n);
		isa[sa[i]] = i;
	}

	const Vector<len_t> lpf = gen_lpf<Vector<len_t>,decltype(lcp),decltype(isa)>(lcp, isa);
	check_lpf(text, lpf);
}

/**
 * Checks with data structures of the SDSL if
 * the built data structures like LCP are correct.
 */
void check_toolchain(const std::string& text) {
	typedef sdsl::cst_sct3<> cst_t;
    cst_t cst;
    construct_im(cst, text, 1);
	const Vector<len_t> lpf = gen_lpf<Vector<len_t>,decltype(cst.lcp),decltype(cst.csa.isa)>(cst.lcp, cst.csa.isa);

	const len_t n = text.length()+1;
	Vector<len_t> sa(n);
	suffix_sort((const uint8_t*)text.c_str(), sa.data(), n);

	Vector<len_t> plcp { construct_phi_array<Vector<len_t>,decltype(sa)>(sa) };
	phi_algorithm(plcp,text);
	Vector<len_t> lcp { construct_lcp_array<Vector<len_t>,decltype(plcp),decltype(sa)>(plcp,sa) };
	Vector<len_t>().swap(plcp); // delete plcp
	const rmq_succinct_sct<> rmqlcp { &lcp };

	Vector<len_t> isa(n);
	for(size_t i = 0; i < n; ++i) {
		DCHECK_LT(sa[i], n);
		isa[sa[i]] = i;
	}
	Vector<len_t>().swap(sa); //delete sa

	const Vector<len_t> lpf2 = gen_lpf<Vector<len_t>,decltype(lcp),decltype(isa)>(lcp, isa);
	for(size_t i = 0; i < lpf.size(); ++i) {
		DCHECK_EQ(lpf[i], lpf2[i]);
	}
		check_lpf(text, lpf);
}


void check_rmq_on_lcp(std::string text) {
	const len_t n = text.length()+1;
	Vector<len_t> sa(n);
	suffix_sort((const uint8_t*)text.c_str(), sa.data(), n);

	Vector<len_t> plcp { construct_phi_array<Vector<len_t>,decltype(sa)>(sa) };
	phi_algorithm(plcp,text);
	const Vector<len_t> lcp { construct_lcp_array<Vector<len_t>,decltype(plcp),decltype(sa)>(plcp,sa) };
	Vector<len_t>().swap(plcp); // delete plcp
	const sdsl::rmq_succinct_sct<> rmqlcp { &lcp };

	Vector<len_t> isa(n);
	for(size_t i = 0; i < n; ++i) {
		DCHECK_LT(sa[i], n);
		isa[sa[i]] = i;
	}

	const Vector<len_t> lpf = gen_lpf<Vector<len_t>,decltype(lcp),decltype(isa)>(lcp, isa);
	const sdsl::rmq_succinct_sct<> rmqlpf { &lpf };

	//generate LCP^-1
	std::reverse(text.begin(), text.end());
	suffix_sort((const uint8_t*)text.c_str(), sa.data(), n);

	plcp = construct_phi_array<Vector<len_t>,decltype(sa)>(sa);
	phi_algorithm(plcp,text);
	Vector<len_t> lcs = construct_lcp_array<Vector<len_t>,decltype(plcp),decltype(sa)>(plcp,sa);
	Vector<len_t>().swap(plcp); // delete plcp
	const sdsl::rmq_succinct_sct<> rmqlcs { &lcs };

	Vector<len_t> isai(n);
	for(size_t i = 0; i < n; ++i) {
		DCHECK_LT(sa[i], n);
		isai[sa[i]] = i;
	}
	Vector<len_t>().swap(sa); //delete sa
	std::reverse(text.begin(), text.end());



	auto lcpq = [&isa,&lcp,&rmqlcp] (const len_t a, const len_t b) {
		return lcp_rmq(isa,lcp,rmqlcp,a,b);
	};
	auto lcsq = [&n,&isai,&lcs,&rmqlcs] (const len_t a, const len_t b) {
		return lcs_rmq(isai,lcs,rmqlcs,a,b);
	};


	for(len_t i = 0; i < text.length();++i) {
		for(len_t j = 0; j < text.length();++j) {
			if(j==i) continue;
			DCHECK_EQ(lcp_naive(text,i,j), lcpq(i,j));
			DCHECK_EQ(lcs_naive(text,i,j), lcsq(i,j));
		}
	}
}


#include <random>
namespace Ranges {
	std::pair<size_t,size_t> binary(65,66);
	std::pair<size_t,size_t> ternary(65,67);
	std::pair<size_t,size_t> numbers(48,57);
	std::pair<size_t,size_t> printable(33,123);
};

std::string random_uniform(const size_t length, const std::pair<size_t,size_t> range = Ranges::numbers, size_t seed = 0) {
	std::string s(length,0);
	std::default_random_engine engine(seed);
	std::uniform_int_distribution<char> dist(range.first,range.second);
	for(size_t i = 0; i < length; ++i) {
		s[i] = dist(engine);
	}
	return s;
}

#include <gtest/gtest.h>

TEST(algo, test_algo) {
	for(size_t i = 0; i < 10; ++i) {
	for(size_t j = 2; j < 40; ++j) {
		check_square_algo(random_uniform(j,Ranges::binary, 0));
		check_square_algo(random_uniform(j,Ranges::ternary, 0));
		check_square_algo(random_uniform(j,Ranges::numbers, 0));
	}
	}
}

TEST(DS, check_toolchain) {
	for(size_t j = 2; j < 50; ++j) {
		check_toolchain(random_uniform(j,Ranges::binary, 0));
		check_toolchain(random_uniform(j,Ranges::ternary, 0));
		check_toolchain(random_uniform(j,Ranges::numbers, 0));
	}
}

TEST(DS, check_lpf) {
	for(size_t i = 0; i < 10; ++i) {
		for(size_t j = 2; j < 100; ++j) {
			check_lpf(random_uniform(j,Ranges::binary, 0));
			check_lpf(random_uniform(j,Ranges::ternary, 0));
			check_lpf(random_uniform(j,Ranges::numbers, 0));
		}
	}
}


TEST(DS, test_rmq) {
	for(size_t j = 2; j < 50; ++j) {
		check_rmq_on_lcp(random_uniform(j,Ranges::binary, 0));
		check_rmq_on_lcp(random_uniform(j,Ranges::ternary, 0));
		check_rmq_on_lcp(random_uniform(j,Ranges::numbers, 0));
	}
}


