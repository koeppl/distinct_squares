#pragma once

#include <iostream>
#include "/scripts/code/vector.hpp"
#include <string>
#include <fstream>
#include <streambuf>
#include <divsufsort.h>
#include <divsufsort64.h>
#include <sdsl/int_vector.hpp>
#include <sdsl/rmq_succinct_sada.hpp>
#include <sdsl/rmq_succinct_sct.hpp>
#include <glog/logging.h>


#define DEBUG
#ifdef DEBUG
#define DODEBUG(x) x
#else
#define DODEBUG(x)
#endif

#ifdef DNDEBUG
#ifdef DCHECK 
#undef DCHECK
#define DCHECK(x) 
#endif//DCHECK
#ifdef VLOG_IS_ON
#undef VLOG_IS_ON
#define VLOG_IS_ON(x) false
#endif//VLOG_IS_ON
#endif//DNDEBUG



typedef uint32_t len_t; //! this type should be large enough to address all text positions

constexpr len_t RMQ_THRES_LCP = 10000000;
constexpr len_t RMQ_THRES_TRYLCP = 1000000;
constexpr len_t RMQ_THRES_LPF = 1000000;

/**
 * Given a not yet reported leftmost occurence of a square we right-rotate it
 * to find all leftmost occurrences of squares that can be formed by right-rotating the currently found one
 * Runs in O(p) time
 * @text the input text
 * @param s the starting position of the found square
 * @param p the period
 * @param lpf the LPF array
 * @param marker marks the starting position of a found square, helps to prevent reporting the same occurrence twice
 * @param function of type (len_t pos, len_t period) that is called when a square at position pos is found
 */
template<class lpf_t, class report_t>
void linear_rightrotate(const std::string& text, len_t s, len_t p, const lpf_t& lpf, sdsl::bit_vector& marker, const report_t& report) {
	for(len_t j = 1; j < p; ++j) {
		if(text[s+j-1] != text[s+j+2*p-1]) { break; }
		if(lpf[s+j] >= 2*p) { continue; } 	// actually: break! TODO, this will cause slowdown!
		if(!marker[s+j]) {
			marker[s+j]=true;
			report(s+j,p);
		}
	}
}

/**
 * @see linear_rightrotate
 * @param rmqlpf an RMQ data structure on the LPF array
 * @lcpq a function that returns the length of the longest common prefix of two suffixes of T
 */
template<class lpf_t, class report_t, class lcpq_t>
void rmq_rightrotate(len_t s, len_t p, const lpf_t& lpf, sdsl::bit_vector& marker, const report_t& report,
		const sdsl::rmq_succinct_sct<>& rmqlpf, const lcpq_t& lcpq) {
				const len_t replength = lcpq(s,s+p,2*p);
				const len_t right_range = std::min(s+p-1, s+replength-p);
				if(s+replength >= p && right_range > s) {
					std::stack<std::pair<len_t,len_t>> staple;
					staple.emplace(s+1, right_range);
					while(!staple.empty()) {
						const std::pair<len_t,len_t> range = staple.top(); staple.pop();
						DCHECK_LE(range.first, range.second);
						const len_t minlpf_idx =  rmqlpf(range.first, range.second);
						DCHECK_GE(minlpf_idx, range.first);
						DCHECK_LE(minlpf_idx, range.second);
						if(lpf[minlpf_idx] >= 2*p) { continue; }
						if(!marker[minlpf_idx]) {
							marker[minlpf_idx]=true;
							report(minlpf_idx,p);
						}
						if(range.first < minlpf_idx)
							staple.emplace(range.first, minlpf_idx-1);
						if(range.second > minlpf_idx)
							staple.emplace(minlpf_idx+1, range.second);
					}
				}
			}


/**
 * Constructs the suffix array
 * @param SA the suffix array
 * @param T the input text
 * @param n the length of the input text including its null-byte
 * @author Yuta Mori, "divsufsort" '15
 */
inline int32_t suffix_sort(const uint8_t* T, uint32_t* SA, int32_t n) {
	return divsufsort(T,(int32_t*)SA,n);
}
inline int32_t suffix_sort(const uint8_t* T, uint64_t* SA, int64_t n) {
	return divsufsort64(T,(int64_t*)SA,n);
}

// #if defined(DEBUG) && !defined(NDEBUG) //functions that cost more than constant time to check
// template<class T>
// void assert_permutation(const T& p, size_t n) {
//     for(size_t i = 0; i < n; ++i)
//     for(size_t j = 0; j < n; ++j) {
//         if(i == j) continue;
//         DCHECK_NE(p[i],p[j]); // << "at positions " << i << " and " << j;
//         DCHECK_LT(p[i],n);
//     }
// }
// #else
// template<class T> inline void assert_permutation(const T&, size_t) {}
// #endif

/**
 * Constructs the LPF array
 * @param sa the suffix array
 * @param lcp the LCP array
 * @return the LPF array
 * @author Crochemore et al., "LPF computation revisited", IWOCA'09
 */
template<class lpf_t, class lcp_t, class isa_t>
lpf_t gen_lpf(const lcp_t& lcp, const isa_t& isa) {
	const len_t n = lcp.size();
	lpf_t lcpcopy(n+1);
	for(len_t i = 0; i < n; ++i) { lcpcopy[i] = lcp[i]; }
	lcpcopy[n] = 0;
	lpf_t prev(n);
	lpf_t next(n);
	lpf_t lpf(n);
	constexpr len_t undef { static_cast<len_t>(-1) };
	for(len_t r = 0; r < n; ++r) {
		prev[r] = r-1;
		next[r] = r+1;
	}
	DCHECK_EQ(prev[0], undef);
	for(len_t j = n; j > 0; --j) {
		const len_t i = j-1; DCHECK_GT(j,0);
		const len_t r = isa[i];
		lpf[i] = std::max(lcpcopy[r], lcpcopy[next[r]] );
		DCHECK_LT(next[r],n+1);

		lcpcopy[next[r]] = std::min(lcpcopy[r], lcpcopy[next[r]]);
		if(prev[r] != undef) { next[prev[r]] = next[r]; }
		if(next[r] < n) { prev[next[r]] = prev[r]; }
	}
	return lpf;
}


/**
 * Constructs the phi array with phi[sa[i]] = sa[i-1]
 * @param sa the suffix array
 * @return the phi array
 */
template<class phi_t, class sa_t>
	inline phi_t construct_phi_array(const sa_t& sa) {
		const len_t n = sa.size();
		//assert_permutation(sa,n);

		phi_t phi(n);
		for(size_t i = 1, prev = sa[0]; i < n; i++) {
			phi[sa[i]] = prev;
			prev = sa[i];
		}
		phi[sa[0]] = sa[n-1];
		//assert_permutation(phi,n);
		return phi;
	}
/**
 * Constructs the PLCP array
 * @param phi the phi-array. Will be overwritten with PLCP
 * @text the input text
 * @author Kärkkäinen et. al, "Permuted Longest-Common-Prefix Array", CPM'09
 */
template<typename text_t, typename phi_t>
	inline void phi_algorithm(phi_t& phi, const text_t& text) {
		const size_t n = phi.size();
//			phi_t plcp(std::move(phi.data()));
		for(len_t i = 0, l = 0; i < n - 1; ++i) {
			const len_t phii = phi[i];
			DCHECK_LT(i+l, n);
			DCHECK_LT(phii+l, n);
			DCHECK_NE(i, phii);
			while(text[i+l] == text[phii+l]) { 
				l++;
				DCHECK_LT(i+l, n);
				DCHECK_LT(phii+l, n);
			}
			phi[i] = l;
			if(l) {
				--l;
			}
		}
	}

/**
 * Constructs the LCP array
 * @param PLCP the plcp array
 * @sa the suffix array
 */
    template <typename lcp_t, typename plcp_t, typename sa_t>
	inline lcp_t construct_lcp_array(const plcp_t& plcp, const sa_t& sa) {
        const len_t n = sa.size();
		lcp_t lcp(n);
		lcp[0] = 0;
		for(len_t i = 1; i < n; i++) { //TODO: start at 0, see line 149
			DCHECK_LT(sa[i], n);
			lcp[i] = plcp[sa[i]];
		}
		DCHECK_EQ( ([&lcp,&plcp,&sa] () {
				for(size_t i = 1; i < lcp.size(); ++i) { //TODO: start at 0, see line 149
				DCHECK_EQ(lcp[i], plcp[sa[i]]);
				}
				return true;
				}()), true);
		return lcp;
	}

/**
 * Compute the longst common prefix of T[a..] and T[b..]
 * @param isa the inverse suffix array
 * @param lcp the LCP array
 * @param rmqlcp an RMQ data structure on LCP
 */
template<class isa_t, class lcp_t>
inline len_t lcp_rmq(const isa_t& isa, const lcp_t& lcp, const sdsl::rmq_succinct_sct<>& rmqlcp, const len_t a, const len_t b) {
	DCHECK_NE(a,b);
	const len_t& ia = isa[a];
	const len_t& ib = isa[b];
	DCHECK_NE(ia,ib);
	const len_t idx = rmqlcp(std::min(ia,ib)+1, std::max(ia,ib));
	return lcp[idx];
}

inline len_t lcs_naive(const std::string& text, const len_t a, const len_t b, len_t upper_bound) {
	upper_bound = std::min(std::min(a,b),upper_bound);
	len_t i=0;
	while(upper_bound >= i && text[a-i] == text[b-i]) { ++i;}
	return i;
}

inline len_t lcp_naive(const std::string& text, const len_t a, const len_t b, const len_t upper_bound) {
	len_t i=0;
	while(upper_bound >= i && text[a+i] == text[b+i]) { ++i;}
	return i;
}

/**
 * Compute the longst common suffix of T[..a] and T[..b]
 * @param isai the inverse suffix array of the reversed text
 * @param lcs the LCP array of the reversed text
 * @param rmqlcs an RMQ data structure on lcs
 */
template<class isa_t, class lcp_t>
inline len_t lcs_rmq(const isa_t& isai, const lcp_t& lcs, const sdsl::rmq_succinct_sct<>& rmqlcs, const len_t a, const len_t b) {
	const len_t& n = isai.size();
	DCHECK_NE(a,b);
	DCHECK_LE(a, n-2);
	DCHECK_LE(b, n-2);
	const len_t& ia = isai[n-2-a];
	const len_t& ib = isai[n-2-b];
	DCHECK_NE(ia,ib);
	const len_t idx = rmqlcs(std::min(ia,ib)+1, std::max(ia,ib));
	return lcs[idx];
};


/**
 * @param text the input text. We revert the text in place (and back again), so it cannot be made constant
 * @param report_square a function of type void (*report_square)(len_t start, len_t period)) 
 */
template<class report_t>
void compute_distinct_squares(std::string& text, const report_t& report_square) {
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



	auto lcpq = [&text, &isa,&lcp,&rmqlcp] (const len_t a, const len_t b, const len_t upper_bound) {
		if(upper_bound < RMQ_THRES_LCP) {
			return lcp_naive(text, a,b, upper_bound);
		} else {
			const len_t ret = lcp_naive(text, a,b, RMQ_THRES_TRYLCP);
			if(ret <= RMQ_THRES_TRYLCP) return ret;
		}
		return lcp_rmq(isa,lcp,rmqlcp,a,b);
	};
	auto lcsq = [&text, &n,&isai,&lcs,&rmqlcs] (const len_t a, const len_t b) {
		const len_t ret = lcs_naive(text, a,b, RMQ_THRES_TRYLCP);
		if(ret <= RMQ_THRES_TRYLCP) return ret;
		return lcs_rmq(isai,lcs,rmqlcs,a,b);
	};



	DODEBUG(std::set<std::string> check;)

	auto report = [&] (const len_t pos, const len_t period) {
		report_square(pos,period);
		if(VLOG_IS_ON(1)) {
			DVLOG(1) << "T[" << pos << "," << (pos+period*2-1) << "] = " << text.substr(pos,period) << "," << text.substr(pos+period,period) << " | ";
			DVLOG(1) << "lpf=" << lpf[pos] << " ";

			std::stringstream ss;
			for(len_t i = pos; i < pos+period*2; ++i) { 
				if(i == pos+period) DVLOG(1) << ","; 
				ss << "(" << ((size_t)text[i]) << ")"; 
			}
			DVLOG(1) << ss;
		}
		//checks
		DCHECK_EQ([&] () {
				for(len_t i = 0; i < period; ++i) { 
				DCHECK_EQ(text[pos+i],text[pos+period+i]);
				}
				DODEBUG(
				const std::string& sub = text.substr(pos,period*2);
				DCHECK(check.find(sub) == check.end());
				check.insert(sub);
				);


				return true;
		}(),true);
	};


	DVLOG(1) << "LZ-Fact";
	len_t maxperiod = 0;
	len_t num_factors = 1;
	for(len_t i = 1; i < text.length(); ) {
		const len_t factor_length = std::max(static_cast<len_t>(1),lpf[i]);
		const len_t next_factor_begin = i+factor_length;
		const len_t next_factor_length = std::max(static_cast<len_t>(1),lpf[next_factor_begin]);
		maxperiod = std::max(next_factor_length+factor_length, maxperiod);
		i = next_factor_begin;
		++num_factors;
		DVLOG(1) << i << ", ";
	}
	DVLOG(1) << "Longest period: " << maxperiod;

	len_t* Z = new len_t[num_factors]; //factor_id of the next long factor
	len_t* P = new len_t[num_factors]; //starting position of the next long factor
	memset(Z, 0, sizeof(len_t)*num_factors);

	for(len_t p = 1; p < maxperiod; ++p) {
		sdsl::bit_vector marker(n);

		auto report_and_rotate  = [&] (const len_t s, const len_t p, char type) {
			if(lpf[s] < 2*p && !marker[s]) { 
				DVLOG(1) << type;
				report(s,p); 
				marker[s]=true;
			}
			DCHECK([&] () { // test if RMQ and linear return the same states
				std::vector<len_t> pos;
				sdsl::bit_vector marker1(n);
				linear_rightrotate(text,s,p,lpf,marker1, [&pos] (len_t a, len_t) { pos.push_back(a); });
				std::sort(pos.begin(), pos.end());
				std::vector<len_t> pos2;
				sdsl::bit_vector marker2(n);
				rmq_rightrotate(s,p,lpf,marker2, [&pos,&pos2] (len_t a, len_t) { 
						pos2.push_back(a);
						DCHECK(std::find(pos.begin(), pos.end(), a) != pos.end());
						}, rmqlpf, lcpq);
				std::sort(pos2.begin(), pos2.end());
				DCHECK(pos == pos2);
				DCHECK(marker1 == marker2);
				return (pos == pos2) && (marker1 == marker2);
				}
				());

			if(p < RMQ_THRES_LPF) {
				linear_rightrotate(text,s,p,lpf,marker, report);
			} else {
				rmq_rightrotate(s,p,lpf,marker, report, rmqlpf, lcpq);
			}
		};

		for(len_t i = 0, factor_id = 0; i < text.length();) { DCHECK_LT(factor_id, num_factors);
			if(Z[factor_id] != 0) {
				const len_t old_factor_id = factor_id;
				do {
					i = P[factor_id];
					factor_id = Z[factor_id];
					if(i >= text.length()) break;
				}
				while (Z[factor_id] != 0);
				Z[old_factor_id] = factor_id; DCHECK_LE(factor_id, num_factors);
				P[old_factor_id] = i;
				if(i >= text.length()) break;
			}

			len_t factor_length = std::max(static_cast<len_t>(1),lpf[i]);
			len_t next_factor_begin = i+factor_length;
			len_t next_factor_length = std::max(static_cast<len_t>(1),lpf[next_factor_begin]);
			
			if(factor_length+next_factor_length < p) {
				const len_t old_factor_id = factor_id;
				while(factor_length+next_factor_length < p) {
					if(Z[factor_id] != 0) {
						do {
						i = P[factor_id];
						factor_id = Z[factor_id];
						if(i >= text.length()) break;
						} while(Z[factor_id] != 0);
					} else {
						i+=factor_length;
						++factor_id;
					}
					if(i >= text.length()) break;
					factor_length = std::max(static_cast<len_t>(1),lpf[i]);
					next_factor_begin = i+factor_length;
					next_factor_length = std::max(static_cast<len_t>(1),lpf[next_factor_begin]);
				}
				Z[old_factor_id] = factor_id; DCHECK_LE(factor_id, num_factors);
				P[old_factor_id] = i;
				if(i >= text.length()) break;
			}

			//backward
			if(factor_length >= p) {
				const len_t q = next_factor_begin-p; DCHECK_GE(next_factor_begin,p);
				const len_t length_r = lcpq(next_factor_begin,q,p);
				if(length_r > 0 ) {
					const len_t length_l = (q == 0) ? 0 : lcsq(next_factor_begin-1,q-1);
					if(length_l + length_r >= p) {
						const len_t s = std::max(q - length_l, (q +1< p) ? 0 : q - p + 1); DCHECK_GE(q, length_l);
						if(!marker[s]) {
							report_and_rotate(s,p,'A');
						}
					}
				}
			}
			// forward
			if(factor_length+next_factor_length >= p && i+p < n) {
				const len_t q = i+p;
				const len_t length_l = (i == 0) ? 0 : lcsq(i-1,q-1); DCHECK_GE(i, length_l); 
				const len_t s = std::max(i - length_l, (i < p+1) ? 0 : i - p + 1); 
				if(length_l > 0 && s+p <= next_factor_begin) {
					const len_t length_r = lcpq(i,q,p);
					if(length_r > 0 && length_l + length_r >= p) { 
						if(!marker[s]) {
							report_and_rotate(s,p,'B');
						}
					}
				}
			
			}

			i+=factor_length;
			++factor_id;
		}
	}
	delete [] Z;
	delete [] P;

}



