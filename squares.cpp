#include <sdsl/suffix_trees.hpp>
#include <iostream>
#include "/scripts/code/vector.hpp"
#include <string>
#include <fstream>
#include <streambuf>
#include <divsufsort.h>
#include <divsufsort64.h>

//#define DEBUG
#ifdef DEBUG
#define tdc_hdebug(x) x
#else
#define tdc_hdebug(x)
#endif

using namespace std;
using namespace sdsl;

typedef cst_sct3<> cst_t;
typedef cst_t::char_type char_type;

typedef uint32_t len_t;


inline int32_t suffix_sort(const uint8_t* T, uint32_t* SA, int32_t n) {
	return divsufsort(T,(int32_t*)SA,n);
}
inline int32_t suffix_sort(const uint8_t* T, uint64_t* SA, int64_t n) {
	return divsufsort64(T,(int64_t*)SA,n);
}

#if defined(DEBUG) && !defined(NDEBUG) //functions that cost more than constant time to check
template<class T>
void assert_permutation(const T& p, size_t n) {
    for(size_t i = 0; i < n; ++i)
    for(size_t j = 0; j < n; ++j) {
        if(i == j) continue;
        DCHECK_NE(p[i],p[j]); // << "at positions " << i << " and " << j;
        DCHECK_LT(p[i],n);
    }
}
#else
template<class T> inline void assert_permutation(const T&, size_t) {}
#endif

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


template<class phi_t, class sa_t>
	inline phi_t construct_phi_array(const sa_t& sa) {
		const len_t n = sa.size();
		assert_permutation(sa,n);

		phi_t phi(n);
		for(size_t i = 1, prev = sa[0]; i < n; i++) {
			phi[sa[i]] = prev;
			prev = sa[i];
		}
		phi[sa[0]] = sa[n-1];
		assert_permutation(phi,n);
		return phi;
	}
/*
 * Constructs PLCP inplace in phi
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

    template <typename lcp_t, typename plcp_t, typename sa_t>
	inline lcp_t construct_lcp_array(const plcp_t& plcp, const sa_t& sa) {
        const len_t n = sa.size();
		//DCHECK_EQ(plcp[sa[0]],0);
		//plcp[sa[0]] = 0;
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

void check_toolchain(const std::string& text) {
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
	len_t lcs_naive(const std::string& text, const len_t a, const len_t b) {
		const len_t m = std::min(a,b);
		len_t i=0;
		while(m >= i && text[a-i] == text[b-i]) { ++i;}
		return i;
	};

	len_t lcp_naive(const std::string& text, const len_t a, const len_t b) {
		len_t i=0;
		while(text[a+i] == text[b+i]) { ++i;}
		return i;
	};


struct Square {
	len_t start;
	len_t period;
	Square(len_t _start, len_t _period) : start(_start), period(_period) {}
} __attribute__((__packed__));



Vector<Square> find_naively(const std::string& text) {
	const len_t n = text.length()+1;
	std::set<std::string> check;
	Vector<Square> squares;

	auto report = [&] (const len_t pos, const len_t period) {
		const std::string sub = text.substr(pos,period*2);
		if(check.find(sub) != check.end()) return;
		check.insert(sub);

		cout << "T[" << pos << "," << (pos+period*2-1) << "] = " << text.substr(pos,period) << "," << text.substr(pos+period,period) << " | ";
		cout << endl;
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


int main(int argc, char* argv[])
{
    if (argc < 2) {
        cout << "Usage: "<< argv[0] << " file" << endl;
        return 1;
    }


	std::ifstream t(argv[1]);
	std::string text((std::istreambuf_iterator<char>(t)),
			std::istreambuf_iterator<char>());

	Vector<Square> naive_squares { find_naively(text) };
//	return 1;


	tdc_hdebug(check_toolchain(text));


	const len_t n = text.length()+1;
	Vector<len_t> sa(n);
	suffix_sort((const uint8_t*)text.c_str(), sa.data(), n);

	Vector<len_t> plcp { construct_phi_array<Vector<len_t>,decltype(sa)>(sa) };
	phi_algorithm(plcp,text);
	const Vector<len_t> lcp { construct_lcp_array<Vector<len_t>,decltype(plcp),decltype(sa)>(plcp,sa) };
	Vector<len_t>().swap(plcp); // delete plcp
	const rmq_succinct_sct<> rmqlcp { &lcp };

	Vector<len_t> isa(n);
	for(size_t i = 0; i < n; ++i) {
		DCHECK_LT(sa[i], n);
		isa[sa[i]] = i;
	}

	const Vector<len_t> lpf = gen_lpf<Vector<len_t>,decltype(lcp),decltype(isa)>(lcp, isa);
	tdc_hdebug(check_lpf(text, lpf));
	const rmq_succinct_sct<> rmqlpf { &lpf };

	//generate LCP^-1
	std::reverse(text.begin(), text.end());
	suffix_sort((const uint8_t*)text.c_str(), sa.data(), n);

	plcp = construct_phi_array<Vector<len_t>,decltype(sa)>(sa);
	phi_algorithm(plcp,text);
	Vector<len_t> lcs = construct_lcp_array<Vector<len_t>,decltype(plcp),decltype(sa)>(plcp,sa);
	Vector<len_t>().swap(plcp); // delete plcp
	const rmq_succinct_sct<> rmqlcs { &lcs };

	Vector<len_t> isai(n);
	for(size_t i = 0; i < n; ++i) {
		DCHECK_LT(sa[i], n);
		isai[sa[i]] = i;
	}
	Vector<len_t>().swap(sa); //delete sa
	std::reverse(text.begin(), text.end());



	auto lcpq = [&isa,&lcp,&rmqlcp] (const std::string&, const len_t a, const len_t b) {
		DCHECK_NE(a,b);
		const len_t& ia = isa[a];
		const len_t& ib = isa[b];
		DCHECK_NE(ia,ib);
		const len_t idx = rmqlcp(std::min(ia,ib)+1, std::max(ia,ib));
		return lcp[idx];
	};
	auto lcsq = [&n,&isai,&lcs,&rmqlcs] (const std::string&, const len_t a, const len_t b) {
		DCHECK_NE(a,b);
		DCHECK_LE(a, n-2);
		DCHECK_LE(b, n-2);
		const len_t& ia = isai[n-2-a];
		const len_t& ib = isai[n-2-b];
		DCHECK_NE(ia,ib);
		const len_t idx = rmqlcs(std::min(ia,ib)+1, std::max(ia,ib));
		return lcs[idx];
	};

//	tdc_hdebug(
			std::set<std::string> check;
//			)

	Vector<Square> squares;
	auto report = [&] (const len_t pos, const len_t period) {
		squares.emplace_back(pos,period);

		cout << "T[" << pos << "," << (pos+period*2-1) << "] = " << text.substr(pos,period) << "," << text.substr(pos+period,period) << " | ";
		cout << "lpf=" << lpf[pos] << " ";

		for(len_t i = pos; i < pos+period*2; ++i) { 
			if(i == pos+period) cout << ","; 
			cout << "(" << ((size_t)text[i]) << ")"; 
		}
		cout <<	endl;

		//checks
		DCHECK_EQ([&] () {
				for(len_t i = 0; i < period; ++i) { 
				DCHECK_EQ(text[pos+i],text[pos+period+i]);
				}
//				tdc_hdebug(
				const std::string& sub = text.substr(pos,period*2);
				DCHECK(check.find(sub) == check.end());
				check.insert(sub);
//				);


				return true;
		}(),true);
	};

	tdc_hdebug(
	for(len_t i = 0; i < text.length();++i) {
		for(len_t j = 0; j < text.length();++j) {
			if(j==i) continue;
			DCHECK_EQ(lcp_naive(text,i,j), lcpq(text,i,j));
			DCHECK_EQ(lcs_naive(text,i,j), lcsq(text,i,j));
		}
	} );

	std::cout << "LZ-Fact" << std::endl;
	len_t maxperiod = 0;
	for(len_t i = 1; i < text.length(); ) {
		const len_t factor_length = std::max(static_cast<len_t>(1),lpf[i]);
		const len_t next_factor_begin = i+factor_length;
		const len_t next_factor_length = std::max(static_cast<len_t>(1),lpf[next_factor_begin]);
		maxperiod = std::max(next_factor_length+factor_length, maxperiod);
		i = next_factor_begin;
		std::cout << i << ", ";
	}
	std::cout << "Longest period: " << maxperiod << std::endl; 


	for(len_t p = 1; p < maxperiod; ++p) {
		sdsl::bit_vector marker(n);

		auto report_and_rotate  = [&] (const len_t s, const len_t p, char type) {
			if(lpf[s] < 2*p && !marker[s]) { 
				cout << type;
				report(s,p); 
				marker[s]=true;
			}

			{//test
				std::vector<len_t> pos;
				for(len_t j = 1; j < p; ++j) {
					if(text[s+j-1] != text[s+j+2*p-1]) { break; }
					if(lpf[s+j] >= 2*p) { continue; } 	// actually: break! TODO, this will cause slowdown!
					if(!marker[s+j]) {
						pos.push_back(s+j);
					}
				}
				std::sort(pos.begin(), pos.end());
				std::vector<len_t> pos2;
				const len_t replength = lcpq(text,s,s+p);
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
							pos2.push_back(minlpf_idx);
							if(std::find(pos.begin(), pos.end(), minlpf_idx) == pos.end()) { DCHECK(false); }
						}
						if(range.first < minlpf_idx)
							staple.emplace(range.first, minlpf_idx-1);
						if(range.second > minlpf_idx)
							staple.emplace(minlpf_idx+1, range.second);
					}
				}
					std::sort(pos2.begin(), pos2.end());
					DCHECK(pos == pos2);


				//end test
			}
			for(len_t j = 1; j < p; ++j) {
				if(text[s+j-1] != text[s+j+2*p-1]) { break; }
				if(lpf[s+j] >= 2*p) { continue; } 	// actually: break! TODO, this will cause slowdown!
				if(!marker[s+j]) {
					cout << "R";
					report(s+j,p);
					marker[s+j]=true;
				}
			}
		};


		for(len_t i = 0; i < text.length();) {
			const len_t factor_length = std::max(static_cast<len_t>(1),lpf[i]);
			const len_t next_factor_begin = i+factor_length;
			//backward
			if(factor_length >= p) {
				const len_t q = next_factor_begin-p; DCHECK_GE(next_factor_begin,p);
				const len_t length_r = lcpq(text,next_factor_begin,q);
				const len_t length_l = (q == 0) ? 0 : lcsq(text,next_factor_begin-1,q-1);
				if(length_l + length_r >= p && length_r > 0) {
					const len_t s = std::max(q - length_l, (q +1< p) ? 0 : q - p + 1); DCHECK_GE(q, length_l);
					if(!marker[s])
						report_and_rotate(s,p,'A');
				}
			}
			// forward
			const len_t next_factor_length = std::max(static_cast<len_t>(1),lpf[next_factor_begin]);
			if(factor_length+next_factor_length >= p && i+p < n) {
				const len_t q = i+p;
				const len_t length_r = lcpq(text,i,q);
				const len_t length_l = (i == 0) ? 0 : lcsq(text,i-1,q-1); DCHECK_GE(i, length_l); 
				const len_t s = std::max(i - length_l, (i < p+1) ? 0 : i - p + 1); 
				if(length_r > 0 && length_l > 0 && length_l + length_r >= p && s+p <= next_factor_begin) { // TODO: actually s+p < next_factor_begin
					if(!marker[s])
						report_and_rotate(s,p,'B');
					//bordercase: if a factor divides a square exactly in the middle
					// if(i - length_l < s && s == i - p + 1) report_and_rotate(i-p,p,'C');
				}
			
			}

			i+=factor_length;
		}
	}


	std::sort(squares.begin(), squares.end(), [] (const Square& a, const Square& b) 
			{
			if(a.start == b.start) return a.period < b.period;
			return a.start < b.start; 
			});

	for(size_t s = 0; s < std::min(naive_squares.size(),squares.size()); ++s) {
		if(squares[s].start != naive_squares[s].start || squares[s].period != naive_squares[s].period) {
			std::cout << "Naive : " << naive_squares[s].start << "," << naive_squares[s].period << endl;
			std::cout << "Own   : " << squares[s].start << "," << squares[s].period << endl;
		}
		// DCHECK_EQ(squares[s].start, naive_squares[s].start);
		// DCHECK_EQ(squares[s].period, naive_squares[s].period);
	}
	DCHECK_EQ(squares.size(), naive_squares.size());



}
