#ifndef MY_COMPLEX_DEF_H_
#define MY_COMPLEX_DEF_H_

#include<complex>

//

typedef std::complex<float> scomplex;
typedef std::complex<double> dcomplex;

//

template<typename T>
struct complex_tool {
	typedef T base_type;
	static inline T conj(const T& x) { return x; }
	static inline bool is_complex() { return false; }
	static inline T real(const T& x) { return x; }
	static inline T imag(const T& x) { return 0; }
	static inline T norm(const T& x) { return x*x; }
	static inline T make(const T& re, const T& im) { return re; }
};

template<>
struct complex_tool<scomplex> {
	typedef float base_type;
	static inline scomplex conj(const scomplex& x) { return std::conj<float>(x); }
	static inline bool is_complex() { return true; }
	static inline float real(const scomplex& x) { return x.real(); }
	static inline float imag(const scomplex& x) { return x.imag(); }
	static inline float norm(const scomplex& x) { return std::norm(x); }
	static inline scomplex make(const float& re, const float& im) { return scomplex(re,im); }
};

template<>
struct complex_tool<dcomplex> {
	typedef double base_type;
	static inline dcomplex conj(const dcomplex& x) { return std::conj<double>(x); }
	static inline bool is_complex() { return true; }
	static inline double real(const dcomplex& x) { return x.real(); }
	static inline double imag(const dcomplex& x) { return x.imag(); }
	static inline double norm(const dcomplex& x) { return std::norm(x); }
	static inline dcomplex make(const double& re, const double& im) { return dcomplex(re,im); }
};

template<typename T>
T dirac_bracket(const int n, const T* x, const T* y) {
	T sum=0;
	for (int i=0;i<n;i++) sum+=complex_tool<T>::conj(x[i])*y[i];
	return sum;
}

template<typename OP, typename T>
T dirac_bracket(const OP& op, const int& n, const T* v1, const T* v2) 
{
	// calculate <v1|op|v2>
	T sum=0;
	for (int i=0;i<n;++i) {
		int j=i;
		int p=op.apply(j);
		if (p!=0) sum+=complex_tool<T>::conj(v1[j])*v2[i]*static_cast<T>(p);
	}
	return sum;
}

template<typename OP, typename BASIS1, typename BASIS2, typename T>
T dirac_bracket(const OP& op, const BASIS1& bL, const T* vL, 
		const BASIS2& bR, const T* vR) 
{
	T sum=0;
	if (bL.size()<bR.size()) {
		// (<R|op^+|L>)^*
		OP opcj=op.conj();
		for (int i=0;i<bL.size();++i) {
			int s=bL[i];
			int p=opcj.apply(s); 
			if (p==0) continue;
			int j=bR.lookup(s);
			if (j==-1) continue;
			// std::cout << i << "\t" << j << std::endl;
			sum+=complex_tool<T>::conj(vR[j])*vL[i]*static_cast<T>(p);
		}
		return complex_tool<T>::conj(sum);
	}
	for (int i=0;i<bR.size();++i) {
		int s=bR[i];
		int p=op.apply(s); 
		if (p==0) continue;
		int j=bL.lookup(s);
		if (j==-1) continue;
		sum+=complex_tool<T>::conj(vL[j])*vR[i]*static_cast<T>(p);
	}
	return sum;
}

template<typename OP, typename BASIS1, typename BASIS2>
int dirac_bracket_TF(const OP& op, const BASIS1& bL, const BASIS2& bR) 
{
	if (bL.size()<bR.size()) {
		// (<R|op^+|L>)^*
		OP opcj=op.conj();
		for (int i=0;i<bL.size();++i) {
			int s=bL[i];
			int p=opcj.apply(s); 
			if (p==0) continue;
			int j=bR.lookup(s);
			if (j==-1) continue;
			return 1;
		}
		return 0;
	}
	// <L|op|R>
	for (int i=0;i<bR.size();++i) {
		int s=bR[i];
		int p=op.apply(s); 
		if (p==0) continue;
		int j=bL.lookup(s);
		if (j==-1) continue;
		return 1;
	}
	return 0;
}

#endif
