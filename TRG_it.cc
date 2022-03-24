// Last modified date : 2019-12-20

#include "itensor/all.h"
#include <iostream>
#include <cstring>
#include <vector>
#include <list>
#include <sys/time.h>
#include <cmath>
#include <fstream>
#include <complex_bessel.h>
#include <mkl_lapack.h>
#include <mkl_blas.h>
#include <mkl_lapacke.h>
#define MKL_Complex16 std::complex<double>
#define MKL_Complex8 std::complex<float>
#define _MKL_Complex16 std::complex<double>


typedef std::complex<double> dcomplex;

using namespace std;
using namespace itensor;
/*
void zheev( const char* jobz, const char* uplo, const MKL_INT* n,
            MKL_Complex16* a, const MKL_INT* lda, double* w,
            MKL_Complex16* work, const MKL_INT* lwork, double* rwork,
            MKL_INT* info );
*/

complex<double> Z_HOTRG(complex<double> K, complex<double> h, int dcut, int L);

//void heev_f(const char jobz, const char uplo, int n, dcomplex* a, int lda, double* w);

int main(int argc, char* argv[])
{
	int L = 3;      //L : size 2^L x 2^L
	int Lp = pow(4,L);      // number of sites
	int dcut = 110;  //dcut : Bond dimension cutoff

	cout.precision(14);

	Cplx K (1.1199 + 0.00*atoi(argv[1]),0);

	double hr = 0;
	double hi = atof(argv[2]);

	Cplx h0 (hr, 0.001);
	Cplx h (hr, hi);

	Cplx o0 = Z_HOTRG(K, h0, dcut, L)*Lp;
	Cplx o1 = Z_HOTRG(K, h, dcut, L)*Lp;

	Cplx Z = exp(o1-o0);

	if(imag(Z)>=0)
		printf("%.16g+%.16gj\n",real(Z),imag(Z));
	else
		printf("%.16g%.16gj\n",real(Z),imag(Z));


        return 0;
}


complex<double> Z_HOTRG(complex<double> K, complex<double> h, int dcut, int L)
{
	Cplx B = K;
	vector<int> in1, in2;
	double pi = acos(-1.0);

	complex<double> Ren[2*L];

	int ind = 0;


	vector<double> MM1, MM2, Us1, Us2;
	vector<double> Ss1, Ss2;
	auto d0 = dcut-1;		//Bond dimension
	complex<double> *bes = new complex<double>[d0];
	// q : q-states. Large q ~ XY model. 
	// For our available bond dimensions, q = 1000 is effectively same as the XY limit. 
	auto q = 1000;
	auto dx = d0, dy = d0;
	auto x1 = Index("x1",d0,Xtype);
	auto y1 = Index("y1",d0,Ytype);
	auto x2 = prime(x1,1);
	auto y2 = prime(y1,1);
	auto x1p = Index("x1p",dx);
	auto x2p = prime(x1p,1);
	auto x3 = x1;
	auto x4 = x1;
	auto y3 = x1;
	auto y4 = x1;
	auto i = x1p;
	auto c1 = combiner(i,x1p);
	auto c2 = c1;
	auto j = i;
	auto x = i;
	auto xp = i;
	auto y = i;
	auto yp = i;
	auto ui = i;
	auto uii = i;
	auto y1p = i;
	auto y2p = i;
	auto xx1 = i;


	// Defining the local tensor for XY model. 
	for(auto s1: range(d0))
	{
		bes[s1] = sqrt(sp_bessel::besselI(0.5*(d0-1)-(double)s1,B));
	}


	auto T = ITensor(x1,x2,y1,y2);
	auto T1 = T;
	auto A = T;
	auto UU1 = T;

	Cplx P;
	for(auto s1 : range1(d0))
	for(auto s2 : range1(d0))
	for(auto s3 : range1(d0))
	for(auto s4 : range1(d0))
	{
		P = bes[s1-1]*bes[s2-1]*bes[s3-1]*bes[s4-1]*sp_bessel::besselI(double(s2+s4-s1-s3),B*h);
		T.set(x1(s1),x2(s2),y1(s3),y2(s4),P);
	}

	delete [] bes;

	auto M = ITensor(x,xp,y1,y2);

	double s1 = 0, s2 = 0;

	double dxn, dyn;
	
	// End of definition of local tensor


	auto Tn = ITensor(x1,y1,x2,y2);

	for(auto iter : range1(L))		// RG iteration starts here. 
	{
		dxn = (dcut > dx*dx) ? dx*dx : dcut;
		x1 = noprime(findtype(T,Xtype));
		x2 = prime(x1,1);
		y1 = noprime(findtype(T,Ytype));
		y2 = prime(y1,1);
		x3 = Index("x3",dx);
		x4 = Index("x4",dx);

		T1 = T;


		x1p = Index("x1p",dx);
		x2p = prime(x1p,1);

		i = Index("i",dy);
		j = Index("j",dy);

		T1 = T1*delta(x2,x2p);
		T1 = T1*delta(x1,x1p);

		T *=delta(y2,i);	T1 *=delta(y1,i);	T1 *=delta(y2,y2);


		// Given y-directional contraction of T and T1, 
		// eigenvalues of MT1 (or MT2) represents 
		// the higher-order singular values of (T*T1). 
		// 
		// Graphical representation of T*T1 is as below.
		//
		//		|
		//		|
		//        ------T-------
		//              |
		//              |
		//        ------T1------
		//              |
		//              |
		//

		ITensor MT1 = T;
		ITensor temp = T1;
		MT1 = (MT1*delta(x1,x3)*delta(i,j));
		temp = (temp*delta(x1p,x4)*delta(i,j));
		MT1.conj();	temp.conj();

		MT1 = MT1*T;
		temp = temp*T1;
		MT1 = MT1*temp;
		MT1 = MT1.takeReal();


		ITensor MT2 = T;
		temp = T1;
		MT2 = (MT2*delta(x2,x3)*delta(i,j));
		temp = (temp*delta(x2p,x4)*delta(i,j));
		MT2.conj();	temp.conj();

		MT2= MT2*T;
		temp = temp*T1;


		MT2 = MT2*temp;
		MT2 = MT2.takeReal();


		MT2 = MT2*delta(x2,x1)*delta(x2p,x1p);


		MT2 = MT2*delta(x1,x2)*delta(x1p,x2p);
		ITensor UU;

		ui = Index("A",dx*dx);
		uii = ui;

		int sindex;

		MM1.resize(dx*dx*dx*dx);
		MM2.resize(dx*dx*dx*dx);
		Ss1.resize(dx*dx);
		Ss2.resize(dx*dx);
		in1.resize(dx*dx);
		in2.resize(dx*dx);

		x = Index("x",dx*dx,Atype);
		auto ii1 = Index("ii1",dx);
		auto ii2 = Index("ii2",dx);


		if(dxn == dcut)
		{
			
			for(auto i1 : range1(dx))
			for(auto i2 : range1(dx))
			for(auto i3 : range1(dy))
			for(auto i4 : range1(dy))
			{
				MM1[(i2-1+(i1-1)*dx)+(i4-1+(i3-1)*dx)*dx*dx] = MT1.real(x1(i1),x1p(i2),x3(i3),x4(i4));
				MM2[(i2-1+(i1-1)*dx)+(i4-1+(i3-1)*dx)*dx*dx] = MT2.real(x2(i1),x2p(i2),x3(i3),x4(i4));
			}
			MT1 = ITensor(0);
			MT2 = ITensor(0);
			LAPACKE_dsyev( LAPACK_COL_MAJOR,'V', 'U', dx*dx,MM1.data(),dx*dx,Ss1.data());
			LAPACKE_dsyev( LAPACK_COL_MAJOR,'V', 'U', dx*dx,MM2.data(),dx*dx,Ss2.data());
		
			
	
			for(int q=0;q<dx*dx;q++)
			{
				in1[q] = q; in2[q] = q;
			}

			for(int q=0;q<dx*dx-1;q++)
			{
				for(int p=0;p<dx*dx-1;p++)
				{
					if(fabs(Ss1[in1[p]])>fabs(Ss1[in1[p+1]]))
					{
						int tmp = in1[p];
						in1[p] = in1[p+1];
						in1[p+1] = tmp;
					}
				}
			}

			for(int q=0;q<dx*dx-q;q++)
			for(int p=0;p<dx*dx-1;p++)
			if(fabs(Ss2[in2[p]])>fabs(Ss2[in2[p+1]]))
			{
				int tmp = in2[p];
				in2[p] = in2[p+1];
				in2[p+1] = tmp;
			}


			for(auto i1 : range(dx*dx-dcut))
			{

				s1 += fabs(Ss1[in1[i1]]);
				s2 += fabs(Ss2[in2[i1]]);
			}

			if(s1<s2)
			{
				ui = Index("A",dcut);
				UU = ITensor(ii1,ii2,ui);
				for(auto i1 : range1(dx))
				for(auto i2 : range1(dx))
				for(auto i3 : range1(dcut))
					UU.set(ii1(i1),ii2(i2),ui(i3),MM1[(i2-1+(i1-1)*dx) + (in1[dx*dx-i3])*dx*dx]);
				vector<double>().swap(Ss1);
				vector<double>().swap(MM2);
				vector<double>().swap(Ss2);
				sindex = 1;
				vector<int>().swap(in2);
			}
			else
			{
//				cout << "y-direction : choose UR" << endl;
				ui = Index("A",dcut);
				UU = ITensor(ii1,ii2,ui);
				for(auto i1 : range1(dx))
				for(auto i2 : range1(dx))
				for(auto i3 : range1(dcut))
					UU.set(ii1(i1),ii2(i2),ui(i3),MM2[(i2-1+(i1-1)*dx) + (in2[dx*dx-i3])*dx*dx]);
				vector<double>().swap(MM1);
				vector<double>().swap(Ss1);
				vector<double>().swap(Ss2);
				sindex = 0;
				vector<int>().swap(in1);
			}
		}
		else
		{
			for(auto i1 : range1(dx))
			for(auto i2 : range1(dx))
			for(auto i3 : range1(dy))
			for(auto i4 : range1(dy))
			{
				MM1[(i2-1+(i1-1)*dx)+(i4-1+(i3-1)*dx)*dx*dx] = MT1.real(x1(i1),x1p(i2),x3(i3),x4(i4));
			}
			MT1 = ITensor(0);
			MT2 = ITensor(0);
			LAPACKE_dsyev( LAPACK_COL_MAJOR,'V', 'U', dx*dx,MM1.data(),dx*dx,Ss1.data());
			ui = Index("A",dx*dx);
			UU = ITensor(ii1,ii2,ui);
			for(auto i1 : range1(dx))
			for(auto i2 : range1(dx))
			for(auto i3 : range1(dx*dx))
			{
				UU.set(ii1(i1),ii2(i2),ui(i3),MM1[(i2-1+(i1-1)*dx) + (dx*dx-i3)*dx*dx]);
			}
			vector<double>().swap(MM1);
			vector<double>().swap(Ss1);
			vector<double>().swap(MM2);
			vector<double>().swap(Ss1);
		}

	
		x = Index("x",dxn,Atype);
		xp = Index("xp",dxn,Atype);
		auto Tk = ITensor(x,xp,y1,y2);
		if(dcut==dxn)
		{
			if(sindex==0)
			{
			auto j1 = Index("j1",dx);
			auto j2 = Index("j2",dx);
			UU1 = UU;
			UU1 = UU1 * delta(ii2,j2);
			UU1 = UU1 * delta(ii1,j1);
			UU1 *= delta(ui,xp);
			T = T*delta(x2, j1);	T = T*delta(x1,ii1);
			T1 = T1*delta(x2p,j2);	T1 = T1*delta(x1p,ii2);
			for(int s1=1; s1<=dxn; s1++)
			{

			auto Uk = ITensor(ii1,ii2);
			
			for(auto i1 : range1(dx))
			for(auto i2 : range1(dx))
				Uk.set(ii1(i1),ii2(i2),MM2[(i2-1+(i1-1)*dx) + in2[dx*dx-s1]*dx*dx]);


		
			Tn = Uk*T;

			Tn = Tn*T1;

			Tn = Tn*UU1.conj();
				for(auto i2: range1(int(dxn)))
				for(auto i3: range1(dy))
				for(auto i4: range1(dy))
				{
					Tk.set(x(s1),xp(i2),y1(i3),y2(i4),Tn.cplx(xp(i2),y1(i3),y2(i4)));
				}			
	
			}
			vector<double>().swap(MM2);
			vector<int>().swap(in2);
			Tn = Tk;
			Tk = ITensor(0);
			}
			else
			{
			UU1 = UU;
			auto j1 = Index("j1",dx);
			auto j2 = Index("j2",dx);
			UU1 = UU1 * delta(ii2,j2);
			UU1 = UU1 * delta(ii1,j1);
			UU1 *= delta(ui,xp);

			T = T*delta(x2,ii1); T = T*delta(x1,j1);
			T1 = T1*delta(x2p,ii2);  T1 = T1*delta(x1p,j2);


			for(int s1=1; s1<=dxn; s1++)
			{

			auto Uk = ITensor(ii1,ii2);
			for(auto i1 : range1(dx))
			for(auto i2 : range1(dx))
				Uk.set(ii1(i1),ii2(i2),MM1[(i2-1+(i1-1)*dx) + in1[dx*dx-s1]*dx*dx]);

 
			Tn = Uk*T;
 
			Tn = Tn*T1;
            
			Tn = Tn*UU1.conj();

				for(auto i2: range1(int(dxn)))
				for(auto i3: range1(dy))
				for(auto i4: range1(dy))
				{
					Tk.set(x(s1),xp(i2),y1(i3),y2(i4),Tn.cplx(xp(i2),y1(i3),y2(i4)));
				}

			}
			vector<double>().swap(MM1);
			vector<int>().swap(in1);
			Tn = Tk;
			Tk = ITensor(0);
			}

		}


		else
		{
			Tn = T*T1;
			auto C1 = combiner(x1,x1p);
			auto C2 = combiner(x2,x2p);
			Tn = Tn*C1;
			Tn = Tn*C2;
			auto c1 = commonIndex(Tn,C1);
			auto c2 = commonIndex(Tn,C2);
			Tn = Tn*delta(c1,x);
			Tn = Tn*delta(c2,xp);
		}


		UU = ITensor(0);
		UU1 = ITensor(0);
		M = ITensor(0);


		P = norm(Tn);
		Ren[ind] = P;
		ind++;

		Tn = Tn / P;

		T = Tn;

		Tn = ITensor(0);

		s1 = 0; s2 = 0;

		dx = dxn;

		x1 = Index("x1",dx,Xtype);
		x2 = prime(x1,1);

		T = T * delta(xp,x2);
		T = T * delta(x,x1);

		Print(T);

/*				HOTRG : x-direction starts here.             */


		x1 = noprime(findtype(T,Xtype));
		x2 = prime(x1,1);
		y1 = noprime(findtype(T,Ytype));
		y2 = prime(y1,1);
		y3 = Index("y3",dy);
		y4 = Index("y4",dy);

		dyn = (dcut > dy*dy) ? dy*dy : dcut;	
		
		T1 = T;


		y1p = Index("y1p",dy);
		y2p = prime(y1p,1);
		xx1 = Index("xx1",dx,Xtype);

		i = Index("i",dx);
		j = Index("j",dx);

		T1 *= delta(y1,y1p);
		T1 *= delta(y2,y2p);

//		T = T * delta(x2,i); T = T * delta(x1,x1); T1 = T1 * delta(x1,i);

		T = T * delta(x2,i);
		T = T * delta(x1,xx1);
		T1 = T1 * delta(x1,i);

		MT1 = T;
		temp = T1;
		MT1 = (MT1*delta(y1,y3)*delta(i,j));
		temp = (temp*delta(y1p,y4)*delta(i,j));
		MT1.conj();	temp.conj();

		MT1 = MT1*T;
		temp = temp*T1;
		MT1 = MT1*temp;
		MT1 = MT1.takeReal();


		MT2 = T;
		temp = T1;
		MT2 = (MT2*delta(y2,y3)*delta(i,j));
		temp = (temp*delta(y2p,y4)*delta(i,j));
		MT2.conj();	temp.conj();

		MT2 = MT2*T;
		temp = temp*T1;
		MT2 = MT2*temp;
		MT2 = MT2.takeReal();


//		cout << "x-d step 1 end : " << double(clock() - tt)/CLOCKS_PER_SEC << " sec" << endl;
		
		ui = Index("A",dy*dy);
		uii = ui;

		MT2 = MT2*delta(y2,y1)*delta(y2p,y1p);

//		cout << "Norm: " << norm(MT1-MT2) << endl;

		MT2 = MT2*delta(y1,y2)*delta(y1p,y2p);


		MM1.resize(dy*dy*dy*dy);
		MM2.resize(dy*dy*dy*dy);
		Ss1.resize(dy*dy);
		Ss2.resize(dy*dy);
		in1.resize(dy*dy);
		in2.resize(dy*dy);

		y = Index("y",dy*dy,Btype);
		ii1 = Index("ii1",dy);
		ii2 = Index("ii2",dy);


		if(dyn == dcut)
		{
			for(auto i1 : range1(dy))
			for(auto i2 : range1(dy))
			for(auto i3 : range1(dy))
			for(auto i4 : range1(dy))
			{
				MM1[(i2-1+(i1-1)*dy)+(i4-1+(i3-1)*dy)*dy*dy] = MT1.real(y1(i1),y1p(i2),y3(i3),y4(i4));
				MM2[(i2-1+(i1-1)*dy)+(i4-1+(i3-1)*dy)*dy*dy] = MT2.real(y2(i1),y2p(i2),y3(i3),y4(i4));
			}
//			cout << "ok" << endl;	
			MT1 = ITensor(0);
			MT2 = ITensor(0);
//			heev_f('V','U',dy*dy,MM1.data(),dy*dy,Ss1.data());
//			heev_f('V','U',dy*dy,MM2.data(),dy*dy,Ss2.data());
			LAPACKE_dsyev( LAPACK_COL_MAJOR,'V', 'U', dy*dy,MM1.data(),dy*dy,Ss1.data());
			LAPACKE_dsyev( LAPACK_COL_MAJOR,'V', 'U', dy*dy,MM2.data(),dy*dy,Ss2.data());


			for(int q=0;q<dy*dy;q++)
			{
				in1[q]=q; in2[q]=q;
			}

			for(int q=0;q<dy*dy-1;q++)
			for(int p=0;p<dy*dy-1;p++)
			if(fabs(Ss1[in1[p]])>fabs(Ss1[in1[p+1]]))
			{
				int tmp = in1[p];
				in1[p] = in1[p+1];
				in1[p+1] = tmp;
			}

			for(int q=0;q<dy*dy-1;q++)
			for(int p=0;p<dy*dy-1;p++)
			if(fabs(Ss2[in2[p]])>fabs(Ss2[in2[p+1]]))
			{
				int tmp = in2[p];
				in2[p] = in2[p+1];
				in2[p+1] = tmp;
			}
			

			for(auto i1 : range1(dy*dy-dcut))
			{
//				cout << "s3\t" << Ss1[in1[i1]] << endl;
//				cout << "s4\t" << Ss2[in2[i1]] << endl << endl;
				s1 += fabs(Ss1[in1[i1]]);
				s2 += fabs(Ss2[in2[i1]]);
			}
//			cout << "Sum ? " << endl;
//			cout << "s3 = " << s1 << endl;
//			cout << "s4 = " << s2 << endl;
			if(s1<s2)
			{
//				cout << "x-direction : choose UU" << endl;
				ui = Index("A",dcut);
				UU = ITensor(ii1,ii2,ui);
				for(auto i1 : range1(dy))
				for(auto i2 : range1(dy))
				for(auto i3 : range1(dcut))
					UU.set(ii1(i1),ii2(i2),ui(i3),MM1[(i2-1+(i1-1)*dy)+(in1[dy*dy-i3])*dy*dy]);
//				vector<double>().swap(MM1);
				vector<double>().swap(Ss1);
				vector<double>().swap(MM2);
				vector<double>().swap(Ss2);
				sindex=1;
			}
			else
			{
//				cout << "x-direction : choose UD" << endl;
				ui = Index("A",dcut);
				UU = ITensor(ii1,ii2,ui);
				for(auto i1 : range1(dy))
				for(auto i2 : range1(dy))
				for(auto i3 : range1(dcut))
					UU.set(ii1(i1),ii2(i2),ui(i3),MM2[(i2-1+(i1-1)*dy) + (in2[dy*dy-i3])*dy*dy]);
				vector<double>().swap(MM1);
				vector<double>().swap(Ss1);
//				vector<double>().swap(MM2);
				vector<double>().swap(Ss2);
				sindex=0;
			}
//
		}
		else
		{
			for(auto i1 : range1(dy))
			for(auto i2 : range1(dy))
			for(auto i3 : range1(dy))
			for(auto i4 : range1(dy))
			{
				MM1[(i2-1+(i1-1)*dy)+(i4-1+(i3-1)*dy)*dy*dy] = MT1.real(y1(i1),y1p(i2),y3(i3),y4(i4));
			}
			MT1 = ITensor(0);
			MT2 = ITensor(0);
			LAPACKE_dsyev( LAPACK_COL_MAJOR,'V', 'U', dy*dy,MM1.data(),dy*dy,Ss1.data());
			ui = Index("A",dy*dy);
			UU = ITensor(ii1,ii2,ui);
			for(auto i1 : range1(dy))
			for(auto i2 : range1(dy))
			for(auto i3 : range1(dy*dy))
				UU.set(ii1(i1),ii2(i2),ui(i3),MM1[(i2-1+(i1-1)*dy) + (dy*dy-i3)*dy*dy]);
			vector<double>().swap(MM1);
			vector<double>().swap(Ss1);
			vector<double>().swap(MM2);
			vector<double>().swap(Ss2);
		}



		y = Index("y",dyn,Btype);
		yp = Index("yp",dyn,Btype);
		Tk = ITensor(y,yp,xx1,x2);
		if(dyn==dcut)
		{
		if(sindex==0)
		{
		UU1 = UU;
		auto j1 = Index("j1",dy);
		auto j2 = Index("j2",dy);
		UU1 = UU1 * delta(ii2,j2);
		UU1 = UU1 * delta(ii1,j1);


//		UU *= delta(ui,y);
		UU1 *= delta(ui,yp);

		T = T*delta(y2,j1);	T = T*delta(y1,ii1);
		T1 = T1*delta(y2p,j2);	T1 = T1*delta(y1p,ii2);

		for(int s1=1;s1<=dyn;s1++)
		{
		auto Uk = ITensor(ii1,ii2);

		for(auto i1: range1(dy))
		for(auto i2: range1(dy))
			Uk.set(ii1(i1),ii2(i2),MM2[(i2-1+(i1-1)*dy) + (in2[dy*dy-s1])*dy*dy]);

		Tn = Uk*T;
		Tn = Tn*T1;
		Tn = Tn*UU1.conj();

		for(auto i2: range1(int(dyn)))
		for(auto i3: range1(dx))
		for(auto i4: range1(dx))
		{
			Tk.set(y(s1),yp(i2),xx1(i3),x2(i4),Tn.cplx(yp(i2),xx1(i3),x2(i4)));
		}
		}
		vector<double>().swap(MM2);
		vector<int>().swap(in2);
		Tn = Tk;
		Tk = ITensor(0);
		}
		else
		{
		UU1 = UU; 
		auto j1 = Index("j1",dy);
		auto j2 = Index("j2",dy);
		UU1 = UU1 * delta(ii2,j2);
		UU1 = UU1 * delta(ii1,j1);
 
//			UU *= delta(ui,y);
		UU1 *= delta(ui,yp);

		T = T*delta(y2,ii1); T = T*delta(y1,j1);
		T1 = T1*delta(y2p,ii2);  T1 = T1*delta(y1p,j2);
 
		for(int s1=1;s1<=dyn;s1++)
		{
		auto Uk = ITensor(ii1,ii2);

		for(auto i1: range1(dy))
		for(auto i2: range1(dy))
			Uk.set(ii1(i1),ii2(i2),MM1[(i2-1+(i1-1)*dy) + (in1[dy*dy-s1])*dy*dy]);


			Tn = Uk*T;
			Tn = Tn*T1;
			Tn = Tn*UU1.conj();

		for(auto i2: range1(int(dyn)))
		for(auto i3: range1(dx))
		for(auto i4: range1(dx))
                {
			Tk.set(y(s1),yp(i2),xx1(i3),x2(i4),Tn.cplx(yp(i2),xx1(i3),x2(i4)));
		}
		}
		vector<double>().swap(MM2);
		vector<int>().swap(in2);
		Tn = Tk;
		Tk = ITensor(0);



		}
		}

		else
		{
			Tn = T*T1;
			auto C3 = combiner(y1,y1p);
			auto C4 = combiner(y2,y2p);
			Tn = Tn*C3;
			Tn = Tn*C4;
			auto c3 = commonIndex(Tn,C3);
			auto c4 = commonIndex(Tn,C4);
			Tn = Tn*delta(c3,y);
			Tn = Tn*delta(c4,yp);
		}
/*
		Print(UU);
		Print(UU1);
		Print(M);
		Print(Tn);
*/

		T = ITensor(0);
		UU = ITensor(0);
		UU1 = ITensor(0);

//		cout << "x-d step 3 end : " << double(clock() - tt)/CLOCKS_PER_SEC << " sec" << endl;

//		Print(Tn);

		P = norm(Tn);
//		P = fabs((Tn.cplx(xx1(1),x2(1),y(1),yp(1))));
		Ren[ind] = P;
		ind++;

		Tn = Tn / P;

		T = Tn;
		
		Tn = ITensor(0);

		dy = dyn;

		y1 = Index("y1",dy,Ytype);
		y2 = prime(y1,1);

		T = T * delta(xx1,x1);
		T *= delta(y,y1);
		T *= delta(yp,y2);

		Print(T);

	}

	complex<double> Z = 0;

	x1 = noprime(findtype(T,Xtype));
	y1 = noprime(findtype(T,Ytype));
	x2 = prime(x1,1);
	y2 = prime(y1,1);

//////		log(Tr(Tn)) ///////
	for(auto s1 : range1(dx))
	for(auto s2 : range1(dy))
		Z = Z + (T.cplx(x1(s1),x2(s1),y1(s2),y2(s2)));

	T = ITensor(0);
	 
	Z = log(Z)/pow(2,2*L);

//	cout << Z << endl;

	int mul = pow(4,L);

	for(int p=0;p<2*L;p++)
	{
//		cout << Ren[p] << endl;

//		mul = mul/2;
//		for(int q=0;q<mul;q++)
//			Z = Z*Ren[p];
		Z = Z + log(Ren[p])/pow(2,p+1);
	}


	return Z;
}


/*
void heev_f(char jobz, char uplo, int n, dcomplex* a, int lda, double* w)
{
    int lwork, info;
    dcomplex tmp;
    double* rwork=new double[3*n-2];
    lwork=-1;
    zheev(&jobz,&uplo,&n,a,&lda,w,&tmp,&lwork,rwork,&info);
    lwork=(int)(tmp.real()+0.1);
    dcomplex* work=new dcomplex[lwork];
    zheev(&jobz,&uplo,&n,a,&lda,w,work,&lwork,rwork,&info);
    delete[] work;
    delete[] rwork;
    assert(info==0);
}
*/
