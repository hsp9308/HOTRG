1. Libraries needed

- ITensor ver.2 (https://github.com/ITensor/ITensor)
 (Used in the Fisher zero calculation and later YLzero calculation.)
- Complex bessel function (https://github.com/joeydumont/complex_bessel)
- Intel MKL 2017, 2018 (Not certain for the later versions.)

2. How to compile 

(1) Install the libraries above (make sure that Intel MKL is correctly linked in the ITensor installation)

(2) modify the Makefile;
correct the variable LIBRARY_DIR to one's itensor installation directory.

(3) Done! run 'make'. 

3. What is opt.py ?

opt.py finds the partition function zero by taking output complex number of './TRG_it' and minimize a function of partition function. 

So you may run 'python opt.py' to find the zeros with Nelder-Mead method.

You need to set the initial position of parameters (here, (x1, y1) as a complex temperature). 
