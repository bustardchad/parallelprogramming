# parallelprogramming
Chad Bustard, CS 759 Final Project

The powerpoint presentiation is CS759_final.pptx, and the final report is CS759_FinalReport.pdf.

Codes: in the folder called 'code' 
	sequential.cpp - sequential implementation
	parallel.cpp - naive parallel implementation
	parallel2.cpp - fastest parallel implementation with only outer for loop parallelized
	parallel3.cpp - same as parallel2.cpp but with temporary storage unew(x,y) for u(x,y)
	wavefront.cpp - wavefront method

Executables: 
	sequential
	parallel1
	parallel2
	parallel3
	wavefront
When run, each executable makes a file with the initial grid set-up, as well as the final solution (uFinal_parallel2.txt, for example). I made a MatLab program named plotU.m that will plot the result. Each of them converges to thesolution using Mathematica, which is plotted in poissonMathematica.png.

Timing Results: all are in the timing analysis folder
	timeAnalysis_wave.png has the timing for each implementation. I used 8 processors and a tolerance of 0.1.
