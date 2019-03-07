#CFLAGS = -Wall  -I ../extern/include/ -L ../extern/lib/  
#LFLAGS = -Wall  -I ../extern/include/ -L ../extern/lib/ -lgsl -lcblas -latlas -llapack -lblas -lF77 -lI77 -lc -lm 

LFLAGS = -Wall -lgsl -lgslcblas -latlas -llapack -lc -lm

test: test.c segSites.c ../hmm/numerical.c ../hmm/popGenTools.c ../hmm/nrutil.c ../hmm/adkGSL.c
	gcc test.c segSites.c ../hmm/numerical.c ../hmm/popGenTools.c  ../hmm/nrutil.c ../hmm/adkGSL.c -o test $(CFLAGS) 

segSiteHMM3: segSites.c segSiteHMM3.c ../hmm/hmm.c ../hmm/popGenTools.c ../hmm/adkGSL.c ../hmm/numerical.c ../hmm/nrutil.c 
	gcc segSites.c segSiteHMM3.c ../hmm/hmm.c ../hmm/popGenTools.c ../hmm/adkGSL.c ../hmm/numerical.c ../hmm/nrutil.c -o segSiteHMM3 $(LFLAGS)
clean:
	rm test
	rm segSiteHMM
	rm segSiteHMM3
