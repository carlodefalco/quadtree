all: __advdiff_local_matrix__.oct __laplacian_local_matrix__.oct __reaction_local_matrix__.oct

%.oct : %.cc
	mkoctfile $<

clean:
	$(RM) *.o

distclean: clean
	$(RM) *.oct
