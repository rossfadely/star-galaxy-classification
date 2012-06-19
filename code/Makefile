
CC = gcc

CFLAGS = -Wall


INCLUDES = -I/usr/local/scisoft/include -I/usr/local/scisoft/packages/gsl/include


LFLAGS = -L/usr/local/scisoft/lib/  -L/usr/local/scisoft/packages/gsl/lib
LIBS = -lgsl -lgslcblas -lcfitsio -lm


SRCS = HBSGsep.c filter_calcs.c get_num_files.c get_filelength.c read_file.c regrid_filter.c calc_normalization.c data_calcs.c count_data.c printfitserror.c calc_datavals.c star_model_calcs.c regrid_sed.c integrate_sed.c gal_model_calcs.c fit_calcs.c coeff_calcs.c fit_star_template.c coeff_meanvar.c fit_gal_template.c sparse_calcs.c calc_P_F_kS.c calc_P_F_kG.c optimize.c read_hyppars.c loglikelihood.c calc_P_F_S.c calc_P_F_G.c write_lnprob.c write_hyppars.c write_lnlikeratio.c star_model_mags.c gal_model_mags.c write_modelmags.c write_minchi2.c

OBJS = $(SRCS:.c=.o)

MAIN = HBSGsep

all:    $(MAIN) 
	@ rm -f *.o

$(MAIN): $(OBJS) 
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

doc: 	code_doc.pdf

%.pdf: %.tex
	pdflatex $<
	bash -c " ( grep Rerun $*.log && pdflatex $< ) || echo noRerun "
	bash -c " ( grep Rerun $*.log && pdflatex $< ) || echo noRerun "
