from distutils.core import setup, Extension

gsl_dir = '/home/rfadely/local/'
src_dir = './'

extmod1 = Extension('_model_maker',
                    include_dirs = [gsl_dir+'include'],
                    libraries = ['gsl'],
                    library_dirs = [gsl_dir+'lib'], 
                    sources = [src_dir+'_model_maker.c',
                               src_dir+'calc_normalization.c',
                               src_dir+'get_filelength.c',
                               src_dir+'get_num_files.c',
                               src_dir+'integrate_sed.c',
                               src_dir+'read_file.c',
                               src_dir+'regrid.c'])
extmod2 = Extension('_fit_models',
                    sources = [src_dir+'_fit_models.c'])
extmod3 = Extension('_coeff_marginalization',
                    sources = [src_dir+'_coeff_marginalization.c'])
extmod4 = Extension('foo',
                    sources = [src_dir+'foo.c',
                               src_dir+'get_num_files.c',
                               src_dir+'calc_normalization.c',
                               src_dir+'read_file.c',
                               src_dir+'get_filelength.c'])

setup(ext_modules=[extmod1,extmod2,extmod3,extmod4])
