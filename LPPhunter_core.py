# -*- coding: utf-8 -*-
# Copyright 2014-16 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release.
# For more info please contact zhixu.ni@uni-leipzig.de

import ConfigParser
from LibLPPhunter.xic import XIC


config = ConfigParser.ConfigParser()
config.read('config.ini')

infile_type = config.get('inputfile', 'filetype')
if infile_type.lower() == 'mzml':
    infile_name = config.get('inputfile', 'filename')
    print infile_name

    xic_spec = XIC(infile_name)
    xic_spec.find_mz(778.560)

else:
    print 'No input mzML'