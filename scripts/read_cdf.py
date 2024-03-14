from cdflib import CDF
import numpy as np


'''
CDF.attget(         CDF.cdf_info(       CDF.epochrange(     CDF.increment       CDF.print_attrs(    CDF.varattsget(     CDF.varinq(         CDF.version
CDF.attinq(         CDF.close(          CDF.globalattsget(  CDF.mro(            CDF.release         CDF.varget(         CDF.vdr_info(
'''

def read_cdf(filename, silent = False):
    '''
    read_cdf() is a generic routine for reading '.cdf' files 
    --- duplicates IDL routine read_cdf.pro by Kosta Horaites
    returns all variables into a structure

    EXAMPLE CALL: a = read_cdf('test.cdf')

    NOTES:
    There are two different kinds of data stored in cdf files: rvariables and zvariables.
    rvariables, or "Regular" variables, have to conform to a rigid structures
    (all rvariables in a file have to have the same dimension, which is defined ahead of time). 
    zvariables are more flexible (no importance in the letter "z", just a name).

    Reading CDF Files
    The following built-in IDL commands are used to read data from a CDF file:
    CDF_OPEN: Open an existing CDF file.
    CDF_INQUIRE: Call this function to find the general information about the contents of the CDF file.
    CDF_CONTROL: Call this function to obtain further information about the CDF file
    CDF_VARINQ: Retrieve the names, types, sizes, and other information about the variables in the CDF file.
    CDF_VARGET: Retrieve the variable values.
    CDF_ATTINQ: Optionally, retrieve the names, scope and other information about the CDFs attributes.
    CDF_ATTGET: Optionally, retrieve the attributes.
    CDF_CLOSE: Close the file.
    If the structure of the CDF file is already known, the inquiry routines do not need to be called--
    only CDF_OPEN, CDF_ATTGET, CDF_VARGET, and CDF_CLOSE would be needed.
    
    '''
    c = CDF(filename)
    str = c.cdf_info()

    nnvars = len(str['rVariables'])
    nnzvars = len(str['zVariables'])

    ndim = np.zeros(nnvars + nnzvars, dtype = int)
    rec_counts = np.zeros(nnvars + nnzvars, dtype = int)  # "record counts" (i.e. number of measurements in the file)
    names = [''] * (nnvars + nnzvars)
    vartypes = [''] * (nnvars + nnzvars)

    if nnvars + nnzvars > 0:
        for i in range(nnvars + nnzvars):
            varstr = c.varinq(i)
            names[i] = varstr['Variable']
            ndim[i] = varstr['Num_Dims']        # or Num_Dims +1 ?
            vartypes[i] = varstr['Var_Type']     # 'rVariable' or 'zVariable'
            if not silent:
                print(varstr['Variable'])
            rec_counts[i] = varstr['Last_Rec']

    output = {}

    #for i in range(len(rec_counts)):
    #    if rec_counts[i] == 0:
    #        rec_counts[i] = 1   #   prevents some kind of error with CDF_VARGET

    for i, name in enumerate(names):
        if name != '':
            output[name] = c.varget(name)

    c.close()

    return output
 
