import analysator as pt
import numpy as np
# import vlsvrs
import os
import hashlib
import pickle
import importlib

datalocation = "/turso/group/spacephysics/analysator/CI/analysator-test-data/vlasiator/"
files=["3D/FID/bulk1/bulk1.0000995.vlsv",
       "3D/FHA/bulk1/bulk1.0000990.vlsv",
]


class Tester:
    def __init__(self,filename=None):
        self.filename=filename
        self.vlsvobj=None
        self.hashes_dict_rust={}
        self.hashes_dict_python={}
    def changeFile(self,filename):
        self.filename=filename
    def loadPickle(self,file):
        self.pickled = pickle.load(file)
    def dumpPickle(self,file):
        pickle.dump(self.hashes_dict,file)

    def load(self,backend=None):
        if not backend or backend.lower()=="rust":
            raise NotImplementedError()
            # self.vlsvobj_rust=vlsvrs.VlsvFile(self.filename)
        if not backend or backend.lower()=="python":
            self.vlsvobj_python=pt.vlsvfile.VlsvReader(self.filename)

    def setHashTarget(self,backend):
        if backend=='rust':
            self.vlsvobj=self.vlsvobj_rust
        elif backend=='python':
            self.vlsvobj=self.vlsvobj_python
        else:
            print("None set, give valid backend")

    def hash(self,func,args,op=None,opargs=None,both=False,loop=False,flatten=True,sort=False):
            
        def update(vlsvobj,op,opargs,args,hashdict,loop=False):
            #If we want to repeat same function func with different arguments
            if loop:
                for arg in args:
                    print(arg,args)
                    update(vlsvobj,op,opargs,arg,hashdict)
                return 0
            argkey=str(args)
            opsname="_"+str(op)+"_"+str(opargs)
            #Get the method of the vlsvobj that matches the given func str
            t=getattr(vlsvobj,func)
            #Handle arguments and call the function with the given args to get return value
            if type(args) is dict:
                retval=t(**args)
            elif type(args) is list:
                retval=t(*args)
            else: 
                raise IOError(f"Wrong args type: {type(args)} {args}")
            #If we want to do operations on the retval for example reshaping, type chaning or sorting
            if op and opargs:
                #Make into list for handling
                if type(op) is not list:
                    op=[op]
                    opargs=[opargs]

                for i,f in enumerate(op):
                    try:
                        fun=getattr(retval,f)
                    except AttributeError:
                        try: 
                            #if given function is not method of retval we make retval the argument of function
                            if '.' in f:
                                funcl=f.split('.')
                                #in case it is inside a module like numpy we need to get instance of the module
                                funcl[0]=importlib.import_module(funcl[0])

                                fun=getattr(funcl[0],funcl[1])
                            else:
                                fun=f
                            opargs[i]=[retval]
                        except AttributeError as e:
                            raise AttributeError(f"Did not find func {func} to operate with: {e}")
                    retval=fun(*opargs[i])

            #save hash of the retval as array
            retval=np.array(retval)
            print(retval.shape)
            if flatten:
                retval.reshape((-1,))
            if sort:
                retval.sort()
            if self.filename not in hashdict.keys():
                hashdict[self.filename]={}
                if func not in hashdict:
                    hashdict[self.filename][func]={}
            hashdict[self.filename][func][argkey]=[hashlib.sha256(retval.tobytes()).hexdigest(),opsname]

        if not both:
            if self.vlsvobj==self.vlsvobj_python:
                hashdict=self.hashes_dict_python
            elif self.vlsvobj==self.vlsvobj_rust:
                hashdict=self.hashes_dict_rust
            update(self.vlsvobj,op,opargs,args,hashdict,loop)
        else:
            update(self.vlsvobj_rust,op,opargs,args,self.hashes_dict_rust,loop)
            update(self.vlsvobj_python,op,opargs,args,self.hashes_dict_python,loop)

    def compare(self,funcpy,argspy,funcrust,argsrust):
        try:
            py=getattr(self.vlsvobj_python,funcpy)
            retval_py=py(**argspy)

            rust=getattr(self.vlsvobj_rust,funcrust)
            retval_rust=rust(**argsrust)

        except Exception as e:
            raise e

        if type(retval_py) is dict and type(retval_rust) is dict :
            stack=list(retval_rust.keys())
            if (len(retval_py)!=len(retval_rust)) and len(list(retval_py.keys()))!=0:
                raise SystemError("one or both of the dictionaries returned by the readers are empty")
            for key in retval_py.keys():
                if retval_rust[key] == retval_py[key]:
                    stack.remove(key) #maybe a some ohter way to remove it is faster? 
                else:
                    raise SystemError("returned dictionary values between vlsvreader and vlsvrs do not match")
            if len(stack)!=0:
                raise KeyError("returned dictionry from vlsvrs contains keys not present in dictonary returned by python.") 
            return True
        else:
            raise NotImplementedError 
        

#read ref from file
#(Cellid as variable) reading single cellid ,reading cellid list, (input and output cellid should be same assertion).
#vector (tensor variables from datareduction), read_variable, cellid 0 is nonexistant, error check.
#some datareduction
#read_vdf() (vlsvrs)
# read_variable, compare against fsgird, vg 
# read_velocitycells (vlsvreader)

# read_vdf_spares 
# 
ciTester = Tester()
# files=["s"]
for file in files:

    #Load data 
    filename=os.path.join(datalocation,file)

    ciTester.changeFile(filename)
    ciTester.load()

    #Test compare
    # cid=ciTester.vlsvobj_python.get_cellid_with_vdf(np.array([0,0,0]))
    # ciTester.compare("read_velocity_cells",{"cellid":cid,"pop":"proton"},"read_vdf_sparse",{"cid":cid,"pop":"proton"})
    
    variables_to_test=["CellID","vg_rhom","vg_v","vg_rhoq","proton/vg_rho","proton/vg_v"] #fg_variable read issue with read_variable
    variables_to_test_nonraw=["fg_b","fg_v"]
    
    pylist=ciTester.vlsvobj_python.get_variables()
    variables=[[var] for var in variables_to_test if (var in pylist)] 
    nonraw_vars=[[var,0] for var in variables_to_test_nonraw if (var in pylist)] 
    if False: #Currently vlsvrs is not being used, this will be updated later when that is dependency to analysator 
        #Make hash rust
        ciTester.setHashTarget("rust")
        #ciTester.hash("read_variable",{"variable":"CellID","op":0},op=["reshape","astype","numpy.sort"],opargs=[[tuple([-1])],[int],[]],sort=False,flatten=False)
        pylist=ciTester.vlsvobj_python.get_variables()
        rustlist=ciTester.vlsvobj_rust.list_variables()
        variables=[[var] for var in variables_to_test if (var in pylist and var in rustlist)] 
        nonraw_vars=[[var,0] for var in variables_to_test_nonraw if (var in pylist and var in rustlist)] 
        ciTester.hash("read_variable_raw",variables,loop=True)
        ciTester.hash("read_variable",nonraw_vars,loop=True)

    #Make hash python
    ciTester.setHashTarget("python")
    variables.extend([[var[0]] for var in nonraw_vars]) #prob some prettier way than looping through it all but it's not a big list
    ciTester.hash("read_variable",variables,loop=True)
    
print(ciTester.hashes_dict_python)

quit()
#Should not be used yet
retval=0
key_map_rust_to_py={"read_variable_raw":"read_variable","read_variable":"read_variable"}
for file in ciTester.hashes_dict_rust.keys():
    print(f"------{file}------")
    for key in ciTester.hashes_dict_rust[file].keys():
        py_dict=ciTester.hashes_dict_python[file][key_map_rust_to_py[key]]
        rust_dict=ciTester.hashes_dict_rust[file][key]
        for argcall in rust_dict.keys():
            if rust_dict[argcall][0]!=py_dict[argcall][0]:
                print(rust_dict[argcall][0],py_dict[argcall][0])
                print(f"Hashes do not match for call {argcall}!")
                retval=1
            else:
                print("Match")

if retval==1:
    raise SystemError("Some hashes did not match")
 
