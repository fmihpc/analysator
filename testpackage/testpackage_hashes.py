import analysator as pt
import numpy as np
# import vlsvrs
import os
import hashlib
import pickle
import importlib
import argparse
datalocation = "/turso/group/spacephysics/analysator/CI/analysator-test-data/vlasiator/"
files=[
        "3D/FID/bulk1/bulk1.0000995.vlsv",
       "3D/FHA/bulk1/bulk1.0000990.vlsv",
       "2D/BGA/zero_ehall_layers_23/bulk.0000380.vlsv"
]

parser=argparse.ArgumentParser(
        prog='ciTester',
        description='Used for generating and comparing hashes from analysator output for CI testing.')
parser.add_argument('--generate','-g',help="Generate hashdump into given path",const=str,nargs='?')
parser.add_argument('--compare','-c',help="Compare against given hashdump",const=str,nargs='?')
args=parser.parse_args()
generate_path=args.generate
compare_path=args.compare

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

    def dumpIntoFile(self,indict,dumpname):
        path_to_file=os.path.join(generate_path,dumpname)
        if os.path.isfile(path_to_file):
            print(f"::warning:: {dumpname} found in {generate_path}, overwriting!")
        with open(path_to_file,"w") as file:
            for (filename,funcdict) in indict.items():
                file.write("File: "+filename+'\n')
                for (funccall,hashdict) in funcdict.items():
                    file.write("\tFunction: "+funccall+'\n')
                    for (arg,hash_and_op) in hashdict.items():
                        hash=hash_and_op[0]
                        op=hash_and_op[1]
                        file.write("\t\t"+f"{arg:<30} {hash} {op}"+'\n')
            file.close()
        assert(self.loadFromFile(dumpname,generate_path)==indict)

    def loadFromFile(self,dumpname,compare_path=compare_path):
        outdict={}
        path_to_file=os.path.join(compare_path,dumpname)
        if not os.path.isfile(path_to_file):
            print(f"::warning::Tried to load file {path_to_file} but no file found, skipping this dictionary")
            return 1
        with open(path_to_file,"r") as file:
            for line in file:
                line=line.rstrip('\n')
                if line[:5]=="File:":
                    filename=line.split(" ")[1]
                    if filename not in outdict:
                        outdict[filename]={}
                    continue
                if "Function:" in line and "[" not in line: # bit stupid but works for now
                    function=line.split(" ")[1]
                    #There should not be multiples of function with same filename!
                    try:
                        outdict[filename][function]={}
                    except KeyError:
                        raise IOError("Invalid format of the input file")

                elif "[" in line or "NOTARG" in line: 
                    #above is bit stupid but should filter it little bit since the hash lines should have a list in them
                    listHashInfo=[item.strip('\t') for item in line.split(" ") if item!=""]
                    outdict[filename][function][listHashInfo[0]]=[listHashInfo[1],listHashInfo[2]]
                else:
                    #We should never end up here, but if we do, the dict was not read correctly likely badly formatted
                    raise IOError("File is likely not formatted correctly")
            file.close()
        return outdict

    def loadobj(self,backend=None):
        #if not backend or backend.lower()=="rust":
        #     self.vlsvobj_rust=vlsvrs.VlsvFile(self.filename)
        if not backend or backend.lower()=="python":
            self.vlsvobj_python=pt.vlsvfile.VlsvReader(self.filename)

    def setHashTarget(self,backend):
        if backend=='rust':
            self.vlsvobj=self.vlsvobj_rust
        elif backend=='python':
            self.vlsvobj=self.vlsvobj_python
        else:
            print("None set, give valid backend")

    def hash(self,func,args,op=None,opargs=None,both=False,loop=False,flatten=True,sort=False,argkey_name=None):
            
        def update(vlsvobj,op,opargs,args,hashdict,loop=False):
            #If we want to repeat same function func with different arguments
            if loop:
                for arg in args:
                    update(vlsvobj,op,opargs,arg,hashdict)
                return 0
            if argkey_name:
                argkey=str(argkey_name+"_NOTARG")
            else:
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
            if flatten:
                retval.reshape((-1,))
            if sort:
                retval.sort()
          
            if self.filename not in hashdict.keys():
                hashdict[self.filename]={}
            if func not in hashdict[self.filename]:
                hashdict[self.filename][func]={}

            hashdict[self.filename][func][argkey]=[hashlib.sha256(retval.tobytes()).hexdigest(),opsname]

        if not both:
            #of course fails if one or the other is not defined
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
    def interpolationtest2d(self,varname):
        N = 1000#int(np.sqrt(800))
        delta = 60e6
        xmin = 45.0e6
        xcoords = np.linspace(xmin,xmin+delta,N)
        ymin = -37.51e6 - 1e7*0
        ycoords = np.linspace(ymin,delta+ymin,N)
        X,Y,Z = np.meshgrid(xcoords, ycoords,np.array([-0.25e6]))
        ncoords = np.prod(X.shape)

        coords = np.hstack((np.reshape(X,(ncoords))[:,np.newaxis],
                            np.reshape(Y,(ncoords))[:,np.newaxis],
                            np.reshape(Z,(ncoords))[:,np.newaxis]))

        self.hash("read_interpolated_variable",[varname,coords],argkey_name=varname)
    def interpolationtest3(self):
        RE=6371e3
        coords=[[5*RE,RE,0.5*RE],np.array([[10*RE,RE,0.5*RE],[5*RE,RE,0.1*RE]]),np.array([[5*RE,RE,0.5*RE],[8*RE,RE,0.1*RE]])]
        for i,coord in enumerate(coords):
            self.hash("read_interpolated_variable",["proton/vg_rho", coord],argkey_name=f"proton/vg_rho_{i}")
            self.hash("read_interpolated_variable",["proton/vg_v", coord],argkey_name=f"proton/vg_v_{i}")
            self.hash("read_interpolated_variable",["proton/vg_ptensor",coord],argkey_name=f"proton/vg_ptensor_{i}") 
#read ref from file
#(Cellid as variable) reading single cellid ,reading cellid list, (input and output cellid should be same assertion).
#vector (tensor variables from datareduction), read_variable, cellid 0 is nonexistant, error check.
#some datareduction
#read_vdf() (vlsvrs)
# read_variable, compare against fsgird, vg 
# read_velocitycells (vlsvreader)

# read_vdf_spares 
# 
#
if __name__=="__main__":
    ciTester = Tester()
    for file in files:

        #Load data 
        filename=os.path.join(datalocation,file)

        ciTester.changeFile(filename)
        ciTester.loadobj()
        
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
        if "proton/vg_v" in pylist:
            ciTester.interpolationtest2d("proton/vg_v")
            ciTester.interpolationtest3();
        variables.extend([[var[0]] for var in nonraw_vars]) #prob some prettier way than looping through it all but it's not a big list
        ciTester.hash("read_variable",variables,loop=True)

    if generate_path:
        if ciTester.hashes_dict_python is not None:
            ciTester.dumpIntoFile(ciTester.hashes_dict_python,"hashdump_python.txt")
            os.system(f"cat {os.path.join(generate_path,'hashdump_python.txt')}")
        if ciTester.hashes_dict_rust is not None:
            ciTester.dumpIntoFile(ciTester.hashes_dict_rust,"hashdump_rust.txt")
            os.system(f"cat {os.path.join(generate_path,'hashdump_rust.txt')}")
            
    if compare_path:
        dumps=["hashdump_python.txt","hashdump_rust.txt"]
        for i,hashdump in enumerate(dumps):
            refDict=ciTester.loadFromFile(hashdump)
            if ciTester.hashes_dict_python is not None and i==0:
                hashdict=ciTester.hashes_dict_python
            elif ciTester.hashes_dict_rust is not None and i==1:
                hashdict=ciTester.hashes_dict_rust 
            else:
                print(f"::warning:: ciTester does not have hash dictionary to compare against {hashdump}.")
                continue

            if refDict!=1:
               if refDict!=hashdict: #ordering should not matter for comparison of dictionaries 
                   raise SystemError("Hashes differ")
               elif refDict==hashdict:
                   print("Hash dictionaries match")
            else:
               retval=1

    quit()
#Should not be used yet
    if retval!=1:
        retval=0
    key_map_rust_to_py={"read_variable_raw":"read_variable","read_variable":"read_variable"} #function calls may not match, can be used to map from rust vlsvrs calls to py calls
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
                    continue

    if retval==1:
        raise SystemError("Some hashes did not match")
    
