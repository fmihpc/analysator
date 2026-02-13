import cpphelpers as chp
import analysator as pt
import numpy as np
import time
def buildDescriptor(f):
    now = time.time() 
    f._VlsvReader__read_fileindex_for_cellid()
    print("read",time.time()-now)
    fileindex_for_cellid = f._VlsvReader__fileindex_for_cellid
    xc = f._VlsvReader__xcells
    yc = f._VlsvReader__ycells
    zc = f._VlsvReader__zcells
    max_ref_level = f.get_max_refinement_level()

    cid_offsets = np.zeros(max_ref_level+1, dtype=np.int64)
    isum = 0
    for i in range(0,max_ref_level):
        isum = isum + 2**(3*i) * xc * yc * zc
        cid_offsets[i+1] = isum
    xcells = np.zeros((max_ref_level+1), dtype=np.int64)
    ycells = np.zeros((max_ref_level+1), dtype=np.int64)
    zcells = np.zeros((max_ref_level+1), dtype=np.int64)
    for r in range(max_ref_level+1):
        xcells[r] = xc*2**(r)
        ycells[r] = yc*2**(r)
        zcells[r] = zc*2**(r)
    delta = np.array([[0,0,0],[1,0,0],[0,1,0],[1,1,0],
            [0,0,1],[1,0,1],[0,1,1],[1,1,1]],dtype=np.int32)


    idxToFileIndex = {}
    from io import StringIO
    descr = StringIO()
    print("Building descriptor")
    subdivided = [[] for l in range(max_ref_level+1)]
    idx = 0
    for c in range(1,int(np.prod(f.get_spatial_mesh_size()))+1):
        if c in fileindex_for_cellid.keys():
            descr.write(".")
            idxToFileIndex[idx] = fileindex_for_cellid[c]
        else:
            descr.write("R")
            subdivided[0].append(c)
        idx += 1
    # @jit(nopython=True)
    def children(cid, level):
        # cellind = get_cell_indices(cid, level) # get_cellind here for compilation
        cellids = cid - 1 - cid_offsets[level]
        cellind = np.full(3, -1,dtype=np.int32)
        cellind[0] = (cellids)%(xcells[level])
        cellind[1] = ((cellids)//(xcells[level]))%(ycells[level])
        cellind[2] = (cellids)//(xcells[level]*ycells[level])
        cellind = cellind*2
        

        out = np.zeros((8,),dtype=np.int64)
        out[:] = cid_offsets[level+1] + (cellind[0] + delta[:,0]) + xcells[level+1]*(cellind[1] + delta[:,1]) + (cellind[2] + delta[:,2])*xcells[level+1]*ycells[level+1] + 1
        return out           
    descr.write("|")
    for l in range(1,max_ref_level+1):
        for c in subdivided[l-1]:
            for child in children(c, l-1):
                if child in fileindex_for_cellid.keys():
                    descr.write(".")
                    idxToFileIndex[idx] = fileindex_for_cellid[child]
                else:
                    descr.write("R")
                    subdivided[l].append(child)
                idx += 1
        if l < max_ref_level:
            descr.write("|")
    return descr.getvalue(), idxToFileIndex

file="/wrk-vakka/group/spacephysics/vlasiator/3D/FID/bulk1/bulk1.0001213.vlsv"
f=pt.vlsvfile.VlsvReader(file)

now=time.time()
var=buildDescriptor(f)
print(time.time()-now)

now=time.time()

f._VlsvReader__read_fileindex_for_cellid()
print(time.time()-now)
fileindex_for_cellid = f._VlsvReader__fileindex_for_cellid
xc= f._VlsvReader__xcells
yc= f._VlsvReader__ycells
zc= f._VlsvReader__zcells
max_ref_level = f.get_max_refinement_level()     
#print(f.get_spatial_mesh_size()) just [xc,yc,zc]
test=chp.test(fileindex_for_cellid,xc,yc,zc,max_ref_level)
print(time.time()-now)
print(test[0]==var[0])
print(test[1]==var[1])

