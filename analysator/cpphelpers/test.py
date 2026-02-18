import cpphelpers as chp
import analysator as pt
import time

file="/home/siclasse/bulk.0000110.vlsv"
f=pt.vlsvfile.VlsvReader(file)
test=pt.vlsvfile.VlsvVtkReader()
test.SetFileName(file)

var=test.buildDescriptorPython()

now=time.time()
var=test.buildDescriptorPython()
print(time.time()-now)

now=time.time()
test=test.buildDescriptor()
print(time.time()-now)
print(test[0]==var[0])
print(test[1]==var[1])

