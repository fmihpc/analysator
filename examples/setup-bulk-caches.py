import analysator as pt
import sys
import glob

if len(sys.argv) != 2:
   print("Need globbable argument, eg. python setup-bulk-caches.py /wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1/bulk1.*.vlsv")
   sys.exit(1)

fns = glob.glob(sys.argv[1])

fns.sort()


for fn in fns:
   f = pt.vlsvfile.VlsvReader(fn)
   f.cache_optimization_files(True)
   del f