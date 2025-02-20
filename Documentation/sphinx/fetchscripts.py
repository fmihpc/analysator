import glob
import os

files = glob.glob('C:/Users/samel/OneDrive/Tiedostot/TET/analysator/scripts/*.py')

for file in files:
    file = os.path.basename(file)[:-3]
    print(file)
    with open('scripts.rst', 'a') as f:
        f.write(file + "\n")
        f.write("-"*len(file)+ "\n")
        f.write("\n" + ".. automodule:: " + (file) + "\n")
        f.write("\n" + "------------\n\n")

    # with open(file + '.rst', 'w') as f:
    #    f.write(file + "\n")
    #    f.write("-"*len(file)+ "\n")
    #    f.write("\n" + ".. automodule:: " + (file) + "\n")
    #    f.write("\t:members:")
