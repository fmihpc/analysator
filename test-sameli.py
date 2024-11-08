import pytools as pt


f = pt.vlsvfile.VlsvReader("C:/Users/samel/OneDrive/Tiedostot/TET/egi-like-bulk.0002702.vlsv") # vaihda tähän oikea polku tiedostoon


#f.list()


pt.plot.plot_colormap3dslice(vlsvobj=f, var='vg_b_vol', outputdir='./')

