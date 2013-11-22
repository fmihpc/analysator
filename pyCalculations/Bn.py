# File which includes calculation for B_n for bowshock

def B_n_bowshock( vlsvReader, points ):
   ''' Calculates B_n vector(s) for given points

       .. code-block::

          # Example:
          f = VlsvFile("testfile.vlsv")
          B_n_vectors = B_n_bowshock( vlsvReader=f, points=np.array([[2,3,4], [2,3,5], [5, 7, 7]])
          print B_n_vectors

       .. note:: B_n vector will be of length points-1
   '''
   


