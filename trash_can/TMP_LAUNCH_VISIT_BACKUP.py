visitLaunched = False

def launch_visit(noWindow=True):
   '''
   Launches visit
   :param noWindow=True   Determines whether window is shown or not
   '''
   global visitLaunched
   if visitLaunched == True:
      print "Visit already launched"
      return
   if noWindow == True:
      LaunchNowin(vdir=pathToVisit)
   else:
      Launch(vdir=pathToVisit)
   visitLaunched = True
   print "Visit launched"
   return

