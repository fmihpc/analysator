#!/usr/bin/python -i
# 
# This file is part of Analysator.
# Copyright 2013-2016 Finnish Meteorological Institute
# Copyright 2017-2018 University of Helsinki
# 
# For details of usage, see the COPYING file and read the "Rules of the Road"
# at http://www.physics.helsinki.fi/vlasiator/
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# 


import visit as vis
import loadvisitsettings as visSettings
import subprocess
from useroptions import *


def make_movie( variableName, minValue, maxValue, inputDirectory, inputFileName, outputDirectory, outputFileName, colorTable="hot_desaturated", startFrame=-1, endFrame=-1 ):
   '''
   Function for making a movie
   
   Arguments:
   :param variableName                  Name of the variable
   :param minValue                      Minimum value of the variable
   :param maxValue                      Maximum value of the variable
   :param inputDirectory                Path to input vlsv/silo files
   :param inputFileName                 Name of the file(s) so for example if the filenames are bulk.0000.silo, bulk.0001.silo, .. then inputFileName=\"bulk.*.silo\""
   :param outputDirectory               Path to output directory
   :param outputFileName                Name of the output file
   :param colorTable="hot_desaturated"  Color table for the plots
   :param startFrame=-1                 Starting frame of the movie (-1 equals 0)
   :param endFrame=-1                   Starting frame of the movie (-1 equals last frame)
   '''
   # OPTIONS
   #################################################################
   # Input frame properties
   _startFrame = startFrame # Note: if _startFrame is set to -1 the start frame gets set to 0
   _endFrame = endFrame # Note: if _endFrame is set to -1 the _endFrame is automatically the number of frames in the database
   
   # Input variable
   _variableName = variableName
   minVariableValue = minValue
   maxVariableValue = maxValue
   colorTableName = colorTable
   
   # Input directory and file names
   #_outputDir = "/home/hannukse/MOVINGFRAME_MOVIES/AAJ_BZ_REMAKE/" # Set the output directory (Where .png s are saved)
   _outputDir = outputDirectory
   #_outputFileName = "BZ_FORESHOCK_2_" # The file names for the png files. These for ex. will be saved visit0000.png, visit0001.png, ..
   _outputFileName = outputFileName # The file names for the png files.
   #databaseName = "localhost:/home/hannukse/meteo/stornext/field/vlasiator/2D/AAJ/silo_files/bulk.*.silo database" # For navigating to the silo files
   databaseName = "localhost:" + inputDirectory + inputFileName + " database" # For navigating to the silo files
   # Note: a slice of the plot in z-axis is taken automatically
   #################################################################

   # LaunchNowin(vdir=visitBinDirectory)
   #dx = speedX * frameInSeconds # Note: This is in meters per frame!
   #dy = speedY * frameInSeconds # Note: This is in meters per frame!
   #LaunchNowin(vdir="/usr/local/visit/bin")
   #Set up window and annotations
   #vis.LaunchNowin(vdir="/usr/local/visit/bin")
   vis.OpenDatabase(databaseName, 0)

   #Load settings
   visSettings.load_visit_settings()

   vis.AddPlot("Pseudocolor", _variableName, 1, 1) #CONTINUE
   vis.SetActivePlots(0)
   vis.PseudocolorAtts = vis.PseudocolorAttributes()
   vis.PseudocolorAtts.legendFlag = 1
   vis.PseudocolorAtts.lightingFlag = 1
   vis.PseudocolorAtts.minFlag = 1
   vis.PseudocolorAtts.maxFlag = 1
   vis.PseudocolorAtts.centering = vis.PseudocolorAtts.Natural  # Natural, Nodal, Zonal
   vis.PseudocolorAtts.scaling = vis.PseudocolorAtts.Linear  # Linear, Log, Skew
   vis.PseudocolorAtts.limitsMode = vis.PseudocolorAtts.CurrentPlot  # OriginalData, CurrentPlot
   vis.PseudocolorAtts.min = minVariableValue
   vis.PseudocolorAtts.max = maxVariableValue
   vis.PseudocolorAtts.pointSize = 0.05
   vis.PseudocolorAtts.pointType = vis.PseudocolorAtts.Point  # Box, Axis, Icosahedron, Point, Sphere
   vis.PseudocolorAtts.skewFactor = 1
   vis.PseudocolorAtts.opacity = 1
   vis.PseudocolorAtts.colorTableName = colorTableName
   vis.PseudocolorAtts.invertColorTable = 0
   vis.PseudocolorAtts.smoothingLevel = 0
   vis.PseudocolorAtts.pointSizeVarEnabled = 0
   vis.PseudocolorAtts.pointSizeVar = "default"
   vis.PseudocolorAtts.pointSizePixels = 2
   vis.PseudocolorAtts.lineStyle = vis.PseudocolorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
   vis.PseudocolorAtts.lineWidth = 0
   vis.PseudocolorAtts.opacityType = vis.PseudocolorAtts.Explicit  # Explicit, ColorTable
   vis.SetPlotOptions(vis.PseudocolorAtts)

   
   vis.SetActivePlots(0)
   vis.AddOperator("Slice", 1)
   vis.SetActivePlots(0)
   vis.SliceAtts = vis.SliceAttributes()
   vis.SliceAtts.originType = vis.SliceAtts.Intercept  # Point, Intercept, Percent, Zone, Node
   vis.SliceAtts.originPoint = (0, 0, 0)
   vis.SliceAtts.originIntercept = 0
   vis.SliceAtts.originPercent = 0
   vis.SliceAtts.originZone = 0
   vis.SliceAtts.originNode = 0
   vis.SliceAtts.normal = (0, 0, 1)
   vis.SliceAtts.axisType = vis.SliceAtts.ZAxis  # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
   vis.SliceAtts.upAxis = (0, 1, 0)
   vis.SliceAtts.project2d = 1
   vis.SliceAtts.interactive = 1
   vis.SliceAtts.flip = 0
   vis.SliceAtts.originZoneDomain = 0
   vis.SliceAtts.originNodeDomain = 0
   vis.SliceAtts.meshName = "SpatialGrid"
   vis.SliceAtts.theta = 0
   vis.SliceAtts.phi = 90
   vis.SetOperatorOptions(vis.SliceAtts, 1)
   vis.DrawPlots()
   
   if _endFrame == -1:
      _endFrame = vis.TimeSliderGetNStates() - 1
   
   if _startFrame == -1:
      _startFrame = 0
   
   # Iterate through frames
   for i in range(_startFrame, _endFrame+1):
      vis.SetTimeSliderState(i)
      frame = i - _startFrame
      vis.SaveWindowAtts = vis.SaveWindowAttributes()
      vis.SaveWindowAtts.outputToCurrentDirectory = 0
      vis.SaveWindowAtts.outputDirectory = _outputDir
      vis.SaveWindowAtts.fileName = _outputFileName
      vis.SaveWindowAtts.family = 1
      vis.SaveWindowAtts.format = vis.SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
      vis.SaveWindowAtts.width = 3000
      vis.SaveWindowAtts.height = 3000
      vis.SaveWindowAtts.screenCapture = 0
      vis.SaveWindowAtts.saveTiled = 0
      vis.SaveWindowAtts.quality = 100
      vis.SaveWindowAtts.progressive = 0
      vis.SaveWindowAtts.binary = 0
      vis.SaveWindowAtts.stereo = 0
      vis.SaveWindowAtts.compression = vis.SaveWindowAtts.PackBits  # None, PackBits, Jpeg, Deflate
      vis.SaveWindowAtts.forceMerge = 0
      vis.SaveWindowAtts.resConstraint = vis.SaveWindowAtts.ScreenProportions  # NoConstraint, EqualWidthHeight, ScreenProportions
      vis.SaveWindowAtts.advancedMultiWindowSave = 0
      vis.SetSaveWindowAttributes(vis.SaveWindowAtts)
      vis.SaveWindow()
   vis.DeleteActivePlots()
   vis.CloseDatabase(databaseName)
   # Make the movie:
   #subprocess.call("./moviecompilescript.sh " + _outputDir + " " + _outputFileName)
   pyVisitPath = "pyVisit/"
   #subprocess.call(pythonLibDirectoryPath + pyVisitPath + "moviecompilescript.sh")
   #subprocess.call(pythonLibDirectoryPath + pyVisitPath + "moviecompilescript.sh " + _outputDir + " " + _outputFileName)
   framerate = "10"
   subprocess.call([pythonLibDirectoryPath + pyVisitPath + "moviecompilescript.sh", _outputDir, _outputFileName, framerate])
