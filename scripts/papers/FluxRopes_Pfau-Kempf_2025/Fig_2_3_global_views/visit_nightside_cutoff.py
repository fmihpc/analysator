tindex = 1600

filepath = "/wrk-vakka/users/ykempf/FHA_fluxropes/fluxrope_global_views/"
db = "/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1/bulk1.000"+str(tindex)+".vlsv"

# Define expressions
DefineScalarExpression('curvature_radius_re', 'curvature_radius / 6371000.0') # noqa

TransformAtts_RE = TransformAttributes() # noqa
TransformAtts_RE.doScale = 1
TransformAtts_RE.scaleOrigin = (0, 0, 0)
TransformAtts_RE.scaleX = 1.56961e-07
TransformAtts_RE.scaleY = 1.56961e-07
TransformAtts_RE.scaleZ = 1.56961e-07



for cutoff in [3.0, 5.0, 7.0]:
   filename = "FHA_nightside_fluxropes_cutoff_"+str(tindex)+"_"+str(cutoff)+".png"
   dbvtu = "/wrk-vakka/users/ykempf/FHA_fluxropes/streamlines/streamlines_"+str(tindex)+"_"+str(cutoff)+".vtk.vtu"

   # Create plots
   # Create plot 1
   OpenDatabase(db) # noqa
   AddPlot("Contour", "vg_connection", 0, 0) # noqa
   atts = ContourAttributes() # noqa
   atts.defaultPalette.equalSpacingFlag = 1
   atts.defaultPalette.discreteFlag = 1
   atts.defaultPalette.tagNames = ("Default", "Discrete")
   atts.changedColors = (0)
   atts.colorType = atts.ColorBySingleColor  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
   atts.colorTableName = "Default"
   atts.invertColorTable = 0
   atts.legendFlag = 0
   atts.lineWidth = 0
   atts.singleColor = (192, 192, 192, 255)
   atts.contourMethod = atts.Value  # Level, Value, Percent
   atts.contourNLevels = 10
   atts.contourValue = (0.1)
   atts.contourPercent = ()
   atts.minFlag = 0
   atts.maxFlag = 0
   atts.min = 0
   atts.max = 1
   atts.scaling = atts.Linear  # Linear, Log
   atts.wireframe = 0
   SetPlotOptions(atts) # noqa

   AddOperator("Transform", 0) # noqa
   SetOperatorOptions(TransformAtts_RE, 1, 0) # noqa

#   AddOperator("Box", 1)
#   opatts = BoxAttributes()
#   opatts.amount = opatts.Some  # Some, All
#   opatts.minx = -1e9
#   opatts.maxx = 0
#   opatts.miny = -1e+09
#   opatts.maxy = 1e+09
#   opatts.minz = -1e+09
#   opatts.maxz = 1e+09
#   opatts.inverse = 0
#   SetOperatorOptions(opatts)

   # Create plot 4
   OpenDatabase(dbvtu) # noqa
   AddPlot("Pseudocolor", "curvature_radius_re", 0, 0) # noqa
   atts = PseudocolorAttributes() # noqa
   atts.scaling = atts.Linear  # Linear, Log, Skew
   atts.skewFactor = 1
   atts.limitsMode = atts.OriginalData  # OriginalData, ActualData
   atts.minFlag = 1
   atts.min = 0
   atts.useBelowMinColor = 0
   atts.belowMinColor = (0, 0, 0, 255)
   atts.maxFlag = 1
   atts.max = 6
   atts.useAboveMaxColor = 0
   atts.aboveMaxColor = (0, 0, 0, 255)
   atts.centering = atts.Natural  # Natural, Nodal, Zonal
   atts.colorTableName = "plasma"
   atts.invertColorTable = 1
   atts.opacityType = atts.FullyOpaque  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
   atts.opacityVariable = ""
   atts.opacity = 1
   atts.opacityVarMin = 0
   atts.opacityVarMax = 1
   atts.opacityVarMinFlag = 0
   atts.opacityVarMaxFlag = 0
   atts.pointSize = 0.05
   atts.pointType = atts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
   atts.pointSizeVarEnabled = 0
   atts.pointSizeVar = "default"
   atts.pointSizePixels = 2
   atts.lineType = atts.Line  # Line, Tube, Ribbon
   atts.lineWidth = 3
   atts.tubeResolution = 10
   atts.tubeRadiusSizeType = atts.FractionOfBBox  # Absolute, FractionOfBBox
   atts.tubeRadiusAbsolute = 0.125
   atts.tubeRadiusBBox = 0.005
   atts.tubeRadiusVarEnabled = 0
   atts.tubeRadiusVar = ""
   atts.tubeRadiusVarRatio = 10
   atts.tailStyle = atts.Spheres  # NONE, Spheres, Cones
   atts.headStyle = atts.NONE  # NONE, Spheres, Cones
   atts.endPointRadiusSizeType = atts.Absolute  # Absolute, FractionOfBBox
   atts.endPointRadiusAbsolute = 0.15
   atts.endPointRadiusBBox = 0.05
   atts.endPointResolution = 40
   atts.endPointRatio = 5
   atts.endPointRadiusVarEnabled = 0
   atts.endPointRadiusVar = ""
   atts.endPointRadiusVarRatio = 10
   atts.renderSurfaces = 1
   atts.renderWireframe = 0
   atts.renderPoints = 0
   atts.smoothingLevel = 0
   atts.legendFlag = 1
   atts.lightingFlag = 1
   atts.wireframeColor = (0, 0, 0, 0)
   atts.pointColor = (0, 0, 0, 0)
   SetPlotOptions(atts) # noqa

   AddOperator("Transform", 0) # noqa
   SetOperatorOptions(TransformAtts_RE, 1, 0) # noqa

   DrawPlots() # noqa

   plotName = GetPlotList().GetPlots(1).plotName # noqa
   legend = GetAnnotationObject(plotName) # noqa
   # moving the legend
   legend.managePosition = 0
   legend.position = (0.075,0.5)
   legend.xScale = 0.7
   legend.yScale = 1.82
   legend.textColor = (0, 0, 0, 255)
   legend.useForegroundForTextColor = 1
   legend.drawBoundingBox = 0
   legend.boundingBoxColor = (0, 0, 0, 50)
   legend.numberFormat = "%# -9.2g"
   legend.fontFamily = legend.Times  # Arial, Courier, Times
   legend.fontBold = 0
   legend.fontItalic = 0
   legend.fontShadow = 0
   legend.fontHeight = 0.06
   legend.drawLabels = legend.Values # None, Values, Labels, Both
   legend.drawTitle = 1
   legend.useCustomTitle = 1
   legend.customTitle = 'Curvature radius (RE)'
   legend.drawMinMax = 1
   legend.orientation = legend.HorizontalBottom  # VerticalRight, VerticalLeft, HorizontalTop, HorizontalBottom
   legend.controlTicks = 1
   legend.numTicks = 5
   legend.minMaxInclusive = 1
   legend.suppliedValues = ()
   legend.suppliedLabels = ()

   # Set the view
   view = View3DAttributes() # noqa
   view.viewNormal = (0.001, 0, 1)
   view.focus = (0, 0, 0)
   view.viewUp = (0, 1, 0)
   view.viewAngle = 30
   view.parallelScale = 30
   view.nearPlane = -40
   view.farPlane = 40
   view.imagePan = (0.35, 0)
   view.imageZoom = 1.15
   view.perspective = 0
   view.eyeAngle = 2
   view.centerOfRotationSet = 0
   view.centerOfRotation = (-50, 0, 0)
   view.axis3DScaleFlag = 0
   view.axis3DScales = (1, 1, 1)
   view.shear = (0, 0, 1)
   view.windowValid = 1
   SetView3D(view) # noqa

   # Set the annotation attributes
   AnnotationAtts = AnnotationAttributes() # noqa
   AnnotationAtts.axes3D.visible = 1
   AnnotationAtts.axes3D.autoSetTicks = 0
   AnnotationAtts.axes3D.autoSetScaling = 1
   AnnotationAtts.axes3D.lineWidth = 2
   AnnotationAtts.axes3D.tickLocation = AnnotationAtts.axes3D.Inside  # Inside, Outside, Both
   AnnotationAtts.axes3D.axesType = AnnotationAtts.axes3D.FurthestTriad  # ClosestTriad, FurthestTriad, OutsideEdges, StaticTriad, StaticEdges
   AnnotationAtts.axes3D.triadFlag = 1
   AnnotationAtts.axes3D.bboxFlag = 1
   AnnotationAtts.axes3D.xAxis.title.visible = 0
   AnnotationAtts.axes3D.xAxis.title.font.font = AnnotationAtts.axes3D.xAxis.title.font.Times  # Arial, Courier, Times
   AnnotationAtts.axes3D.xAxis.title.font.scale = 0.5
   AnnotationAtts.axes3D.xAxis.title.font.useForegroundColor = 1
   AnnotationAtts.axes3D.xAxis.title.font.color = (0, 0, 0, 255)
   AnnotationAtts.axes3D.xAxis.title.font.bold = 0
   AnnotationAtts.axes3D.xAxis.title.font.italic = 0
   AnnotationAtts.axes3D.xAxis.title.userTitle = 1
   AnnotationAtts.axes3D.xAxis.title.userUnits = 1
   AnnotationAtts.axes3D.xAxis.title.title = "X"
   AnnotationAtts.axes3D.xAxis.title.units = "RE"
   AnnotationAtts.axes3D.xAxis.label.visible = 1
   AnnotationAtts.axes3D.xAxis.label.font.font = AnnotationAtts.axes3D.xAxis.label.font.Times  # Arial, Courier, Times
   AnnotationAtts.axes3D.xAxis.label.font.scale = 1
   AnnotationAtts.axes3D.xAxis.label.font.useForegroundColor = 1
   AnnotationAtts.axes3D.xAxis.label.font.color = (0, 0, 0, 255)
   AnnotationAtts.axes3D.xAxis.label.font.bold = 0
   AnnotationAtts.axes3D.xAxis.label.font.italic = 0
   AnnotationAtts.axes3D.xAxis.label.scaling = 0
   AnnotationAtts.axes3D.xAxis.tickMarks.visible = 1
   AnnotationAtts.axes3D.xAxis.tickMarks.majorMinimum = -150
   AnnotationAtts.axes3D.xAxis.tickMarks.majorMaximum = 100
   AnnotationAtts.axes3D.xAxis.tickMarks.minorSpacing = 2
   AnnotationAtts.axes3D.xAxis.tickMarks.majorSpacing = 10
   AnnotationAtts.axes3D.xAxis.grid = 0
   AnnotationAtts.axes3D.yAxis.title.visible = 0
   AnnotationAtts.axes3D.yAxis.title.font.font = AnnotationAtts.axes3D.yAxis.title.font.Times  # Arial, Courier, Times
   AnnotationAtts.axes3D.yAxis.title.font.scale = 0.5
   AnnotationAtts.axes3D.yAxis.title.font.useForegroundColor = 1
   AnnotationAtts.axes3D.yAxis.title.font.color = (0, 0, 0, 255)
   AnnotationAtts.axes3D.yAxis.title.font.bold = 0
   AnnotationAtts.axes3D.yAxis.title.font.italic = 0
   AnnotationAtts.axes3D.yAxis.title.userTitle = 1
   AnnotationAtts.axes3D.yAxis.title.userUnits = 1
   AnnotationAtts.axes3D.yAxis.title.title = "Y"
   AnnotationAtts.axes3D.yAxis.title.units = "RE"
   AnnotationAtts.axes3D.yAxis.label.visible = 1
   AnnotationAtts.axes3D.yAxis.label.font.font = AnnotationAtts.axes3D.yAxis.label.font.Times  # Arial, Courier, Times
   AnnotationAtts.axes3D.yAxis.label.font.scale = 1
   AnnotationAtts.axes3D.yAxis.label.font.useForegroundColor = 1
   AnnotationAtts.axes3D.yAxis.label.font.color = (0, 0, 0, 255)
   AnnotationAtts.axes3D.yAxis.label.font.bold = 0
   AnnotationAtts.axes3D.yAxis.label.font.italic = 0
   AnnotationAtts.axes3D.yAxis.label.scaling = 0
   AnnotationAtts.axes3D.yAxis.tickMarks.visible = 1
   AnnotationAtts.axes3D.yAxis.tickMarks.majorMinimum = -150
   AnnotationAtts.axes3D.yAxis.tickMarks.majorMaximum = 100
   AnnotationAtts.axes3D.yAxis.tickMarks.minorSpacing = 2
   AnnotationAtts.axes3D.yAxis.tickMarks.majorSpacing = 10
   AnnotationAtts.axes3D.yAxis.grid = 0
   AnnotationAtts.axes3D.zAxis.title.visible = 0
   AnnotationAtts.axes3D.zAxis.title.font.font = AnnotationAtts.axes3D.zAxis.title.font.Times  # Arial, Courier, Times
   AnnotationAtts.axes3D.zAxis.title.font.scale = 1
   AnnotationAtts.axes3D.zAxis.title.font.useForegroundColor = 1
   AnnotationAtts.axes3D.zAxis.title.font.color = (0, 0, 0, 255)
   AnnotationAtts.axes3D.zAxis.title.font.bold = 0
   AnnotationAtts.axes3D.zAxis.title.font.italic = 0
   AnnotationAtts.axes3D.zAxis.title.userTitle = 0
   AnnotationAtts.axes3D.zAxis.title.userUnits = 0
   AnnotationAtts.axes3D.zAxis.title.title = "Z"
   AnnotationAtts.axes3D.zAxis.title.units = "RE"
   AnnotationAtts.axes3D.zAxis.label.visible = 0
   AnnotationAtts.axes3D.zAxis.label.font.font = AnnotationAtts.axes3D.zAxis.label.font.Times  # Arial, Courier, Times
   AnnotationAtts.axes3D.zAxis.label.font.scale = 1
   AnnotationAtts.axes3D.zAxis.label.font.useForegroundColor = 1
   AnnotationAtts.axes3D.zAxis.label.font.color = (0, 0, 0, 255)
   AnnotationAtts.axes3D.zAxis.label.font.bold = 0
   AnnotationAtts.axes3D.zAxis.label.font.italic = 0
   AnnotationAtts.axes3D.zAxis.label.scaling = 0
   AnnotationAtts.axes3D.zAxis.tickMarks.visible = 0
   AnnotationAtts.axes3D.zAxis.tickMarks.majorMinimum = 0
   AnnotationAtts.axes3D.zAxis.tickMarks.majorMaximum = 1
   AnnotationAtts.axes3D.zAxis.tickMarks.minorSpacing = 0.02
   AnnotationAtts.axes3D.zAxis.tickMarks.majorSpacing = 0.2
   AnnotationAtts.axes3D.zAxis.grid = 0
   AnnotationAtts.axes3D.setBBoxLocation = 1
   AnnotationAtts.axes3D.bboxLocation = (-110, 11, -23, 21, -30, 30)
   AnnotationAtts.axes3D.triadColor = (0, 0, 0)
   AnnotationAtts.axes3D.triadLineWidth = 3
   AnnotationAtts.axes3D.triadFont = AnnotationAtts.axes3D.zAxis.label.font.Times
   AnnotationAtts.axes3D.triadBold = 1
   AnnotationAtts.axes3D.triadItalic = 0
   AnnotationAtts.axes3D.triadSetManually = 0
   AnnotationAtts.userInfoFlag = 0
   AnnotationAtts.userInfoFont.font = AnnotationAtts.userInfoFont.Arial  # Arial, Courier, Times
   AnnotationAtts.userInfoFont.scale = 1
   AnnotationAtts.userInfoFont.useForegroundColor = 1
   AnnotationAtts.userInfoFont.color = (0, 0, 0, 255)
   AnnotationAtts.userInfoFont.bold = 0
   AnnotationAtts.userInfoFont.italic = 0
   AnnotationAtts.databaseInfoFlag = 0
   AnnotationAtts.timeInfoFlag = 1
   AnnotationAtts.databaseInfoFont.font = AnnotationAtts.databaseInfoFont.Arial  # Arial, Courier, Times
   AnnotationAtts.databaseInfoFont.scale = 1
   AnnotationAtts.databaseInfoFont.useForegroundColor = 1
   AnnotationAtts.databaseInfoFont.color = (0, 0, 0, 255)
   AnnotationAtts.databaseInfoFont.bold = 0
   AnnotationAtts.databaseInfoFont.italic = 0
   AnnotationAtts.databaseInfoExpansionMode = AnnotationAtts.File  # File, Directory, Full, Smart, SmartDirectory
   AnnotationAtts.databaseInfoTimeScale = 1
   AnnotationAtts.databaseInfoTimeOffset = 0
   AnnotationAtts.legendInfoFlag = 1
   AnnotationAtts.backgroundColor = (255, 255, 255, 255)
   AnnotationAtts.foregroundColor = (0, 0, 0, 255)
   AnnotationAtts.gradientBackgroundStyle = AnnotationAtts.Radial  # TopToBottom, BottomToTop, LeftToRight, RightToLeft, Radial
   AnnotationAtts.gradientColor1 = (0, 0, 255, 255)
   AnnotationAtts.gradientColor2 = (0, 0, 0, 255)
   AnnotationAtts.backgroundMode = AnnotationAtts.Solid  # Solid, Gradient, Image, ImageSphere
   AnnotationAtts.backgroundImage = ""
   AnnotationAtts.imageRepeatX = 1
   AnnotationAtts.imageRepeatY = 1
   AnnotationAtts.axesArray.visible = 1
   AnnotationAtts.axesArray.ticksVisible = 1
   AnnotationAtts.axesArray.autoSetTicks = 1
   AnnotationAtts.axesArray.autoSetScaling = 1
   AnnotationAtts.axesArray.lineWidth = 0
   AnnotationAtts.axesArray.axes.title.visible = 1
   AnnotationAtts.axesArray.axes.title.font.font = AnnotationAtts.axesArray.axes.title.font.Arial  # Arial, Courier, Times
   AnnotationAtts.axesArray.axes.title.font.scale = 1
   AnnotationAtts.axesArray.axes.title.font.useForegroundColor = 1
   AnnotationAtts.axesArray.axes.title.font.color = (0, 0, 0, 255)
   AnnotationAtts.axesArray.axes.title.font.bold = 0
   AnnotationAtts.axesArray.axes.title.font.italic = 0
   AnnotationAtts.axesArray.axes.title.userTitle = 0
   AnnotationAtts.axesArray.axes.title.userUnits = 0
   AnnotationAtts.axesArray.axes.title.title = ""
   AnnotationAtts.axesArray.axes.title.units = ""
   AnnotationAtts.axesArray.axes.label.visible = 1
   AnnotationAtts.axesArray.axes.label.font.font = AnnotationAtts.axesArray.axes.label.font.Arial  # Arial, Courier, Times
   AnnotationAtts.axesArray.axes.label.font.scale = 1
   AnnotationAtts.axesArray.axes.label.font.useForegroundColor = 1
   AnnotationAtts.axesArray.axes.label.font.color = (0, 0, 0, 255)
   AnnotationAtts.axesArray.axes.label.font.bold = 0
   AnnotationAtts.axesArray.axes.label.font.italic = 0
   AnnotationAtts.axesArray.axes.label.scaling = 0
   AnnotationAtts.axesArray.axes.tickMarks.visible = 1
   AnnotationAtts.axesArray.axes.tickMarks.majorMinimum = 0
   AnnotationAtts.axesArray.axes.tickMarks.majorMaximum = 1
   AnnotationAtts.axesArray.axes.tickMarks.minorSpacing = 0.02
   AnnotationAtts.axesArray.axes.tickMarks.majorSpacing = 0.2
   AnnotationAtts.axesArray.axes.grid = 0
   SetAnnotationAttributes(AnnotationAtts) # noqa

   SaveWindowatts = SaveWindowAttributes() # noqa
   SaveWindowatts.outputToCurrentDirectory = 0
   SaveWindowatts.outputDirectory = filepath
   SaveWindowatts.fileName = filename
   SaveWindowatts.family = 0
   SaveWindowatts.format = SaveWindowatts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY, EXR
   SaveWindowatts.width = 3200
   SaveWindowatts.height = 1300
   SaveWindowatts.screenCapture = 0
   SaveWindowatts.saveTiled = 0
   SaveWindowatts.quality = 80
   SaveWindowatts.progressive = 0
   SaveWindowatts.binary = 0
   SaveWindowatts.stereo = 0
#   SaveWindowatts.compression = SaveWindowatts.NONE  # NONE, PackBits, Jpeg, Deflate, LZW
   SaveWindowatts.forceMerge = 0
   SaveWindowatts.resConstraint = SaveWindowatts.NoConstraint  # NoConstraint, EqualWidthHeight, ScreenProportions
   SaveWindowatts.pixelData = 1
   SaveWindowatts.opts.types = ()
   SaveWindowatts.opts.help = ""
   SetSaveWindowAttributes(SaveWindowatts) # noqa
   SaveWindow() # noqa
   
   DeleteAllPlots() # noqa
   CloseDatabase(db) # noqa
   CloseDatabase(dbvtu) # noqa


exit()
