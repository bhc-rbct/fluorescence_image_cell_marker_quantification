import qupath.ext.stardist.StarDist2D
import qupath.lib.objects.PathDetectionObject
import qupath.lib.roi.ROIs
import qupath.lib.regions.ImagePlane
import qupath.lib.io.PathIO
import qupath.lib.objects.PathObjects
import qupath.lib.common.GeneralTools
import java.nio.file.Paths
import qupath.lib.regions.ImageRegion
import qupath.opencv.ops.ImageOps
import java.math.BigDecimal
// Thresholds for marker positivity
def markerThresholds = [
    "panCK": 0.026829,
    "KRT56": 0.0770,
    //"p40": 0.1,
   // "GATA6": 0.1,
    //"HNF4A": 0.0955,
    // "aSMA": 0.1
]

def generalOutputDir = "C:/Users/Schacherer/Documents/projects/Tom/image_analysis_18_06_25/results/"
def starDistPathModel = "C:/Users/Schacherer/Documents/projects/Tom/image_analysis_18_06_25/dsb2018_heavy_augment.pb"

def newChannelNames = ["DAPI", "p40", "KLF5", "KRT56", "HNF4A", "panCK", "Bcl-xL"]

def cellExpansion = 5.0

setChannelNames(*newChannelNames)



channelMins = [
     "DAPI": 0,
     "p40": 0,
     "KLF5": 0,
     "KRT56": 0,
     "HNF4A": 0,
     "panCK": 0,
     "Bcl-xL": 0,
 ]

channelMaxs = [
     "DAPI": 255,
     "p40": 255,
     "KLF5": 255,
     "KRT56": 255,
     "HNF4A": 255,
     "panCK": 255,
     "Bcl-xL": 255,
 ]

def panCK_plus_KRT56_plus = getPathClass("KRT56+")
def panCK_plus_KRT56_minus = getPathClass("KRT56-")
def panCK_minus = getPathClass("panCK-")
def panCK_plus = getPathClass("panCK+")

def stardist = StarDist2D.builder(starDistPathModel)
    .threshold(0.5)
    .channels('DAPI')
    .normalizePercentiles(0.2, 99.8)
    .pixelSize(0.5)
    .cellExpansion(cellExpansion)
    .cellConstrainScale(1.5)
    .measureShape()
    .measureIntensity()
    .includeProbability(true)
    .build()
    
// Get current image data and project entry
def imageData = getCurrentImageData()
def entry = getProjectEntry()
def scanFileName = GeneralTools.stripExtension(entry.getImageName().split(" - resolution")[0])
def annotations = getAnnotationObjects().collect()

println "Using Image " + scanFileName

def scanResDir = generalOutputDir + scanFileName + "/"
def outputDir = new File(scanResDir, "quantifications")
outputDir.mkdirs()

//import qupath.lib.roi.ROIs
//import qupath.lib.regions.ImagePlane
// Create ROI from the image dimensions
//def plane = ImagePlane.getPlane(0, 0) // z=0, t=0
//def roi = ROIs.createRectangleROI(
//    0, 0, 
//    imageData.getServer().getWidth(), 
//    imageData.getServer().getHeight(), 
//    plane
//)

// Create annotation with the ROI
//def annotation = PathObjects.createAnnotationObject(roi)



annotations.eachWithIndex { annotation, tissueIdx ->
def tissueIdxplus1 = tissueIdx + 1
def tissueName = "Tissue_" + tissueIdxplus1
def tStartDetect = System.currentTimeMillis()
def detections = getCurrentHierarchy().getObjectsForROI(PathDetectionObject, annotation.getROI())
    removeObjects(detections, true)
println "Running for Tissue " + tissueIdxplus1
annotation.setName("T" + tissueIdxplus1)
// Run StarDist to detect nuclei within the annotation
stardist.detectObjects(imageData, [annotation])
    def tEndDetect = System.currentTimeMillis()
    println "StarDist detection took ${(tEndDetect - tStartDetect)/1000.0} seconds"
    def tStartCells = System.currentTimeMillis()
    // Get all cells (PathCellObject) within the annotation
    def cells = getCurrentHierarchy().getObjectsForROI(qupath.lib.objects.PathDetectionObject, annotation.getROI())
    def tEndCells = System.currentTimeMillis()
    println "Detected cells: " + cells.size()
    println "Cell extraction took ${(tEndCells - tStartCells)/1000.0} seconds"
    // Classify cells and collect measurements

    def allCells = []
    
    // Initialize comboCounts map
    def categoryCounts = [:]

     def tStartClass = System.currentTimeMillis()
    cells.each { cell ->
        def panCK_cell = getNormalizedValue(cell, "panCK",measurement = "Cell: Mean")
        def KRT56_cell = getNormalizedValue(cell, "KRT56", measurement = "Cell: Mean")

        if (panCK_cell > markerThresholds["panCK"]) {
            cell.setPathClass(panCK_plus)
            if (KRT56_cell > markerThresholds["KRT56"]) {
                cell.setPathClass(panCK_plus_KRT56_plus)
            } else {
                cell.setPathClass(panCK_plus_KRT56_minus)
            }
                   
        } else {
            cell.setPathClass(panCK_minus)
        }
        allCells << cell
    }
    def tEndClass = System.currentTimeMillis()
    println "Cell classification took ${(tEndClass - tStartClass)/1000.0} seconds"
    // Export per-group measurements for this tissue

    def tStartExport = System.currentTimeMillis()
    ["DAPI", "p40", "KLF5", "KRT56", "HNF4A", "panCK", "Bcl-xL"]

    // panCK+KRT56– group
    if (!allCells.isEmpty()) {
        def groupFile = new File(outputDir, tissueName + "_allCells_allMarkers.csv")
        groupFile.withWriter { writer ->
            writer.writeLine("panCK cMean,KRT56 cMean,KLF5 nMean,HNF4A nMean,Bcl-xL nMean, p40 nMean")
            allCells.each { cell ->
                writer.writeLine([
                    getNormalizedValue(cell, "panCK", measurement = "Cell: Mean"),
                    getNormalizedValue(cell, "KRT56", measurement = "Cell: Mean"),
                    getNormalizedValue(cell, "KLF5"),
                    getNormalizedValue(cell, "HNF4A"),
                    getNormalizedValue(cell, "Bcl-xL", measurement = "Cell: Mean"),
                    getNormalizedValue(cell, "p40")
                ].join(","))
            }
        }
    }
    

    def tEndExport = System.currentTimeMillis()
    println "Export took ${(tEndExport - tStartExport)/1000.0} seconds"

    
    fireHierarchyUpdate()
    cells.clear()
    System.gc() // Request garbage collection
    Thread.sleep(1000) // Pause briefly
}

stardist.close()

println "\n"


// Normalization function
/**
 * Normalizes a cell measurement value for a given marker using global min/max.
 * @param cell          The cell object
 * @param marker        Marker name (e.g., "GATA6")
 * @param measurement   Measurement suffix (e.g., "Cell: Mean")
 * @return              Normalized value (0-1), or null if measurement is missing
 */
def getNormalizedValue(cell, marker, measurement = "Nucleus: Mean") {
    def key = "${marker}: ${measurement}"
    def rawValue = cell.getMeasurementList().getMeasurementValue(key)
    if (rawValue == null) {
        print "Measurement '${key}' not found for cell ${cell}"
        return null
    }
    def minVal = channelMins[marker]
    def maxVal = channelMaxs[marker]
    if (minVal == null || maxVal == null) {
        print "Global min/max not found for marker '${marker}'"
        return null
    }
    def normValue = (rawValue - minVal) / (maxVal - minVal)
    return Math.max(0.0, Math.min(1.0, normValue)) // Clamp to [0,1]
}
