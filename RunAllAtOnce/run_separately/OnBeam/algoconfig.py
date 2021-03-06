import sys,os

import ROOT

from larlite import larlite as fmwk

VTXSMEAR = 0.0

# ROI REMOVAL
def ROI():

    vtx_producer = 'numuCC_vertex'
    
    algo = fmwk.ROIRemoval()
    algo.setClusterProducer("pandoraCosmic")
    algo.setVertexProducer(vtx_producer)
    algo.setVerbose(False)
    algo.setVertexSmearing(VTXSMEAR)
    algo.setROI(100.)

    return algo

# VTX SLOPE CORRELATION REMOVAL

def VSC():

    vtx_producer = 'numuCC_vertex'
    
    algo = fmwk.VertexSlopeCorrelation()
    algo.setClusterProducer("pandoraCosmic")
    algo.setVertexProducer(vtx_producer)
    algo.setVerbose(False)
    algo.setCutFunction(80,-10,30,5)
    algo.setMinNHits(10)
    algo.setVertexSmearing(VTXSMEAR)
    algo.setROIRadius(100.)

    return algo

# VTX ANGLE CORRELATION REMOVAL

def VAC():

    vtx_producer = 'numuCC_vertex'

    algo = fmwk.VertexAngleCorrelation()
    algo.setClusterProducer("pandoraCosmic")
    algo.setVertexProducer(vtx_producer)
    algo.setVerbose(False)
    algo.setCutFunction(100,-5,15,5)
    algo.setVertexSmearing(VTXSMEAR)
    algo.setMaxAngle(160.)
    
    return algo


# DELTA RAY REMOVAL

def RDR():

    vtx_producer = 'numuCC_vertex'
    
    algo = fmwk.RemoveDeltaRays()
    algo.setClusterProducer("pandoraCosmic")
    algo.setVertexProducer(vtx_producer)
    algo.setVerbose(False)
    algo.setDeltaRayDistMin(1.0);
    algo.setDeltaRayDistMax(15.0);
    algo.setMaxDeltaHits(50);
    algo.setVertexSmearing(VTXSMEAR);
    algo.setROI(120.)

    return algo

def loadAlgo(algoname):

    if (algoname == "ROIRemoval"):
        return ROI()
    if (algoname == "VertexSlopeCorrelation"):
        return VSC()
    if (algoname == "VertexAngleCorrelation"):
        return VAC()
    if (algoname == "RemoveDeltaRays"):
        return RDR()

    print
    print 'ALGORITHM NOT FOUND. QUIT'
    print

    sys.exit()
