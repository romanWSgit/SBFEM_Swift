//
//  sbfemDriver.swift
//  SBFEM
//
//  Created by Roman Wallner- Silberhuber on 03.04.19.
//  Copyright © 2019 Roman Wallner- Silberhuber. All rights reserved.
//

import Foundation
import Accelerate
import simd
import Surge


func rHatC(circumferentialCoordianteEta eta: Double, radialCoordinateXi xi: Double, coordianteOfElementsPoints coordArray: [Double], withScalingCentre scalingCentre: [Double] ) -> [Double] {
    let rHatCVec =  scalarTimesVector(ofScalar: xi, withVector: matrixVectorProduct(ofMatrixA: shapeFunctionsAndDerivatives(fromEta: eta, andPolyOrd: polyOrd, calcDeriv: false).shapeMat, withVectorX: coordArray, plusScaling: 1.0, timesVector: scalingCentre))
    return rHatCVec
}
