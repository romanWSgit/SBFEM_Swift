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

enum sbfemDriverError: Error {
    case invalidShapeFunctionSelection
    case notYetImplementedInFunction
//    case insufficientFunds(coinsNeeded: Int)
//    case outOfStock
}

// MARK: - scaled boundary transformation

func rHatC(radCoordXi xi: Double, circumCoordEta eta: Double,
           coordOfElem coordArray: [Double],
           withScalingCentre scalingCentre: [Double] = [0.0, 0.0] ) -> [Double] {
    return  vectorAdd(ofVectorX: scalarTimesVector(ofScalar: xi, withVector:
        (matrixVectorProduct(ofMatrixA:
            shapeFunctionsAndDerivatives(fromEta: eta, andPolyOrd: polyOrd, calcDeriv: false).shapeMat,
                             withVectorX: coordArray))), withVectorY: scalingCentre)
}

func rHat(radCoordXi xi: Double, circumCoordEta eta: Double,
          coordOfElem coordArray: [Double]) -> [Double] {
    return  scalarTimesVector(ofScalar: xi, withVector:
        (matrixVectorProduct(ofMatrixA:
            shapeFunctionsAndDerivatives(fromEta: eta, andPolyOrd: polyOrd, calcDeriv: false).shapeMat,
                             withVectorX: coordArray)))
}

func rC(circumCoordEta eta: Double, coordOfElem coordArray: [Double],
        withScalingCentre scalingCentre: [Double] = [0.0, 0.0]) throws -> [Double] {
    if shapeFct ==  "standard shape functions" {
        return vectorAdd(ofVectorX: rHatC(radCoordXi: 1.0 , circumCoordEta: eta,
                                          coordOfElem: coordArray, withScalingCentre: scalingCentre),
                         withVectorY: scalarTimesVector(ofScalar: -1.0, withVector: scalingCentre))
    } else if shapeFct == "hierarchical shape functions" {
        print("yet to implement")
        throw sbfemDriverError.notYetImplementedInFunction
    } else {
        throw sbfemDriverError.invalidShapeFunctionSelection
    }
}

func r(circumCoordEta eta: Double, coordOfElem coordArray: [Double])  -> [Double] {
    if shapeFct ==  "standard shape functions" {
        return rHatC(radCoordXi: 1.0 , circumCoordEta: eta, coordOfElem: coordArray)
    } else if shapeFct == "hierarchical shape functions" {
        print("yet to implement")
        return rHatC(radCoordXi: 1.0 , circumCoordEta: eta, coordOfElem: coordArray)
    } else {
        print("shits gone wrong")
        return rHatC(radCoordXi: 1.0 , circumCoordEta: eta, coordOfElem: coordArray)
    }
}

func jMat(circumCoordEta eta: Double, coordOfElem coordArray: [Double],
          andPolyOrd pO: Int) -> [[Double]] {
    return [[r(circumCoordEta: eta, coordOfElem: coordArray)[0],
             r(circumCoordEta: eta, coordOfElem: coordArray)[1]],
            [(matrixVectorProduct(ofMatrixA:
                shapeFunctionsAndDerivatives(fromEta: eta, andPolyOrd: pO, calcDeriv: true).shapeMat,
                                  withVectorX: coordArray))[0],
             (matrixVectorProduct(ofMatrixA:
                shapeFunctionsAndDerivatives(fromEta: eta, andPolyOrd: pO, calcDeriv: true).shapeMat,
                                  withVectorX: coordArray))[1]]]
}

func jHatMat(radCoordXi xi: Double, circumCoordEta eta: Double,
             coordOfElem coordArray: [Double], andPolyOrd pO: Int) -> [[Double]] {
    return [[r(circumCoordEta: eta, coordOfElem: coordArray)[0],
             r(circumCoordEta: eta, coordOfElem: coordArray)[1]],
            [xi * (matrixVectorProduct(ofMatrixA:
                shapeFunctionsAndDerivatives(fromEta: eta, andPolyOrd: pO, calcDeriv: true).shapeMat,
                                       withVectorX: coordArray))[0],
             xi * (matrixVectorProduct(ofMatrixA:
                shapeFunctionsAndDerivatives(fromEta: eta, andPolyOrd: pO, calcDeriv: true).shapeMat,
                                       withVectorX: coordArray))[1]]]
}

func detJ(circumCoordEta eta: Double, coordOfElem coordArray: [Double],
          andPolyOrd pO: Int) -> Double {
    return det2x2(of2x2MatrixA: jMat(circumCoordEta: eta, coordOfElem: coordArray,
                                     andPolyOrd: pO))
}

func gXi(circumCoordEta eta: Double, coordOfElem coordArray: [Double],
         andPolyOrd pO: Int) -> [Double] {
    return [matrixVectorProduct(ofMatrixA:
        shapeFunctionsAndDerivatives(fromEta: eta, andPolyOrd: pO, calcDeriv: true).shapeMat,
                                withVectorX: coordArray)[1], -1 * matrixVectorProduct(ofMatrixA:
                                    shapeFunctionsAndDerivatives(fromEta: eta, andPolyOrd: pO, calcDeriv: true).shapeMat,
                                                                                      withVectorX: coordArray)[0]]
}

func gEta(circumCoordEta eta: Double, coordOfElem coordArray: [Double],
          andPolyOrd pO: Int) -> [Double] {
    return [-1 * r(circumCoordEta: eta, coordOfElem: coordArray)[1],
            r(circumCoordEta: eta, coordOfElem: coordArray)[0]]
}

func nXi(circumCoordEta eta: Double, coordOfElem coordArray: [Double],
         andPolyOrd pO: Int) -> [Double] {
    let vec = gXi(circumCoordEta: eta, coordOfElem: coordArray, andPolyOrd: pO)
    return  scalarTimesVector(ofScalar: (1 / sqrt(dotProduct(ofVectorX: vec, withVectorY: vec))), withVector: vec)
}

func nEta(circumCoordEta eta: Double, coordOfElem coordArray: [Double],
          andPolyOrd pO: Int) -> [Double] {
    let vec = gEta(circumCoordEta: eta, coordOfElem: coordArray, andPolyOrd: pO)
    return  scalarTimesVector(ofScalar: (1 / sqrt(dotProduct(ofVectorX: vec, withVectorY: vec))), withVector: vec)
}

// MARK: - b- and B- matrices
                                                                                
func b1(circumCoordEta eta: Double, coordOfElem coordArray: [Double], andPolyOrd pO: Int,
        withScalingCentre scalingCentre: [Double] = [0.0, 0.0]) -> [[Double]] {
    let invDet = 1 / detJ(circumCoordEta: eta, coordOfElem: coordArray, andPolyOrd: pO)
    return scalarTimesMatrix(ofScalar: invDet, withMatrix:
        [[matrixVectorProduct(ofMatrixA: shapeFunctionsAndDerivatives(fromEta: eta, andPolyOrd: pO,
                                                                      calcDeriv: true).shapeMat,
                                                                        withVectorX: coordArray)[1], 0],
         [0,-matrixVectorProduct(ofMatrixA: shapeFunctionsAndDerivatives(fromEta: eta, andPolyOrd: pO,
                                                                            calcDeriv: true).shapeMat,
                                                                            withVectorX: coordArray)[0]],
         [-matrixVectorProduct(ofMatrixA: shapeFunctionsAndDerivatives(fromEta: eta, andPolyOrd: pO,
                                                                            calcDeriv: true).shapeMat,
                                                                            withVectorX: coordArray)[0],
          matrixVectorProduct(ofMatrixA: shapeFunctionsAndDerivatives(fromEta: eta, andPolyOrd: pO,
                                                                            calcDeriv: true).shapeMat,
                                                                            withVectorX: coordArray)[1]]])
}

func b2(circumCoordEta eta: Double, coordOfElem coordArray: [Double], andPolyOrd pO: Int,
        withScalingCentre scalingCentre: [Double] = [0.0, 0.0]) -> [[Double]] {
    let invDet = 1 / detJ(circumCoordEta: eta, coordOfElem: coordArray, andPolyOrd: pO)
    return scalarTimesMatrix(ofScalar: invDet, withMatrix:
        [[-r(circumCoordEta: eta, coordOfElem: coordArray)[1], 0],
         [0, r(circumCoordEta: eta, coordOfElem: coordArray)[0]],
        [r(circumCoordEta: eta, coordOfElem: coordArray)[0],
         -r(circumCoordEta: eta, coordOfElem: coordArray)[1]]])
}

func B1(circumCoordEta eta: Double,coordOfElem coordArray: [Double], andPolyOrd pO: Int,
        withScalingCentre scalingCentre: [Double] = [0.0, 0.0]) -> [[Double]] {
    return matrixTimesMatrix(ofMatrix: b1(circumCoordEta: eta, coordOfElem: coordArray,
                                          andPolyOrd: pO),
                             withMatrix: shapeFunctionsAndDerivatives(fromEta: eta, andPolyOrd: pO,
                                                                      calcDeriv: false).shapeMat,
                             andScalarScalingFactor: 1.0)
}

func B2(circumCoordEta eta: Double,coordOfElem coordArray: [Double], andPolyOrd pO: Int,
        withScalingCentre scalingCentre: [Double] = [0.0, 0.0]) -> [[Double]] {
    return matrixTimesMatrix(ofMatrix: b2(circumCoordEta: eta, coordOfElem: coordArray,
                                          andPolyOrd: pO),
                             withMatrix: shapeFunctionsAndDerivatives(fromEta: eta, andPolyOrd: pO,
                                                                      calcDeriv: true).shapeMat,
                             andScalarScalingFactor: 1.0)
}

func zMatFct(ofCoeffMatrixE0 e0: [[Double]], andCoeffMatrixE1 e1: [[Double]], andCoeffMatrixE3 e2: [[Double]]) -> [[Double]] {
    let corrFactor = -1e-12
    let e0Inv = invertMatrix(ofMatrix: e0)
    let z12: [[Double]] = scalarTimesMatrix(ofScalar: -1.0, withMatrix: e0Inv)
    let z22: [[Double]] = matrixTimesMatrix(ofMatrix: scalarTimesMatrix(ofScalar: -1.0, withMatrix: e1), withMatrix: e0Inv, andScalarScalingFactor: 1.0)
    
    let z211: [[Double]] = scalarTimesMatrix(ofScalar: -1.0, withMatrix: e2)
    let z212: [[Double]] = matrixTimesMatrix(ofMatrix: scalarTimesMatrix(ofScalar: -1.0, withMatrix: z22), withMatrix: transposeMatrix(ofMatrix: e1), andScalarScalingFactor: 1.0)
    let z213: [[Double]] = scalarTimesMatrix(ofScalar: corrFactor, withMatrix: e0)
    let z21: [[Double]] = matrixAdd(ofMatrix: z211, withMatrix: matrixAdd(ofMatrix: z212, withMatrix: z213))
    
    let z11: [[Double]] = scalarTimesMatrix(ofScalar: -1.0, withMatrix: transposeMatrix(ofMatrix: z22))
    
    let zMat: [[Double]] = superMatrix(ofSubMat1: z11, andSubMat2: z12, andSubMat3: z21, andSubMat4: z22)
//    print()
//    print("z11:")
//    printMatrix(ofMatrix: z11)
//    print("z12:")
//    printMatrix(ofMatrix: z12)
//    print("z21:")
//    printMatrix(ofMatrix: z21)
//    print("z22:")
//    printMatrix(ofMatrix: z22)
//    print()
    return zMat
}
