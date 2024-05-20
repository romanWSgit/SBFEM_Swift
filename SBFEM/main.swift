//
//  main.swift
//  SBFEM
//
//  Created by Roman Wallner- Silberhuber on 27.03.19.
//  Copyright © 2019 Roman Wallner- Silberhuber. All rights reserved.
//

// MARK: - IMPORT SECTION
import Foundation
import Accelerate
import simd
import Darwin

// MARK: - CLASS STRUCTURE

/// The main SBFEM model class
class model {
    
    /// node class
    class node {
        
    }
    
    /// element class
    class element {
        
    }
}





// MARK: - MATH TESTS
print("------------------")
print("---MATH TEST------")
mathTest()
print("------------------")
// MARK: - PARAMS
/// Machine Double Precission Epsilon
print("------------------")
let eps64 = Double.ulpOfOne
print("EPS 64 = \(eps64)")
print("------------------")


/// spatial Dimension
let sd = 2
/// polynominal Order
let polyOrd = 1
/// number of elements
let ielem = 10
/// number of nodes per element 
let nodeElem = 2
/// number of all nodes in the domain
let nodeDim = 10
/// number of elements in 'a' - direction (horizontal)
let dim = nodeDim * sd
/// Dimesion
let a = 4
/// number of elements in 'b' - direction (vertical)
let b = 1
// lengths in 'a' - direction (horizontal)
let l = 24.0
// width in 'b' - direction (horizontal)
let d = 6.0

let integrationOrder = 2

let offset: [Double] = [0,0]
let centre: [Double] = [0,0]
let shapeFct = "standard shape functions"
// Material
let Em: Double = 10000
let nu: Double = 0.2
let rho: Double = 1
let Er: Double = Em
let rr: Double = 1
let rhor: Double = rho

// MARK: - LOADS

let fECantilever: [Double] = [0,-6]




//func printMatrix( _ matrix:[[Double]] ) {
//    for array in matrix {
//        print( array )
//    }
//}

// MARK: - RUN


//printMatrix(ofMatrix: b1(circumCoordEta: 0.5, coordOfElem: [-12.0,  -3.0,  -6.0,  -3.0], andPolyOrd: 1))
//print("\n")
//printMatrix(ofMatrix: b2(circumCoordEta: 0.5, coordOfElem: [-12.0,  -3.0,  -6.0,  -3.0], andPolyOrd: 1))
//print("\n")
//printMatrix(ofMatrix: B1(circumCoordEta: 0.5, coordOfElem: [-12.0,  -3.0,  -6.0,  -3.0], andPolyOrd: 1))
//print("\n")
//printMatrix(ofMatrix: B2(circumCoordEta: 0.5, coordOfElem: [-12.0,  -3.0,  -6.0,  -3.0], andPolyOrd: 1))
//print("\n")
//printMatrix(ofMatrix: invertMatrix(ofMatrix: [[1,2],[3,4]]))
//print("\n")
//printMatrix(ofMatrix: flatDiskMesh(forNumberOfElements: ielem, withPolyOrd: polyOrd, withNodeDimension: nodeDim, withHorizontalDivision: a, withVerticalDivision: b, withHeightD: d, withLengthL: l).elementTable)


print(getRootsAndWeights(forIntegrationOrder: 2).Roots)

print(fCondListCantilever(forSbfemElement: 5, withPolyOrd: 1, withHorizontalDivision: a, withVerticalDivision: b))
var fECantileverMat: [[Double]] = [fECantilever]
print("\n")
print(transposeMatrix(ofMatrix: matrixTimesMatrix(ofMatrix: transposeMatrix(ofMatrix: shapeFunctionsAndDerivatives(fromEta: 0.5, andPolyOrd: polyOrd, calcDeriv: false).shapeMat) , withMatrix: transposeMatrix(ofMatrix: fECantileverMat), andScalarScalingFactor: 1.0))[0])

print("\n")
print("start")
let LTG = flatDiskGeom(forNumberOfElements: ielem, withPolyOrd: polyOrd).LTG
let xs = getRootsAndWeights(forIntegrationOrder: integrationOrder).Roots
let ws = getRootsAndWeights(forIntegrationOrder: integrationOrder).Weights
let elements = flatDiskMesh(forNumberOfElements: ielem, withPolyOrd: polyOrd, withNodeDimension: nodeDim, withHorizontalDivision: a, withVerticalDivision: b, withHeightD: d, withLengthL: l).elementTable
let nd = nodeDim * sd
// MARK: - Math Test

let testVectorAdd = vectorAdd(ofVectorX: [1, 2, 3], withVectorY: [4, 5, -1] )
print(testVectorAdd)
// seems to work

let testScalarTimesVector = scalarTimesVector(ofScalar: 9, withVector: [4, 5, -1])
print(testScalarTimesVector)
// seems to works

let testScalarTimesMatrix = scalarTimesMatrix(ofScalar: 2, withMatrix: [[1, 2], [3, 4]])
print(testScalarTimesMatrix)
// seems to work

let testDotProduct = dotProduct(ofVectorX: [1, 2, -3], withVectorY: [2, 2, 1])
print(testDotProduct)
// seems to work

let testDyadicP = dyadicP(withVectorX:  [1, 2, -3], andVectorY: [2, 2, 1])
print(testDyadicP)
// seems to work

let testMatrixVectorProduct1 = matrixVectorProduct(ofMatrixA: [[1, 2], [3, 4]], withVectorX: [1, 2])
print(testMatrixVectorProduct1)
// seems to work

let testMatrixVectorProduct2 = matrixVectorProduct(ofMatrixA: [[1, 2], [3, 4]], withVectorX: [1, 2], plusScaling: 0.5, timesVector: [2, 3])
print(testMatrixVectorProduct2)
// seems to work

let testMatrixTimesMatrix = matrixTimesMatrix(ofMatrix: [[1, 2], [3, 4]], withMatrix: [[3, 2], [3, 4]], andScalarScalingFactor: 0.5)
print(testMatrixTimesMatrix)
// seems to work


let testMatrixTimesMatrixPlusC = matrixTimesMatrixPlusC(ofMatrix: [[1, 2], [3, 4]], withMatrix: [[3, 2], [3, 4]], addMatrixC: [[1, 2], [1, 2]], andScalarScalingFactor: 2)
print(testMatrixTimesMatrixPlusC)



// MARK: - Coefficient Matrices Assembly
var m0: [[Double]] = Array(repeating: Array(repeating: 0, count: nd), count: nd)
var e0: [[Double]] = Array(repeating: Array(repeating: 0, count: nd), count: nd)
var e1: [[Double]] = Array(repeating: Array(repeating: 0, count: nd), count: nd)
var e2: [[Double]] = Array(repeating: Array(repeating: 0, count: nd), count: nd)
var f: [Double] = Array(repeating: 0, count: nd)
for el in 0..<ielem {
    let ndel = sd*(polyOrd+1)
    var m0el: [[Double]] = Array(repeating: Array(repeating: 0, count: ndel), count: ndel)
    var e0el: [[Double]] = Array(repeating: Array(repeating: 0, count: ndel), count: ndel)
    var e1el: [[Double]] = Array(repeating: Array(repeating: 0, count: ndel), count: ndel)
    var e2el: [[Double]] = Array(repeating: Array(repeating: 0, count: ndel), count: ndel)
    var fel: [Double] = Array(repeating: 0, count: ndel)

    for j in 0..<xs.count {
        let coordVec = elements[el]
        let eta = xs[j]
        let detJScalar = detJ(circumCoordEta: eta, coordOfElem: coordVec, andPolyOrd: polyOrd)
        let gXiVec = gXi(circumCoordEta: eta, coordOfElem: coordVec, andPolyOrd: polyOrd)
        let B1Mat = B1(circumCoordEta: eta, coordOfElem: coordVec, andPolyOrd: polyOrd)
        let B2Mat = B2(circumCoordEta: eta, coordOfElem: coordVec, andPolyOrd: polyOrd)
        let NMat = shapeFunctionsAndDerivatives(fromEta: eta, andPolyOrd: polyOrd, calcDeriv: false).shapeMat
        let scale = 1/Er * ws[j] * detJScalar
        let fECantileverMat: [[Double]] = [fECantilever]

        e0el = matrixTimesMatrixPlusC(ofMatrix: transposeMatrix(ofMatrix: B1Mat),withMatrix:
            matrixTimesMatrix(ofMatrix: dMat, withMatrix: B1Mat, andScalarScalingFactor: 1.0),
                                      addMatrixC: e0el, andScalarScalingFactor: scale)
        e1el = matrixTimesMatrixPlusC(ofMatrix: transposeMatrix(ofMatrix: B2Mat), withMatrix:
            matrixTimesMatrix(ofMatrix: dMat, withMatrix: B1Mat, andScalarScalingFactor: 1.0),
                                      addMatrixC: e1el, andScalarScalingFactor: scale)
        e2el = matrixTimesMatrixPlusC(ofMatrix: transposeMatrix(ofMatrix: B2Mat), withMatrix:
            matrixTimesMatrix(ofMatrix: dMat, withMatrix: B2Mat, andScalarScalingFactor: 1.0),
                                      addMatrixC: e2el, andScalarScalingFactor: scale)
        m0el = matrixTimesMatrixPlusC(ofMatrix: transposeMatrix(ofMatrix: NMat), withMatrix:
            NMat, addMatrixC: m0el, andScalarScalingFactor: scale * rho)
        if fCondListCantilever(forSbfemElement: el, withPolyOrd: polyOrd,
                               withHorizontalDivision: a, withVerticalDivision: b) {
            fel = vectorAdd(ofVectorX: fel, withVectorY:
                transposeMatrix(ofMatrix: matrixTimesMatrix(ofMatrix:
                    transposeMatrix(ofMatrix: NMat) ,withMatrix:
                    transposeMatrix(ofMatrix: fECantileverMat), andScalarScalingFactor:
                    1/Er * ws[j] * sqrt(dotProduct(ofVectorX: gXiVec, withVectorY: gXiVec))))[0])
        }
    }
    let dof = (LTG[el].map({$0 - 1}))
    e0 = addAtIx(inMatrix: e0, insertMatrix: e0el, atIndex: dof)
    e1 = addAtIx(inMatrix: e1, insertMatrix: e1el, atIndex: dof)
    e2 = addAtIx(inMatrix: e2, insertMatrix: e2el, atIndex: dof)
    m0 = addAtIx(inMatrix: m0, insertMatrix: m0el, atIndex: dof)
    f = addAtIxVector(inVector: f, insertVector: fel, atIndex: dof)
}
print(LTG[9].map({$0 - 1}))
//print("")
//print("")
//print("e0el")
//print(e0el)
//print("")
//print("")

// MARK: - ZMat creation

let zMat = zMatFct(ofCoeffMatrixE0: e0, andCoeffMatrixE1: e1, andCoeffMatrixE3: e2)

//print("")
//print("")
//print("zMat")
//printMatrix(ofMatrix: zMat)
//print("")
//print("")
//printMatrixDimensions(ofMatrix: zMat)
//print("")
//print("")
if #available(macOS 13.3, *) {
    let (valRe, valIm ,vec2) = eigensystem(ofMatrix: zMat)
} else {
    // Fallback on earlier versions
}
//let valIm = eigensystem(ofMatrix: zMat).eigenvaluesIm
//let vec = transposeMatrix(ofMatrix: eigensystem(ofMatrix: zMat).rightEigenvectors)
//let vec = eigensystem(ofMatrix: zMat).rightEigenvectors
//let vec, = eigensystem(ofMatrix: transposeMatrix(ofMatrix: zMat)).rightEigenvectors
let vec = transposeMatrix(ofMatrix: vec2)
print("")
print("")
print("Eigenvalues")
print("")
print("")
var valSorted = valRe
valSorted.sort()
print(valSorted)
//printMatrix(e1)
print(valRe)
print("EigenvaluesIM")
print(valIm)


//var (Uvec, s, VTvec) = singularValueDecomposition(ofMatrix: zMat)
//print("")
//print("")
//print("Uvec")
//print(Uvec)
//print("")
//print("")
//print("")
//print("")
//print("s")
//print(s)
//print("")
//print("")
//print("")
//print("")
//print("VTvec")
//print(VTvec)
//print("")
//print("")
//let index = vec.firstIndex(of: vecSorted[0])
//print(index ?? "not found")

var indexList = [Int]()


for i in 0..<valRe.count {
    let index = valRe.firstIndex(of: valSorted[i])
    indexList.append(index!)
}
print("indexList: \(indexList)")

var indexListNew = indexList
for i in 0..<indexListNew.count {
    if (i < (indexListNew.count - 1)) && indexListNew[i+1] == indexListNew[i] {
        indexListNew[i+1] = indexListNew[i+1] + 1
    }
}
print("indexListNew: \(indexListNew)")
var vecSorted = [Double]()
var vecUnSorted = Array(vec.joined())
for i in 0..<(dim * 2) {
    let index = indexListNew[i]
    let vecPart = Array(vecUnSorted[(dim * 2 * index)..<((dim * 2)  * (index+1))])
    for j in 0..<(dim * 2){
        vecSorted.append(vecPart[j])
    }
}

//var vecSortedAlternative = [DSPDoubleComplex]()
//var vecUnSortedAlternative = Array(vec.joined())
//var iV = 0
//while iV < (dim * 2) {
//    let index = indexList[iV]
//
//    if iV < (dim*2 - 1) {
//        if  index == indexList[iV+1]   {
//            var vecPartComplexAtPosI = Array(vecUnSortedAlternative[(dim * 2 * index)..<((dim * 2)  * (index+1))])
//            var vecPartComplexAtPosIPlus1 = Array(vecUnSortedAlternative[(dim * 2 * (index + 1))..<((dim * 2)  * (index + 2))])
//            var vecPartComplexAtPosIPlus1Minus = scalarTimesVector(ofScalar: -1.0, withVector:
//                Array(vecUnSorted[(dim * 2 * (index + 1))..<((dim * 2)  * (index + 2))]))
//            var splitComplexAtPosI = DSPDoubleSplitComplex(realp: &vecPartComplexAtPosI,
//                                                           imagp: &vecPartComplexAtPosIPlus1 )
//            var splitComplexAtPosIPlus1 = DSPDoubleSplitComplex(realp: &vecPartComplexAtPosI,
//                                                                imagp: &vecPartComplexAtPosIPlus1Minus  )
//            
//            var complex = [DSPDoubleComplex](repeating: DSPDoubleComplex(),
//                                             count: vecPartComplexAtPosI .count)
//            var complex2 = [DSPDoubleComplex](repeating: DSPDoubleComplex(),
//                                              count: vecPartComplexAtPosI .count)
//            let n = vDSP_Length(vecPartComplexAtPosI.count)
//            
//            vDSP_ztocD(&splitComplexAtPosI, 1,
//                       &complex, 2,
//                       n)
//            vDSP_ztocD(&splitComplexAtPosIPlus1, 1,
//                       &complex2, 2,
//                       n)
//            for j in 0..<(dim * 2){
//                vecSortedAlternative.append(complex[j])
//            }
//            for j in 0..<(dim * 2){
//                vecSortedAlternative.append(complex2[j])
//            }
//            iV = iV + 2
//        } else {
//            var vecPartComplexAtPosI = Array(vecUnSortedAlternative[(dim * 2 * index)..<((dim * 2)  * (index+1))])
//            var vecPartComplexAtPosIPlus1 = Array(repeating: 0.0, count: dim*2)
//            var splitComplexAtPosI = DSPDoubleSplitComplex(realp: &vecPartComplexAtPosI,
//                                                           imagp: &vecPartComplexAtPosIPlus1)
//            var complex = [DSPDoubleComplex](repeating: DSPDoubleComplex(),
//                                             count: vecPartComplexAtPosI .count)
//            let n = vDSP_Length(vecPartComplexAtPosI.count)
//            
//            vDSP_ztocD(&splitComplexAtPosI, 1,
//                       &complex, 2,
//                       n)
//            for j in 0..<(dim * 2){
//                vecSortedAlternative.append(complex[j])
//            }
//            iV += 1
//        }
//    }
//    var vecPartComplexAtPosI = Array(vecUnSortedAlternative[(dim * 2 * index)..<((dim * 2)  * (index+1))])
//    var vecPartComplexAtPosIPlus1 = Array(repeating: 0.0, count: dim*2)
//    var splitComplexAtPosI = DSPDoubleSplitComplex(realp: &vecPartComplexAtPosI,
//                                                   imagp: &vecPartComplexAtPosIPlus1)
//    var complex = [DSPDoubleComplex](repeating: DSPDoubleComplex(),
//                                     count: vecPartComplexAtPosI .count)
//    let n = vDSP_Length(vecPartComplexAtPosI.count)
//    
//    vDSP_ztocD(&splitComplexAtPosI, 1,
//               &complex, 2,
//               n)
//    for j in 0..<(dim * 2){
//        vecSortedAlternative.append(complex[j])
//    }
//    iV += 1
//   
//}
//
//print(indexListNew)
//print("Eigenvec:")
//printMatrix(ofMatrix: vec)
//print("indexListNewCount:")
//print(indexListNew.count)
//print("indexListNew:")
//
//print()
//let arrSorted: [[Double]] = stride(from: 0, to: dim*4*dim, by: dim*2).map {
//    Array(vecSorted[$0..<$0+dim*2])}
//print("arrSorted :")
//print(arrSorted)
//print()
//printMatrixDimensions(ofMatrix: arrSorted)
//print()
//
//print()
//let arrSortedComplex: [[DSPDoubleComplex]] = stride(from: 0, to: dim*4*dim, by: dim*2).map {
//    Array(vecSortedAlternative[$0..<$0+dim*2])}
//
//// MARK: - LOADS
//
//
//let psiU1 = arrSorted[0..<dim].map { $0[0..<dim].compactMap { $0 } }
//let psiU2 = arrSorted[0..<dim].map { $0[dim..<dim*2].compactMap { $0 } }
//let psiQ1 = arrSorted[dim..<dim*2].map { $0[0..<dim].compactMap { $0 } }
//let psiQ2 = arrSorted[dim..<dim*2].map { $0[dim..<dim*2].compactMap { $0 } }
//let matK = matrixTimesMatrix(ofMatrix:  psiQ1 , withMatrix: invertMatrix(ofMatrix:  psiU1)   , andScalarScalingFactor: 1.0)
//let psiU1Vec = Array(psiU1.joined())
//
//let psiU1Com = arrSortedComplex[0..<dim].map { $0[0..<dim].compactMap { $0 } }
//let psiU2Com = arrSortedComplex[0..<dim].map { $0[dim..<dim*2].compactMap { $0 } }
//let psiQ1Com = arrSortedComplex[dim..<dim*2].map { $0[0..<dim].compactMap { $0 } }
//let psiQ2Com = arrSortedComplex[dim..<dim*2].map { $0[dim..<dim*2].compactMap { $0 } }
//let psiU1ComVec = Array(psiU1Com.joined())
//let psiU2ComVec = Array(psiU2Com.joined())
//let psiQ1ComVec = Array(psiQ1Com.joined())
//let psiQ2ComVec = Array(psiQ2Com.joined())
//print(psiU1ComVec[1].real)
//var psiU1ComVecReal = [Double]()
//var psiU1ComVecImag = [Double]()
//for i in 0..<psiU1ComVec.count  {
//    let real: Double = psiU1ComVec[i].real
//    let imag: Double = psiU1ComVec[i].imag
//    psiU1ComVecReal.append(real)
//    psiU1ComVecImag.append(imag)
//}
//let dU1 = psiU1Com.count
//let psiU1ComReal: [[Double]] = stride(from: 0, to: dU1*dU1, by: dU1).map {
//    Array(psiU1ComVecReal[$0..<$0+dU1])}
//let psiU1ComImag: [[Double]] = stride(from: 0, to: dU1*dU1, by: dU1).map {
//    Array(psiU1ComVecImag[$0..<$0+dU1])}
//print("shit:")
//printMatrix(ofMatrix: invertMatrix(ofMatrix: psiU1ComReal))
//print("shit:")
//printMatrix(ofMatrix: invertMatrix(ofMatrix: psiU1ComReal))
//let (psiU1ComRealInv, psiU1ComImagInv) = invertComplexMatrix(ofRealMatrix: psiU1, andofImagMatrix: psiU1ComImag)
//var psiU1ComRealInvVec = Array(psiU1ComRealInv.joined())
//var psiU1ComImagInvVec = Array(psiU1ComImagInv.joined())
//var splitComplexU1Inv = DSPDoubleSplitComplex(realp: &psiU1ComRealInvVec,
//                                         imagp: &psiU1ComImagInvVec)
//
//
//let nU1inv = vDSP_Length(psiU1ComRealInvVec.count)
//var complexU1Inv = [DSPDoubleComplex](repeating: DSPDoubleComplex(),
//                                 count: psiU1ComRealInvVec.count)
//
//
////DSPSplitComplex tempSplitComplex;
////tempSplitComplex.realp = new float[N/2];
////tempSplitComplex.imagp = new float[N/2];
//
//
////vDSP_ctozD(psiU2ComVec, 2,  , 1, nU1inv)
//////
//
//
//vDSP_ztocD(&splitComplexU1Inv, 1,
//           &complexU1Inv, 2,
//           nU1inv)
//var mresult =  [DSPDoubleComplex](repeating: DSPDoubleComplex(),
//               count: psiU1ComRealInvVec.count*psiU1ComRealInvVec.count)
//
////vDSP_zmmaD(&splitComplexU1Inv, 1, &splitComplexU1Inv, 1, &splitComplexU1Inv, 1, &mresult,
////vDSP_zmmaD(psiQ1Com, 1, complexU1Inv, 1, &mresult, 1, psiU1ComRealInvVec.count, psiU1ComRealInvVec.count, psiU1ComRealInvVec.count)
////print("mresult \(mresult)")    // returns [90.0, 140.0, 190.0, 280.0, 370.0, 270.0, 400.0, 530.0]
//
//let matKcomp = matrixTimesMatrix(ofMatrix:  psiQ1 , withMatrix: invertMatrix(ofMatrix:  psiU1)   , andScalarScalingFactor: 1.0)
//print()
//print()
//print("psiU1:")
//printMatrix(ofMatrix: psiU1)
//print()
//print()
//print()
//
//print()
//print("psiQ1:")
//printMatrix(ofMatrix: psiQ1)
//print()
//print()
//print()
//
//print()
//print("Stiffness Matrix:")
//printMatrix(ofMatrix: matK)
//print()
//print()
//print()
//
//
////print(vecPart2[1])
//
//var kMod = matK
//var fMod = f
//
//let bcskin = [0, 1, 18, 19]
//
//for cntBcs in 0..<bcskin.count {
//    for i in 0..<dim {
//        if i == bcskin[cntBcs] {
//            fMod[i] = 0
//        }
//    }
//    for i in 0..<dim {
//        for j in 0..<dim {
//            if i == bcskin[cntBcs] ||  j == bcskin[cntBcs] {
//                if i == j {
//                    kMod[i][j] = 1
//                } else {
//                    kMod[i][j] = 0
//                }
//            }
//        }
//    }
//}
//print()
//printMatrix(ofMatrix:kMod)
//print()
//print(fMod)
//
//let sol = linearSolve(ofMatrix: kMod, andVector: scalarTimesVector(ofScalar: 1.0, withVector: fMod) )
//
//
//
//
//let solShape: [[Double]] = stride(from: 0, to: sol.count, by: 2).map {
//    Array(sol[$0..<$0+2])}
//
//print(sol)
//print()
//printMatrix(ofMatrix: solShape)
//let ccc0 = matrixAdd(ofMatrix: [[1,2],[3,4]], withMatrix: [[2,2],[2,2]])
//print()
//print()
//printMatrix(ofMatrix: ccc0)
//
//
