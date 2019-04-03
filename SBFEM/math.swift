//
//  math.swift
//  SBFEM
//
//  Created by Roman Wallner- Silberhuber on 03.04.19.
//  Copyright © 2019 Roman Wallner- Silberhuber. All rights reserved.
//

import Foundation
import Accelerate
import simd
import Surge

func scalarTimesVector(ofScalar alpha: Double, withVector X: [Double]) -> [Double] {
    var x:[Double] = X
    let dx = x.count
    cblas_dscal(Int32(dx), alpha, &x, Int32(1))
    let Xnew: [Double] = x
    return Xnew
}

func dotProduct(ofVectorX vectorX: [Double], withVectorY vectorY: [Double]) -> Double {
    var x1 = vectorX
    var x2 = vectorY
    let dx1 = x1.count
    //var X:[Double] = Array(repeating: 0, count: dx1*dx2)
    
    let Xnew = cblas_ddot( Int32(dx1),
                           &x1,
                           1,
                           &x2,
                           1 )
    return Xnew
}

func dyadicP(withVectorX X: [Double], andVectorY Y: [Double]) -> [[Double]] {
    var x1:[Double] = X
    var x2:[Double] = Y
    let dx1 = x1.count
    let dx2 = x2.count
    var X:[Double] = Array(repeating: 0, count: dx1*dx2)
    
    cblas_dger(CblasRowMajor, /* you’re using row-major storage */
        Int32(dx1),           /* the matrix X has dx1 rows ...  */
        Int32(dx2),           /*  ... and dx2 columns.          */
        1.0,           /* scale factor to apply to x1x2' */
        &x1,
        1,             /* stride between elements of x1. */
        &x2,
        1,             /* stride between elements of x2. */
        &X,
        Int32(dx2))
    assert(X.count % dx2 == 0)
    let Xnew: [[Double]] = stride(from: 0, to: X.count, by: dx2).map {
        Array(X[$0..<$0+dx2])
    }
    return Xnew
}

func matrixVectorProduct(ofMatrixA A: [[Double]], withVectorX X: [Double]) -> [Double] {
    var a: [Double] = Array(A.joined())
    var x: [Double] = X
    let dx1 = A.count
    let dx2 = A[0].count
    var y  = [Double](repeating: 0, count: dx1)
    cblas_dgemv(CblasRowMajor, CblasNoTrans,  Int32(dx1), Int32(dx2), 1.0, &a, Int32(dx2), &x, Int32(1), 1.0, &y, Int32(1))
    let Ynew: [Double] = y
    return Ynew
}

func matrixVectorProduct(ofMatrixA A: [[Double]], withVectorX X: [Double], plusScaling scaling: Double, timesVector Y: [Double]) -> [Double] {
    var a: [Double] = Array(A.joined())
    var x: [Double] = X
    var y: [Double] = Y
    let dx1 = A.count
    let dx2 = A[0].count
    cblas_dgemv(CblasRowMajor, CblasNoTrans,  Int32(dx1), Int32(dx2), 1.0, &a, Int32(dx2), &x, Int32(1), scaling, &y, Int32(1))
    let Ynew: [Double] = y
    return Ynew
}
