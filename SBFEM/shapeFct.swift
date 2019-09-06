//
//  shapeFct.swift
//  SBFEM
//
//  Created by Roman Wallner- Silberhuber on 02.04.19.
//  Copyright © 2019 Roman Wallner- Silberhuber. All rights reserved.
//

import Foundation

func lagrangeInterpolate(atPoint x: Double, withPolyOrd i: Int, atDataPoints xm: [Double],
                         calcDeriv deriv: Bool) -> Double {
    let n:Int = xm.count
    var y:Double
    if deriv {
        var k:Double = 0
        y = 0
        for l in 0...n-1 {
            if l != i {
                k = 1 / (xm[i] - xm[l])
                for m in 0...n-1 {
                    if (m != i) && (m != l) {
                        k = k * ((x-xm[m])/(xm[i]-xm[m]))
                    }
                }
                y = y + k
            }
        }
        
    } else {
        y = 1
        for index in 0...n - 1 {
            if i != index {
                y *= (x - xm[index]) / (xm[i] - xm[index])
            }
        }
    }
    return y
}


func shapeFunctionsAndDerivatives(fromEta eta: Double, andPolyOrd polyOrd: Int,
                                  calcDeriv deriv: Bool)-> (shapeVec: [Double], shapeMat: [[Double]]) {
        var nVec: [Double] = []
        var resultMEven: [Double] = []
        var resultMUnEven: [Double] = []
        var nMat: [[Double]] = []
        var result: Double
        for index in 0...polyOrd {
            result = lagrangeInterpolate(atPoint: eta, withPolyOrd: index,
                                         atDataPoints:
                pointTable(generatedFromPolyOrd: polyOrd),
                                         calcDeriv: false)
            if deriv {
                result = lagrangeInterpolate(atPoint: eta, withPolyOrd: index,
                                             atDataPoints:
                    pointTable(generatedFromPolyOrd: polyOrd),
                                             calcDeriv: true)
            }
            resultMEven += [result]
            resultMEven += [0.0]
            resultMUnEven += [0.0]
            resultMUnEven += [result]
            nVec += [result]
        }
        nMat.append(resultMEven)
        nMat.append(resultMUnEven)
        return (nVec, nMat)
}


