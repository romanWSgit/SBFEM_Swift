//
//  numericalIntegration.swift
//  SBFEM
//
//  Created by Roman Wallner- Silberhuber on 02.04.19.
//  Copyright © 2019 Roman Wallner- Silberhuber. All rights reserved.
//

import Foundation
// MARK: - NUMERICAL INTEGRATION
func getRootsAndWeights(forIntegrationOrder N: Int) -> (Roots: [Double], Weights: [Double]) {
    let (roots, weights) = legeRoots(forIntegrationOrder: N, forCoefficentArray: calculateCoefficients(forIntegrationOrder: N))
    return (roots, weights)
}




func calculateCoefficients(forIntegrationOrder N: Int) -> [[Double]]
{
    var lcoef:[[Double]] = Array(repeating: Array(repeating: 0, count: N + 1), count: N + 1)
    lcoef[0][0] = 1.0
    lcoef[1][1] = 1.0
    for n in 2...N {
        lcoef[n][0] = -1 * Double(n - 1) * (Double(lcoef[n - 2][0]) / Double(n))
        for i in 1...n {
            lcoef[n][i] = (Double((2 * n - 1)) * Double(lcoef[n - 1][i - 1]) - Double(n - 1) *  Double(lcoef[n - 2][i])) / Double(n)
        }
    }
    return lcoef
}

func evalEdge(withSeperations n: Int, atPoint x: Double, forCoefficentArray lcoef: [[Double]]) -> Double
{
    var s:Double = lcoef[n][n];
    for index in stride(from: n, to: 0, by: -1) {
        s = s * x + lcoef[n][index - 1]
    }
    return s
}

func evalEdgeDifference(withSeperations n: Int, atPoint x: Double, forCoefficentArray lcoef: [[Double]]) -> Double
{
    let res: Double = Double(n) * (x * evalEdge(withSeperations: n, atPoint: x, forCoefficentArray: lcoef) - evalEdge(withSeperations: n - 1, atPoint: x, forCoefficentArray: lcoef)) / (x * x - 1)
    return res
}

func legeRoots(forIntegrationOrder N: Int, forCoefficentArray lcoef: [[Double]]) -> (weights: [Double], roots: [Double])
{
    var roots:[Double] = Array(repeating: 0, count: N)
    var weights:[Double] = Array(repeating: 0, count: N)
    var x, x1: Double
    for index in 1...N
    {
        x = cos(Double.pi * (Double(index) - 0.25) / (Double(N) + 0.5))
        repeat
        {
            x1 = x
            x -= evalEdge(withSeperations: N, atPoint: x, forCoefficentArray: lcoef) / evalEdgeDifference(withSeperations: N, atPoint: x, forCoefficentArray: lcoef)
        } while ( fdim( x, x1) > 2e-16 )
        /*  fdim( ) was introduced in C99, if it isn't available
         *  on your system, try fabs( ) */
        roots[index - 1] = x
        
        x1 = evalEdgeDifference(withSeperations: N, atPoint: x, forCoefficentArray: lcoef)
        weights[index - 1] = 2 / ((1 - x * x) * x1 * x1)
    }
    return (roots, weights)
}
