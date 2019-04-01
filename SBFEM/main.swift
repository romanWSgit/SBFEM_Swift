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


// MARK: - PARAMS
/// Machine Double Precission Epsilon
let eps64:Double = Double.ulpOfOne

/// ...
// FIXME: What is 'sd' ???
let sd: Int = 2
/// polynominal Order
let poly_ord: Int = 1
/// number of elements
let ielem: Int = 10
/// number of nodes per element
let nodeElem: Int = 2
/// number of all nodes in the domain
let nodedim: Int = 10
/// number of elements in 'a' - direction (horizontal)
let a: Double = 4
/// number of elements in 'b' - direction (vertical)
let b: Double = 1
// lengths in 'a' - direction (horizontal)
let l: Double = 24
// width in 'b' - direction (horizontal)
let d: Double = 6

let offset: [Double] = [0,0]
let centre: [Double] = [0,0]
let shapefct = "standard shape functions"
// Material
let Em: Double = 10000
let nu: Double = 0.2
let rho: Double = 1
let Er: Double = Em
let rr: Double = 1
let rhor: Double = rho

// MARK: - SERVICE FUNCTIONS
/**
 Initializes a new bicycle with the provided parts and specifications.
 
 - Parameters:
 - style: The style of the bicycle
 - gearing: The gearing of the bicycle
 - handlebar: The handlebar of the bicycle
 - frameSize: The frame size of the bicycle, in centimeters
 
 - Returns: A beautiful, brand-new bicycle,
 custom-built just for you.
 */
func def_lagrange(x: Double, i: Int, xm: [Double]) -> Double {
    let n:Int = xm.count
    var y:Double = 1
    for index in 0...n-1 {
        if i != index {
            y *= (x - xm[index]) / (xm[i] - xm[index])
        }
    }
    return y
}

/**
 Initializes a new bicycle with the provided parts and specifications.
 
 - Parameters:
 - style: The style of the bicycle
 - gearing: The gearing of the bicycle
 - handlebar: The handlebar of the bicycle
 - frameSize: The frame size of the bicycle, in centimeters
 
 - Returns: A beautiful, brand-new bicycle,
 custom-built just for you.
 */
func lagrange_diff(x: Double, i: Int, xm: [Double]) -> Double {
    var k:Double = 0
    let n:Int = xm.count
    var y:Double = 0
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
    return y
}

/**
 Initializes a new bicycle with the provided parts and specifications.
 
 - Parameters:
 - style: The style of the bicycle
 - gearing: The gearing of the bicycle
 - handlebar: The handlebar of the bicycle
 - frameSize: The frame size of the bicycle, in centimeters
 
 - Returns: A beautiful, brand-new bicycle,
 custom-built just for you.
 */
func points_table(poly_ord: Int) -> [Double] {
    var pointsTable:[Double] = []
    var result:Double = 0
    for i in 1...(poly_ord + 1) {
        result = -1 + (Double(i) - 1) * 2 / Double(poly_ord)
        pointsTable += [result]
    }
    return pointsTable
}
// Test STUFF
//print(def_lagrange(x: 0.5, i: 1, xm: [-1,1]))
//print(lagrange_diff(x: 0.5, i: 1, xm: [-1,1]))
//print(points_table(poly_ord: 10))


// MARK: Converted Stuff from C to Swift for numerical integration
var Pi:Double = atan2(1, 1) * 4
let N:Int = 2
var lroots:[Double] = Array(repeating: 0, count: N)
var weight:[Double] = Array(repeating: 0, count: N)
var lcoef:[[Double]] = Array(repeating: Array(repeating: 0, count: N + 1), count: N + 1)

func lege_coef()
{
    lcoef[0][0] = 1.0
    lcoef[1][1] = 1.0
    for n in 2...N {
        lcoef[n][0] = -1 * Double(n - 1) * (Double(lcoef[n - 2][0]) / Double(n))
        for i in 1...n {
            lcoef[n][i] = (Double((2 * n - 1)) * Double(lcoef[n - 1][i - 1]) - Double(n - 1) *  Double(lcoef[n - 2][i])) / Double(n)
        }
    }
}

func lege_eval(n: Int, x: Double) -> Double
{
    var s:Double = lcoef[n][n];
    for index in stride(from: n, to: 0, by: -1) {
        s = s * x + lcoef[n][index - 1]
    }
    return s
}

func lege_diff(n: Int, x: Double) -> Double
{
    let res: Double = Double(n) * (x * lege_eval(n: n, x: x) - lege_eval(n: n - 1, x: x)) / (x * x - 1)
    return res
}

func lege_roots()
{
    var x, x1: Double
    for index in 1...N
    {
        x = cos(Pi * (Double(index) - 0.25) / (Double(N) + 0.5))
        repeat
        {
            x1 = x
            x -= lege_eval(n: N, x: x) / lege_diff(n: N, x: x)
        } while ( fdim( x, x1) > 2e-16 )
            /*  fdim( ) was introduced in C99, if it isn't available
             *  on your system, try fabs( ) */
        lroots[index - 1] = x
            
        x1 = lege_diff(n: N, x: x)
        weight[index - 1] = 2 / ((1 - x * x) * x1 * x1)
    }
}
// TODO: How to convert pointers to Swift
//func lege_inte(double (*f)(double), double a, double b) -> Double
//{
//    double c1 = (b - a) / 2, c2 = (b + a) / 2, sum = 0;
//    int i;
//    for (i = 0; i < N; i++)
//    sum += weight[i] * f(c1 * lroots[i] + c2);
//    return c1 * sum;
//}

// MARK: - RUN
lege_coef()
lege_roots()

print("Roots: ");
for index in 0..<N {
    print(lroots[index], terminator:" ")
}
print("\n\nWeight:");
for index in 0..<N {
    print(weight[index], terminator:" ")
}
print("\n")

