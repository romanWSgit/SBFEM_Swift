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
import Surge


// MARK: - PARAMS
/// Machine Double Precission Epsilon
let eps64 = Double.ulpOfOne

/// spatial Dimension
let sd = 2
/// polynominal Order
let polyOrd = 1
/// number of elements
let ielem = 10
/// number of nodes per element
let nodeElem = 2
/// number of all nodes in the domain
let nodedim = 10
/// number of elements in 'a' - direction (horizontal)
let a = 4
/// number of elements in 'b' - direction (vertical)
let b = 1
// lengths in 'a' - direction (horizontal)
let l = 24.0
// width in 'b' - direction (horizontal)
let d = 6.0

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





//func printMatrix( _ matrix:[[Double]] ) {
//    for array in matrix {
//        print( array )
//    }
//}

// MARK: - RUN
let (roots, weights) = getRootsAndWeights(forIntegrationOrder: 2)
print("Roots: ");
print(roots)
print("\n\nWeight:");
print(weights)
print("\n")
printMatrix(shapeFunctionsAndDerivatives(fromEta: 1, andPolyOrd: 2, calcDeriv: false).shapeMat)
printMatrix([rHatC(circumferentialCoordianteEta: 0.5, radialCoordinateXi: 1, coordianteOfElementsPoints: [-12.0, -3.0, -6.0, -3.0], withScalingCentre: [0.0, 0.0])] )


