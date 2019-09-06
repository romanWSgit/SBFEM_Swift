//
//  geometryFunctions.swift
//  SBFEM
//
//  Created by Roman Wallner- Silberhuber on 02.04.19.
//  Copyright © 2019 Roman Wallner- Silberhuber. All rights reserved.
//

import Foundation

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
func pointTable(generatedFromPolyOrd polyOrd: Int) -> [Double] {
    var pointsTableVar:[Double] = []
    var result:Double = 0
    for i in 1...(polyOrd + 1) {
        result = -1 + (Double(i) - 1) * 2 / Double(polyOrd)
        pointsTableVar += [result]
    }
    return pointsTableVar
}


func flatDiskGeom(forNumberOfElements nElem: Int, withPolyOrd pO: Int) -> (EZT: [[Int]], LTG: [[Int]]) {
    var Ezt1D: [Int] = []
    var res: [Int] = []
    for index in 1...nElem {
        res = [index, index + 1]
        Ezt1D += res
    }
    var EZT = Ezt1D.chunked(into: 2)
    EZT[nElem-1][1] = EZT[0][0]
    var Ltg1D: [Int] = []
    res = []
    for i in 1...nElem {
        for j in 0..<(pO + 1) * 2 {
            res = [pO * 2 * i - pO * 2  + 1 + j]
            Ltg1D += res
        }
    }
    var LTG = Ltg1D.chunked(into: 4)
    LTG[nElem-1][2] = EZT[0][0]
    LTG[nElem-1][3] = EZT[0][1]
    return (EZT, LTG)
}

func flatDiskMesh(forNumberOfElements nElem: Int, withPolyOrd pO: Int, withNodeDimension nodeDim: Int,
                  withHorizontalDivision a: Int, withVerticalDivision b: Int, withHeightD d: Double,
                  withLengthL l: Double, withScalingCentre scalingCentre: [Double] = [0.0, 0.0] ,
                  andOffset offS: [Double] = [0.0, 0.0] ) -> (nodes: [[Double]], nodesCentre: [[Double]],elementTable: [[Double]]) {
    var n: [Double] = []
    var res: [Double] = []
    var resC: [Double] = []
    var nC: [Double] = []
    for index in 1...nodeDim {
        let offX = offS[0]
        let offY = offS[1]
        if index <= (a + 1) {
            res = [offX + -1.0 * l / 2.0 + (l / Double(a)) * Double(index - 1), offY + -1.0 * d / 2.0]
          
        } else if (index - 1  >= a + 1) && (index - 1  <= a + b) {
            res = [offX + l / 2.0 , offY + (-d / 2 + Double(index - (a + 1)) * (d / Double(b)))]

        } else if (index > a + b) && (index <= nodeDim - b + 1) {
            res = [offX + (l / 2.0 - (l / Double(a)) * Double(index - (a + b) - 1)), offY + (d / 2.0)]

        } else if (index > nodeDim - b) && (index - 1 <= nodeDim) {
            res = [offX + (l / 2.0), offY + (d / Double(2) - Double(index - (nodeDim - b + 1)) * (d / Double(b)))]

        } else {
            print("wrong shit in Function: flatDiskMesh")
        }
        resC = vectorAdd(ofVectorX: res, withVectorY: scalarTimesVector(ofScalar: -1.0, withVector: scalingCentre))
        n += res
        nC += resC
    }
    let nodes = n.chunked(into: 2)
    let nodesC = nC.chunked(into: 2)
                                                                
    var el: [Double] = []
    for i in 1...nElem {
        for j in 1...pO + 1 {
            if !(i == nElem && j == pO + 1) {
                el += [nodes[j + (i - 1) * pO - 1][0], nodes[j - 1 + (i - 1) * pO][1]]
            } else {
                el += [nodes[0][0], nodes[0][1]]
            }
        }
    }
    let elements = el.chunked(into: 4)
    return (nodes, nodesC, elements)
}

func fCondListCantilever(forSbfemElement elem: Int, withPolyOrd pO: Int,
                         withHorizontalDivision a: Int, withVerticalDivision b: Int) -> Bool {
    return (a+b)/pO <= elem && elem < (2*a+b)/pO
}

