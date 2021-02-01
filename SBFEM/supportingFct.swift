//
//  soportingFunctions.swift
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
func printMatrix(ofMatrix matrix:[[Double]] ) {
    for array in matrix {
        print( array )
    }
}

func printMatrixDimensions (ofMatrix matrix: [[Double]] ) {
    print("Dimension: \(matrix.count) x \(matrix[0].count)")
}

func getIx(fromMatrix a: [[Double]], indexVector idx: [Int]) -> [[Double]] {
    let cnt = idx.count
    //print("idx:",idx)
    //print("count:",cnt)
    //var ret: [[Double]] = Array(repeating: Array(repeating: 0, count: cnt), count: cnt)
    var retFlat: [Double] = Array(repeating: 0, count: cnt * cnt)
    var run = 0
    var row = false
    let len = a.count
    var inci = 0
    var incj = 0
    for i in 0..<len {
        if idx[inci] == i {
            row = true
        }
        for j in 0..<len {
            print ("i:", i)
            print ("inci:", inci)
            print ("run:", run)
            if idx[incj] == j && row == true {
                //print ("j:", j)
                //print("VARIABLE SET")
                retFlat[run] = a[i][j]
                run += 1
                incj += 1
                if incj == cnt {
                    incj = 0
                    row = false
                    inci += 1
                }
            }
        }
    }
    let retNew: [[Double]] = stride(from: 0, to: retFlat.count, by: cnt).map{
        Array(retFlat[$0..<$0+cnt]) }
    return retNew
}

func addAtIx(inMatrix a: [[Double]], insertMatrix insert: [[Double]], atIndex idx: [Int]) -> [[Double]] {
    var aMod = a
    let insertFlat: [Double] = Array(insert.joined())
    var cntN = 0
    for i in idx {
        for j in idx {
            aMod[i][j] = aMod[i][j] + insertFlat[cntN]
            cntN += 1
        }
    }
    return aMod
}

func addAtIxVector(inVector a: [Double], insertVector insert: [Double], atIndex idx: [Int]) -> [Double] {
    var aMod = a
    var cntN = 0
    for i in idx {
        aMod[i] = aMod[i] + insert[cntN]
        cntN += 1
        }
    return aMod
}




func superMatrix(ofSubMat1 mat11: [[Double]], andSubMat2 mat12: [[Double]], andSubMat3 mat21: [[Double]], andSubMat4 mat22: [[Double]]) -> [[Double]] {
    if (mat11.count != mat11[0].count) || (mat12.count != mat12[0].count) ||
        (mat21.count != mat21[0].count) || (mat22.count != mat22[0].count) ||
        (mat12.count != mat12.count) || (mat12.count != mat21.count) ||
        (mat21.count != mat22.count) {
        print("Dimension Error in func: superMatrix(...)")
    }
    let cnt = mat11.count
    print(cnt)
    let superCnt = cnt * 2
    print(superCnt)
    var superMat: [[Double]] = Array(repeating: Array(repeating: 0, count: superCnt), count: superCnt)
    for i in 0..<superCnt {
        for j in 0..<superCnt {
            if i < cnt && j < cnt {
                superMat[i][j] = mat11[i][j]
            } else if (i >= cnt && i < superCnt) && (j < cnt) {
                superMat[i][j] = mat21[i-cnt][j]
            } else if (i < cnt) && (j >= cnt && j < superCnt) {
                superMat[i][j] = mat12[i][j-cnt]
            } else if (i >= cnt && i < superCnt) && (j >= cnt && j < superCnt) {
                superMat[i][j] = mat22[i-cnt][j-cnt]
            }
        }
    }
    return superMat
}
