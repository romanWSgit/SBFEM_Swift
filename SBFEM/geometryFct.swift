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
