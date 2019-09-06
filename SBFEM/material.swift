//
//  material.swift
//  SBFEM
//
//  Created by Roman Wallner- Silberhuber on 09.04.19.
//  Copyright © 2019 Roman Wallner- Silberhuber. All rights reserved.
//

import Foundation
import Accelerate
import simd

// MARK: - D- Matrix

let dMat = scalarTimesMatrix(ofScalar: Em / (1 - pow(nu, 2)), withMatrix: [[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2]])

