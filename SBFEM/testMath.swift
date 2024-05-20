//
//  testMath.swift
//  SBFEM
//
//  Created by Roman Wallner- Silberhuber on 27.04.20.
//  Copyright © 2020 Roman Wallner- Silberhuber. All rights reserved.
//

import Foundation

// MARK: - TestCases

// Test vectorAdd
func mathTest() {
    
    let start = DispatchTime.now() // <<<<<<<<<< Start time
    
    let res:[Double] = vectorAdd(ofVectorX: [1, 2, 3], withVectorY: [3, 4, -5 ])
    
    let end = DispatchTime.now()   // <<<<<<<<<<   end time
    


    let nanoTime = end.uptimeNanoseconds - start.uptimeNanoseconds // <<<<< Difference in nano seconds (UInt64)
    let timeInterval = Double(nanoTime) / 1_000_000_000 // Technically could overflow for long running tests

    print("Time to evaluate problem: \(timeInterval) seconds")
    
   
    print("TEST")
    print(res)
}



