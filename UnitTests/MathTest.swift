//
//  MathTest.swift
//  UnitTests
//
//  Created by Roman Wallner- Silberhuber on 14.01.23.
//  Copyright © 2023 Roman Wallner- Silberhuber. All rights reserved.
//

import XCTest



final class MathTest: XCTestCase {

    override func setUpWithError() throws {
        // Put setup code here. This method is called before the invocation of each test method in the class.
    }

    override func tearDownWithError() throws {
        // Put teardown code here. This method is called after the invocation of each test method in the class.
    }

    func testVectorAdd() {
            
        let start = DispatchTime.now() // <<<<<<<<<< Start time
            
        let res:[Double] = vectorAdd(ofVectorX: [1, 2, 3], withVectorY: [3, 4, -5 ])
            
        let end = DispatchTime.now()   // <<<<<<<<<<   end time
            


        let nanoTime = end.uptimeNanoseconds - start.uptimeNanoseconds // <<<<< Difference in nano seconds (UInt64)
        let timeInterval = Double(nanoTime) / 1_000_000_000 // Technically could overflow for long running tests

        print("Time to evaluate problem: \(timeInterval) seconds")
        // then
        XCTAssertEqual(res, [4.0, 6.0, -2.0], " matrix addition works")
    }
    
    func testScalarTimesVector() {
            
        let start = DispatchTime.now() // <<<<<<<<<< Start time
            
        let res:[Double] = scalarTimesVector(ofScalar: 0.5, withVector: [1.0, 2.0, 3.0, 4.9])
            
        let end = DispatchTime.now()   // <<<<<<<<<<   end time
            


        let nanoTime = end.uptimeNanoseconds - start.uptimeNanoseconds // <<<<< Difference in nano seconds (UInt64)
        let timeInterval = Double(nanoTime) / 1_000_000_000 // Technically could overflow for long running tests

        print("Time to evaluate problem: \(timeInterval) seconds")
        // then
        XCTAssertEqual(res, [0.5, 1.0, 1.5, 2.45], " scalarTimesVector works")
    }


//    func testPerformanceExample() throws {
//        // This is an example of a performance test case.
//        self.measure {
//            let _:[Double] = vectorAdd(ofVectorX: [1, 2, 3], withVectorY: [3, 4, -5 ])
//            
//        }
//    }

}
