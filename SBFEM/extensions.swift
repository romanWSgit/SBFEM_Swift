//
//  extensions.swift
//  SBFEM
//
//  Created by Roman Wallner- Silberhuber on 09.04.19.
//  Copyright © 2019 Roman Wallner- Silberhuber. All rights reserved.
//

import Foundation

extension Array {
    func chunked(into size: Int) -> [[Element]] {
        return stride(from: 0, to: count, by: size).map {
            Array(self[$0 ..< Swift.min($0 + size, count)])
        }
    }
}
