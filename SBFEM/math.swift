//
//  math.swift
//  SBFEM
//
//  Created by Roman Wallner- Silberhuber on 03.04.19.
//  Copyright © 2019 Roman Wallner- Silberhuber. All rights reserved.
//

import Foundation
import Accelerate
import simd

/// Adds two vectors using the cblas_daxpy function.
///
/// This function computes a vector `Ynew` such that:
/// `Ynew = alpha * X + Y` where `alpha` is a scalar (in this case, 1.0).
///
/// - Parxameters:
///   - X: The first vector to be added.
///   - Y: The second vector which will be updated with the result.
/// - Returns: The resulting vector after the addition.
func vectorAdd(ofVectorX X: [Double], withVectorY Y: [Double]) -> [Double] {
    var x:[Double] = X
    var y:[Double] = Y
    let dx = x.count
    if #available(macOS 13.3, *) {
        cblas_daxpy(__LAPACK_int(Int32(dx)), 1.0, &x, __LAPACK_int(Int32(1)), &y, __LAPACK_int(Int32(1)))
    } else {
        // Log error to console or notify user
        print("Error: This application requires macOS 13.3 or later to function properly.")
        fatalError("Unsupported macOS version")
    }
    let Ynew: [Double] = y
    return Ynew
}

func scalarTimesVector(ofScalar alpha: Double, withVector X: [Double]) -> [Double] {
    var x:[Double] = X
    let dx = x.count
    if #available(macOS 13.3, *) {
        cblas_dscal(__LAPACK_int(Int32(dx)), alpha, &x, __LAPACK_int(Int32(1)))
    } else {
        // Log error to console or notify user
        print("Error: This application requires macOS 13.3 or later to function properly.")
        fatalError("Unsupported macOS version")
    }
    let Xnew: [Double] = x
    return Xnew
}

func scalarTimesMatrix(ofScalar beta: Double, withMatrix C: [[Double]]) -> [[Double]] {
    var c: [Double] = Array(C.joined())
    var a = c
    var b = c
    let dx1 = C.count
    let dx2 = C[0].count
    if #available(macOS 13.3, *) {
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, __LAPACK_int(Int32(dx1)), __LAPACK_int(Int32(dx2)), __LAPACK_int(Int32(dx1)), 0.0, &a,  __LAPACK_int(Int32(dx1)), &b, __LAPACK_int(Int32(dx2)), beta, &c, __LAPACK_int(Int32(dx2)))
    } else {
        // Log error to console or notify user
        print("Error: This application requires macOS 13.3 or later to function properly.")
        fatalError("Unsupported macOS version")
    }
    assert(c.count % dx2 == 0)
    let Cnew: [[Double]] = stride(from: 0, to: c.count, by: dx2).map{
        Array(c[$0..<$0+dx2])
    }
    return Cnew
}


func dotProduct(ofVectorX vectorX: [Double], withVectorY vectorY: [Double]) -> Double {
    var x1 = vectorX
    var x2 = vectorY
    let dx1 = x1.count
    //var X:[Double] = Array(repeating: 0, count: dx1*dx2)
    
    if #available(macOS 13.3, *) {
        let Xnew = cblas_ddot( __LAPACK_int(Int32(dx1)),
                               &x1,
                               1,
                               &x2,
                               1 )
        return Xnew
    } else {
        // Log error to console or notify user
        print("Error: This application requires macOS 13.3 or later to function properly.")
        fatalError("Unsupported macOS version")
    }
    
}

func dyadicP(withVectorX X: [Double], andVectorY Y: [Double]) -> [[Double]] {
    var x1:[Double] = X
    var x2:[Double] = Y
    let dx1 = x1.count
    let dx2 = x2.count
    var X:[Double] = Array(repeating: 0, count: dx1*dx2)
    
    if #available(macOS 13.3, *) {
        cblas_dger(CblasRowMajor, /* you’re using row-major storage */
                   __LAPACK_int(Int32(dx1)),           /* the matrix X has dx1 rows ...  */
                   __LAPACK_int(Int32(dx2)),           /*  ... and dx2 columns.          */
                   1.0,           /* scale factor to apply to x1x2' */
                   &x1,
                   1,             /* stride between elements of x1. */
                   &x2,
                   1,             /* stride between elements of x2. */
                   &X,
                   __LAPACK_int(Int32(dx2)))
        assert(X.count % dx2 == 0)
        let Xnew: [[Double]] = stride(from: 0, to: X.count, by: dx2).map {
            Array(X[$0..<$0+dx2])
        }
        return Xnew
    } else {
        // Log error to console or notify user
        print("Error: This application requires macOS 13.3 or later to function properly.")
        fatalError("Unsupported macOS version")
    }
    
}

func matrixVectorProduct(ofMatrixA A: [[Double]], withVectorX X: [Double]) -> [Double] {
    var a: [Double] = Array(A.joined())
    var x: [Double] = X
    let dx1 = A.count
    let dx2 = A[0].count
    var y  = [Double](repeating: 0, count: dx1)
    if #available(macOS 13.3, *) {
        cblas_dgemv(CblasRowMajor, CblasNoTrans,  __LAPACK_int(Int32(dx1)), __LAPACK_int(Int32(dx2)), 1.0, &a, __LAPACK_int(Int32(dx2)), &x, __LAPACK_int(Int32(1)), 1.0, &y, __LAPACK_int(Int32(1)))
        let Ynew: [Double] = y
        return Ynew
    } else {
        // Log error to console or notify user
        print("Error: This application requires macOS 13.3 or later to function properly.")
        fatalError("Unsupported macOS version")
    }
    
}

func matrixVectorProduct(ofMatrixA A: [[Double]], withVectorX X: [Double], plusScaling scaling: Double, timesVector Y: [Double]) -> [Double] {
    var a: [Double] = Array(A.joined())
    var x: [Double] = X
    var y: [Double] = Y
    let dx1 = A.count
    let dx2 = A[0].count
    cblas_dgemv(CblasRowMajor, CblasNoTrans,  __LAPACK_int(Int32(dx1)), __LAPACK_int(Int32(dx2)), 1.0, &a, __LAPACK_int(Int32(dx2)), &x, __LAPACK_int(Double(__LAPACK_int(Int32(1)))), scaling, &y, __LAPACK_int(Int32(1)))
    let Ynew: [Double] = y
    return Ynew
}

func matrixTimesMatrix(ofMatrix A: [[Double]], withMatrix B: [[Double]], andScalarScalingFactor alpha: Double) -> [[Double]] {
    var a: [Double] = Array(A.joined())
    var b: [Double] = Array(B.joined())
    var c: [Double] = Array(repeating: 0.0, count: A.count * B[0].count)
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, __LAPACK_int(Int32(A.count)), __LAPACK_int(Int32(B[0].count)), __LAPACK_int(Int32(A[0].count)), alpha, &a,  __LAPACK_int(Int32(A[0].count)), &b, __LAPACK_int(Int32(B[0].count)), 0.0, &c, __LAPACK_int(Int32(B[0].count)))
    let Cnew: [[Double]] = stride(from: 0, to: c.count, by: B[0].count).map{
        Array(c[$0..<$0+B[0].count])
    }
    return Cnew
}

func matrixTimesMatrixPlusC(ofMatrix A: [[Double]], withMatrix B: [[Double]], addMatrixC C: [[Double]], andScalarScalingFactor alpha: Double) -> [[Double]] {
    var a: [Double] = Array(A.joined())
    var b: [Double] = Array(B.joined())
    var c: [Double] = Array(C.joined())
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, __LAPACK_int(Int32(A.count)), __LAPACK_int(Int32(B[0].count)), __LAPACK_int(Int32(A[0].count)), alpha, &a,  __LAPACK_int(Int32(A[0].count)), &b, __LAPACK_int(Int32(B[0].count)),  1.0, &c, __LAPACK_int(Int32(B[0].count)))
    let Cnew: [[Double]] = stride(from: 0, to: c.count, by: B[0].count).map{
        Array(c[$0..<$0+B[0].count])
    }
    return Cnew
}

func matrixAdd(ofMatrix A: [[Double]], withMatrix B: [[Double]]) -> [[Double]] {
    let a: [Double] = Array(A.joined())
    let b: [Double] = Array(B.joined())
    var c: [Double] = Array(repeating: 0, count: A.count * B[0].count)
    for i in 0..<a.count {
            c[i] = a[i] + b[i]
    }
    
    let Cnew: [[Double]] = stride(from: 0, to: c.count, by: B[0].count).map{
        Array(c[$0..<$0+B[0].count])
    }
    return Cnew
}

func det2x2(of2x2MatrixA A: [[Double]]) -> Double {
    let a = simd_double2x2(rows: [
        simd_double2( A[0][0], A[0][1]),
        simd_double2( A[1][0], A[1][1])
        ])
    return simd_determinant(a)
}

func transposeMatrix(ofMatrix matrix:[[Double]] ) -> [[Double]] {
    
    var result = [[Double]](
        repeating: [Double]( repeating: 0, count: matrix.count ),
        count: matrix[ 0 ].count
    )
    
    for i in 0 ..< matrix.count {
        
        for k in 0 ..< matrix[ 0 ].count {
            
            result[ k ][ i ] = matrix[ i ][ k ]
        }
    }
    
    return result
}

func invertMatrix(ofMatrix matrix : [[Double]]) -> [[Double]] {
    let c: [Double] = Array(matrix.joined())
    var inMatrix: [Double]  = Array(matrix.joined())
    var N = __LAPACK_int(sqrt(Double(c.count)))
    var pivots = [__LAPACK_int](repeating: 0, count: Int(N))
    var workspace = [Double](repeating: 0.0, count: Int(N))
    var error : __LAPACK_int = 0
    
    withUnsafeMutablePointer(to: &N) {
        dgetrf_($0, $0, &inMatrix, $0, &pivots, &error)
        dgetri_($0, &inMatrix, $0, &pivots, &workspace, $0, &error)
    }
    let Cnew: [[Double]] = stride(from: 0, to: c.count, by: matrix[0].count).map{
        Array(inMatrix[$0..<$0+matrix[0].count])
    }
    return Cnew
}

func invertComplexMatrix(ofRealMatrix A : [[Double]], andofImagMatrix B: [[Double]]) -> (real: [[Double]], imag: [[Double]]) {
    let part1 = matrixAdd(ofMatrix: A, withMatrix: matrixTimesMatrix(ofMatrix: matrixTimesMatrix(ofMatrix: B, withMatrix: invertMatrix(ofMatrix: A), andScalarScalingFactor: 1.0), withMatrix: B, andScalarScalingFactor: 1.0) )
    let part2 = matrixAdd(ofMatrix: B, withMatrix: matrixTimesMatrix(ofMatrix: matrixTimesMatrix(ofMatrix: A, withMatrix: invertMatrix(ofMatrix: B), andScalarScalingFactor: 1.0), withMatrix: A, andScalarScalingFactor: 1.0) )
    let real = invertMatrix(ofMatrix: part1)
    let imag = scalarTimesMatrix(ofScalar: -1.0, withMatrix: invertMatrix(ofMatrix: part2))
    return (real, imag)
}

//func conditionNumber(ofMatrix matrix: [[Double]]) -> Double {
//
//    typealias LAInt = __CLPK_integer
//    var mat:[Double] = Array(matrix.joined())
//    let equations = 3
//    var numberOfEquations:LAInt = 3
//    var columnsInA:       LAInt = 3
//    var elementsInB:      LAInt = 3
//    var bSolutionCount:   LAInt = 1
//
//    dgecon_(<#T##__norm: UnsafeMutablePointer<Int8>!##UnsafeMutablePointer<Int8>!#>, <#T##__n: UnsafeMutablePointer<__CLPK_integer>!##UnsafeMutablePointer<__CLPK_integer>!#>, <#T##__a: UnsafeMutablePointer<__CLPK_complex>!##UnsafeMutablePointer<__CLPK_complex>!#>, <#T##__lda: UnsafeMutablePointer<__CLPK_integer>!##UnsafeMutablePointer<__CLPK_integer>!#>, <#T##__anorm: UnsafeMutablePointer<__CLPK_real>!##UnsafeMutablePointer<__CLPK_real>!#>, <#T##__rcond: UnsafeMutablePointer<__CLPK_real>!##UnsafeMutablePointer<__CLPK_real>!#>, <#T##__work: UnsafeMutablePointer<__CLPK_complex>!##UnsafeMutablePointer<__CLPK_complex>!#>, <#T##__rwork: UnsafeMutablePointer<__CLPK_real>!##UnsafeMutablePointer<__CLPK_real>!#>, <#T##__info: UnsafeMutablePointer<__CLPK_integer>!##UnsafeMutablePointer<__CLPK_integer>!#>)
//
//    return
//}

@available(macOS 13.3, *)
func eigensystem(ofMatrix matrix: [[Double]]) -> (eigenvaluesRe: [Double], eigenvaluesIm: [Double], rightEigenvectors: [[Double]]) {
    var mat: [Double] = Array(matrix.joined())
    
    var N = __LAPACK_int(sqrt(Double(mat.count)))
    var N1 = N
    var N2 = N
    var N3 = N
    //var workspaceQuery = [Double](repeating: 0.0, count: 1)
    var error : __LAPACK_int = 0
    var lwork = __LAPACK_int(-1)
    // Real parts of eigenvalues
    var wr = [Double](repeating: 0, count: Int(N)    )
    // Imaginary parts of eigenvalues
    var wi = [Double](repeating: 0, count: Int(N))
    // Left eigenvectors
    var vl = [Double](repeating: 0, count: Int(N*N))
    // Right eigenvectors
    var vr = [Double](repeating: 0, count: Int(N*N))
    var workspaceQuery: Double = 0.0

    let strrr = UnsafeMutablePointer<Int8>(mutating: ("V" as NSString).utf8String)

    //UnsafeMutablePointer(("V" as NSString).utf8CString)
    dgeev_(strrr!, strrr!, &N, &mat, &N1, &wr, &wi, &vl, &N2, &vr, &N3, &workspaceQuery, &lwork, &error)
    
    
    print("workspaceQuery: \(workspaceQuery)")
    // size workspace per the results of the query:
    var workspace = [Double](repeating: 0.0, count: Int(workspaceQuery))
//    var workspace = [Double](repeating: 0.0, count: 1400)
    lwork = __LAPACK_int(workspace.count)

    dgeev_(strrr!, strrr!, &N, &mat, &N1, &wr, &wi, &vl, &N2, &vr, &N3, &workspace, &lwork, &error)
    

   
    let eigWr: [Double] = vl
    let eigVRe: [Double] = wr
    let eigVIm: [Double] = wi
    let dx = eigVRe.count
    let eigWrNew: [[Double]] = stride(from: 0, to: dx*dx, by: dx).map {
            Array(eigWr[$0..<$0+dx])}
    return (eigVRe, eigVIm, eigWrNew)
}


// Compute the singular value decomposition of a matrix.
//
// Factors the matrix A as u * diag(s) * v, where u and v are unitary and s is a 1-d vector of a's singular values.
//
// - Parameter: matrix Matrix for which to calculate singular values.
// - Returns: Tuple containing a matrix U, a vector s, and a matrix V, such that A = U s V.H.
//
// - Notes:
// The decomposition is performed using LAPACK routine sgesvd_. The SVD is commonly written as a = U S V.H. The v
// returned by this function is V.H and u = U. If U is a unitary matrix, it means that it satisfies U.H = inv(U).
// The rows of v are the eigenvectors of a.H a. The columns of u are the eigenvectors of a a.H. For row i in v and
// column i in u, the corresponding eigenvalue is s[i]^2.

func singularValueDecomposition(ofMatrix A: [[Double]]) -> (U: [Double], s: [Double], V:[Double]){
    let AT  = transposeMatrix(ofMatrix: A)
    var ATvec = Array(AT.joined())
    var jobz1: Int8 = 65 // 'A'
    var jobz2 = jobz1
    
    var m = __LAPACK_int(A.count)
    var n = __LAPACK_int(A[0].count)
    
    var lda = m
    var ldu = m
    var ldvt = n
    
    // Allocate workspace size variables. By specifying an lWork of -1, we let LAPACK compute the optimal workspace.
    var wkOpt = Double(0.0)
    var lWork = __LAPACK_int(-1)
    var info = __LAPACK_int(0)
    //var iWork = [__CLPK_integer](repeating: 0, count: Int(8 * min(m, n)))
    
    // Create the output vectors and matrices
    var s: [Double] = Array(repeating: 0.0, count: Int(min(m, n)))
//    var U: [[Double]]  = Array(repeating: 0.0, rows: Int(ldu), columns: Int(m))
    let U = [[Double]](repeating: [Double]( repeating: 0.0, count: Int(ldu) ), count: Int(m))
//    var VT: [[Double]]  = Array(repeating: 0.0, rows: Int(ldvt), columns: Int(n))
    let VT = [[Double]](repeating: [Double]( repeating: 0.0, count: Int(ldvt) ), count: Int(n))
    // Query and allocate the optimal workspace
    var Uvec = Array(U.joined())
    var VTvec = Array(VT.joined())
    dgesvd_(&jobz1, &jobz2, &m, &n, &ATvec, &lda, &s, &Uvec, &ldu, &VTvec, &ldvt,  &wkOpt, &lWork, &info)
    //sgesdd_(&jobz, &m, &n, &_A.grid, &lda, &s, &U.grid, &ldu, &VT.grid, &ldvt, &wkOpt, &lWork, &iWork, &info)
   
    // Create the relevant workspace and update lWork to 0
    lWork = __LAPACK_int(wkOpt)
    var work = [Double](repeating: 0.0, count: Int(lWork))
    
    // Compute SVD. There are two possible lapack functions that work equally well?
    dgesvd_(&jobz1, &jobz2, &m, &n, &ATvec, &lda, &s, &Uvec, &ldu, &VTvec, &ldvt, &work, &lWork, &info)
    //sgesdd_(&jobz, &m, &n, &_A.grid, &lda, &s, &U.grid, &ldu, &VT.grid, &ldvt, &work, &lWork, &iWork, &info)
    
    // Check for convergence.
    
    
    return (Uvec, s, VTvec)
}

func linearSolve(ofMatrix matrix: [[Double]], andVector vector: [Double]) -> [Double] {
    typealias LAInt = __LAPACK_int
    var mat: [Double] = Array(matrix.joined())
    var b:[Double] = vector
    
    let equationsInt =  vector.count
    let equations =  LAInt(vector.count)
    var numberOfEquations:LAInt = equations
    var columnsInA:       LAInt = equations
    var elementsInB:      LAInt = equations
    var bSolutionCount:   LAInt = 1
    var outputOk: LAInt = 0
    var pivot = [LAInt](repeating: 0, count: equationsInt)
    
    dgesv_( &numberOfEquations, &bSolutionCount, &mat, &columnsInA, &pivot, &b, &elementsInB, &outputOk)
    
    return b
}
