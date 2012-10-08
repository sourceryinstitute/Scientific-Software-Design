/*
 * Copyright (c) 2011, Damian Rouson, Jim Xia, and Xiaofeng Xu.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the names of Damian Rouson, Jim Xia, and Xiaofeng Xu nor the
 *       names of any other contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL DAMIAN ROUSON, JIM XIA, and XIAOFENG XU BE LIABLE 
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include "gaussian_elimination.h"
#include <cmath>

crd_t gaussian_elimination (const mat_t& lhs, const crd_t& rhs)
{
    const real_t pivot_tolerance=0.01;
    const int n = lhs.rows();

    if ( n != lhs.cols() || n != rhs.size()) {
        std::cerr << "gaussian_elimination: ill-posed system"
            << std::endl;
        throw gaussian_elimination_error();
    }

    // copy parameters to preserve lhs and rhs
    crd_t b = rhs;
    mat_t A = lhs;

    //-----  Gaussian elimination -----
    // forward elimination
    for (int p = 0; p < n-1; ++p) {
        if (fabs(A(p,p)) < pivot_tolerance) {
            std::cerr << "gaussian_elimination: use pivoting"
                << std::endl;
            throw gaussian_elimination_error();
        }
        for (int row = p+1; row < n; ++row) {
            real_t factor = A(row,p) / A(p,p);

            for (int col = p; col < n; ++col) {
                 A(row,col) -= A(p,col)*factor;
            }

            b.at(row) -= b.at(p)*factor;
        }
    }

    // Back substitution
    crd_t x(n);

    x.at(n-1) = b.at(n-1) / A(n-1, n-1);

    for (int row = n-2; row >=0; --row) {
        real_t sum = 0.0;

        for (int col = row+1; col < n; ++ col) {
            sum += A(row,col) * x.at(col);
        }
        x.at(row) = (b.at(row) - sum) / A(row, row);
    }
    return x;
}

