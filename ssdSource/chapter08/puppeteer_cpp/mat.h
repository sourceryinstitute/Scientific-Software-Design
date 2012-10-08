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
#ifndef _H_MAT_
#define _H_MAT_

#include "globals.h"
#include <iostream>
#include <iomanip>

class mat_t {
public:
  typedef crd_t::value_type value_type;
  typedef crd_t::reference  reference;

  mat_t();
  mat_t(int rows, int cols);
  void clear();
  void resize(int rows, int cols);
  void clear_resize(int rows, int cols, value_type value = 0);
  void identity(int rows);
  int rows() const;
  int cols() const;
  value_type operator()(int r, int c) const;
  reference operator()(int r, int c);
  void set_submat(int r, int c, const mat_t &other);
  mat_t& operator-=(const mat_t&);
  mat_t& operator*=(real_t);

private:
  int r_, c_;
  crd_t data_;
};

inline mat_t operator-(const mat_t &a, const mat_t &b) {
  mat_t retval(a);
  retval -= b;
  return retval;
}

inline mat_t operator*(real_t value, const mat_t &matrix) {
  mat_t retval(matrix);
  retval *= value;
  return retval;
}

struct dim_t {
  const int eqs;
  const int vars;

  dim_t(int eqcnt, int varcnt) :
    eqs(eqcnt), vars(varcnt)
  {}
};

inline std::ostream& operator<<(std::ostream &os, const mat_t &mat) {
  std::ios_base::fmtflags flags = os.flags();
  for(int r = 0; r < mat.rows(); ++r) {
    os << "[";
    for(int c = 0; c < mat.cols(); ++c) {
      os << " " <<std::setw(12) <<std::setprecision(8)
      	 <<std::fixed <<mat(r,c);
    }
    os << "]\n";
  }
  os.flags(flags);
  return os;
}

#endif // !_H_MAT_
