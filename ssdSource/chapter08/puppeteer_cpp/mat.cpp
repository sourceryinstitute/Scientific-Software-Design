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
#include "mat.h"
#include <exception>
#include <iostream>

struct matrix_error : public std::exception {
  virtual ~matrix_error() throw() {}
};

mat_t::mat_t() :
  r_(0), c_(0)
{}

mat_t::mat_t(int rows, int cols) {
  this->resize(rows, cols);
}

void mat_t::clear() {
  r_ = c_ = 0;
  data_.clear();
}

void mat_t::resize(int rows, int cols) {
  if(rows < 0 || cols < 0) {
    std::cerr << "mat_t::resize:  Rows and columns must be >= 0.\n";
    throw matrix_error();
  }
  if(! data_.empty()) {
    // Copy data over.
  }
  else {
    // common case.
    data_.resize(rows*cols);
    r_ = rows;
    c_ = cols;
  }
}

void mat_t::clear_resize(int rows, int cols, value_type value) {
  if(rows < 0 || cols < 0) {
    std::cerr << "mat_t::clear_resize:  Rows and columns must be >= 0\n";
    throw matrix_error();
  }
  data_.resize(rows*cols);
  r_ = rows;
  c_ = cols;
  std::fill(data_.begin(), data_.end(), value);
}

void mat_t::identity(int size) {
  this->clear_resize(size, size, 0);
  for(int i = 0; i < size; ++i) {
    this->operator()(i,i) = 1;
  }
}

int mat_t::rows() const {
  return r_;
}

int mat_t::cols() const {
  return c_;
}

mat_t::value_type mat_t::operator()(int r, int c) const {
  if(r < 0 || r >= r_ || c < 0 || c >= c_) {
    std::cerr << "mat_t::operator():  Invalid index (" << r << ", " << c
	      << ").  Bounds are (" << r_ << ", " << c << ")\n";
    throw matrix_error();
  }
  return data_.at(c*r_ + r);
}

mat_t::reference mat_t::operator()(int r, int c) {
  if(r < 0 || r >= r_ || c < 0 || c >= c_) {
    std::cerr << "mat_t::operator():  Invalid index (" << r << ", " << c
	      << ").  Bounds are (" << r_ << ", " << c << ")\n";
    throw matrix_error();
  }
  return data_.at(c*r_ + r);
}

void mat_t::set_submat(int startrow, int startcol, const mat_t &other) {
  for(int r = 0; r < other.rows(); ++r) {
    for(int c = 0; c < other.cols(); ++c) {
      this->operator()(r+startrow, c+startcol) = other(r,c);
    }
  }
}

mat_t& mat_t::operator-=(const mat_t &other) {
  if(this->rows() != other.rows() || this->cols() != other.cols()) {
    std::cerr <<
     "mat_t::operator-=:  Matrices must be of identical size.\n";
    throw matrix_error();
  }
  const size_t size = data_.size();
  for(size_t i = 0; i < size; ++i)
    data_[i] -= other.data_[i];
  return *this;
}

mat_t& mat_t::operator*=(real_t value) {
  for(crd_t::iterator it = data_.begin(); it != data_.end(); ++it)
    *it *= value;
  return *this;
}
