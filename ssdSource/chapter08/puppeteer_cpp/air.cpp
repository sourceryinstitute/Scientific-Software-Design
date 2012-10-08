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
#include "air.h"

const int air::dim_ = 2;

air::air(real_t x, real_t sigma) {
  x_.push_back(x);
  x_.push_back(sigma);
}

const crd_t& air::coordinate() const {
  return x_;
}

air air::d_dt(const crd_t &y) const {
  return air((x_.at(1) * (y.at(0) - x_.at(0))), 0);
}

mat_t air::d_dAir(const crd_t& y) const {
  mat_t result(dim_, dim_);
  result(0, 0) = -x_.at(1);
  result(0, 1) = y.at(0) - x_.at(0);
  return result;
}

mat_t air::d_dy(const crd_t &y) const {
  mat_t result(dim_, y.size());
  result(0, 0) = x_.at(1);
  return result;
}

air& air::operator+=(const air &other) {
  x_.at(0) += other.x_.at(0);
  x_.at(1)   += other.x_.at(1);
  return *this;
}

air& air::operator-=(const air &other) {
  x_.at(0) -= other.x_.at(0);
  x_.at(1)   -= other.x_.at(1);
  return *this;
}

air& air::operator*=(real_t value) {
  x_.at(0) *= value;
  x_.at(1) *= value;
  return *this;
}
