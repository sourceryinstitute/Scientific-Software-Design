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
#include "ground.h"

const int ground::dim_ = 1;

ground::ground(real_t z, real_t beta) :
  z_(1, z), beta_(beta)
{}

const crd_t& ground::coordinate() const {
  return z_;
}

real_t ground::beta() const {
  return beta_;
}

ground ground::d_dt(const crd_t& x, const crd_t& y) const {
  return ground((x.at(0) * y.at(0) - beta_ * z_.at(0)), 0);
}

mat_t ground::d_dGround(const crd_t&, const crd_t&) const {
  mat_t result(dim_, dim_);
  result(0, 0) = -beta_;
  return result;
}

mat_t ground::d_dx(const crd_t &x, const crd_t &y) const {
  mat_t result(dim_, x.size());
  result(0, 0) = y.at(0);
  return result;
}

mat_t ground::d_dy(const crd_t &x, const crd_t &y) const {
  mat_t result(dim_, y.size());
  result(0, 0) = x.at(0);
  return result;
}

ground& ground::operator+=(const ground &other) {
  z_.at(0) += other.z_.at(0);
  beta_ += other.beta_;
  return *this;
}

ground& ground::operator-=(const ground& other) {
  z_.at(0) -= other.z_.at(0);
  beta_ -= other.beta_;
  return *this;
}

ground& ground::operator*=(real_t value) {
  z_.at(0) *= value;
  beta_ *= value;
  return *this;
}
