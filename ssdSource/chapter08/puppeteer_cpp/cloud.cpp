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
#include "cloud.h"

const int cloud::dim_ = 1;

cloud::cloud(real_t y, real_t rho) :
  y_(1, y), rho_(rho)
{}

const crd_t& cloud::coordinate() const {
  return y_;
}

real_t cloud::rho() const {
  return rho_;
}

cloud cloud::d_dt(const crd_t &x, const crd_t &z) const {
  return cloud((x.at(0) * (rho_ - z.at(0)) - y_.at(0)), 0);
}

mat_t cloud::d_dCloud(const crd_t&, const crd_t&) const {
  mat_t result(dim_, dim_);
  result(0, 0) = -1;
  return result;
}

mat_t cloud::d_dx(const crd_t &x, const crd_t &z) const {
  mat_t result(dim_, x.size());
  result(0, 0) = rho_ - z.at(0);
  return result;
}

mat_t cloud::d_dz(const crd_t &x, const crd_t &z) const {
  mat_t result(dim_, z.size());
  result(0, 0) = -x.at(0);
  return result;
}

cloud& cloud::operator+=(const cloud &other) {
  y_.at(0) += other.y_.at(0);
  rho_ += other.rho_;
  return *this;
}
cloud& cloud::operator-=(const cloud &other) {
  y_.at(0) -= other.y_.at(0);
  rho_ -= other.rho_;
  return *this;
}

cloud& cloud::operator*=(real_t value) {
  y_.at(0) *= value;
  rho_ *= value;
  return *this;
}
