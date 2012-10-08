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
#include "atmosphere.h"
#include <exception>
#include <cmath>

struct atmosphere_error : public std::exception {
  virtual ~atmosphere_error() throw() {}
};

typedef atmosphere::ptr_t ptr_t;

atmosphere::atmosphere(const air &a, const cloud &c, const ground &g) :
  air_(a), cloud_(c), ground_(g)
{}

atmosphere::~atmosphere() {
}

ptr_t atmosphere::d_dt() const {
  return
    ptr_t(new atmosphere(air_.d_dt(cloud_.coordinate()),
       cloud_.d_dt(air_.coordinate(),ground_.coordinate()),
       ground_.d_dt(air_.coordinate(),cloud_.coordinate())));
}

void atmosphere::dRHS_dV(mat_t &result) const {
  // Figure out the required dimensions.
  dim_t adim = air_.dimensions();
  dim_t cdim = cloud_.dimensions();
  dim_t gdim = ground_.dimensions();
  // Resize the result matrix.
  result.clear_resize(adim.eqs  + cdim.eqs  + gdim.eqs,
          adim.vars + cdim.vars + gdim.vars);
  if(result.rows() != result.cols()) {
    std::cerr << "atmosphere::dRHS_dV:  Ill-formed problem:  total of "
        << result.rows() << " equations and " << result.cols()
        << " variables\n";
    throw atmosphere_error();
  }
  // dAir/dAir
  result.set_submat(0, 0, air_.d_dAir(cloud_.coordinate()));
  // dAir/dCloud
  result.set_submat(0, adim.vars, air_.d_dy(cloud_.coordinate()));
  // dAir/dGround is all zero -- skipping that one.
  // dCloud/dAir
  result.set_submat(adim.eqs, 0,
        cloud_.d_dx(air_.coordinate(),ground_.coordinate()));
  // dCloud/dCloud
  result.set_submat(adim.eqs, adim.vars,
        cloud_.d_dCloud(air_.coordinate(),ground_.coordinate()));
  // dCloud/dGround
  result.set_submat(adim.eqs, adim.vars+cdim.vars,
        cloud_.d_dz(air_.coordinate(),ground_.coordinate()));
  // dGround/dAir
  result.set_submat(adim.eqs+cdim.eqs, 0,
        ground_.d_dx(air_.coordinate(),cloud_.coordinate()));
  // dGround/dCloud
  result.set_submat(adim.eqs+cdim.eqs, adim.vars,
        ground_.d_dy(air_.coordinate(),cloud_.coordinate()));
  // dGround/dGround
  result.set_submat(adim.eqs+cdim.eqs, adim.vars+cdim.vars,
        ground_.d_dGround(air_.coordinate(),cloud_.coordinate()));
}

ptr_t atmosphere::clone() const {
  return ptr_t(new atmosphere(*this));
}

crd_t atmosphere::state_vector() const {
  const crd_t &cc = cloud_.coordinate(), &gc = ground_.coordinate();
  crd_t state_space = air_.coordinate();
  state_space.insert(state_space.end(), cc.begin(), cc.end());
  state_space.insert(state_space.end(), gc.begin(), gc.end());
  return state_space;
}


ptr_t atmosphere::inverse_times(const mat_t &lhs) const {
  static const real_t pivot_tolerance = 1e-2;

  const int n = lhs.rows();
  crd_t b = this->state_vector();
  if((n != lhs.cols()) || (size_t(n) != b.size())) {
    std::cerr <<"integrand::inverse_times:  ill-posed matrix problem\n";
    throw atmosphere_error();
  }
  crd_t x(n);
  mat_t A(lhs);
  for(int p = 0; p < n-1; ++p) { // forward elimination
    if(fabs(A(p,p)) < pivot_tolerance) {
      std::cerr << "integrand::inverse_times:  "
    << "use an algorithm with pivoting\n";
      throw atmosphere_error();
    }
    for(int row = p+1; row < n; ++row) {
      real_t factor = A(row,p) / A(p,p);
      for(int col = p; col < n; ++col) {
  A(row,col) = A(row,col) - A(p,col)*factor;
      }
      b.at(row) = b.at(row)- b.at(p)*factor;
    }
  }
  x.at(n-1) = b.at(n-1) / A(n-1,n-1); // back substitution
  for(int row = n-1; row >= 0; --row) {
    real_t the_sum = 0;
    for(int col = row+1; col < n; ++col) {
      the_sum += A(row,col) * x.at(col);
    }
    x.at(row) = (b.at(row) - the_sum) / A(row,row);
  }
  return ptr_t(new atmosphere(air(x.at(0), x.at(1)),
            cloud(x.at(2), cloud_.rho()),
            ground(x.at(3), ground_.beta())));
}

ptr_t atmosphere::operator+=(ptr_t other) {
  self_t added = cast<atmosphere>(other);
  if(other == NULL) {
    std::cerr << "atmosphere::operator+=:  Invalid input type\n";
    throw atmosphere_error();
  }
  air_    += added->air_;
  cloud_  += added->cloud_;
  ground_ += added->ground_;
  return ptr_t(this);
}

ptr_t atmosphere::operator-=(ptr_t other) {
  self_t subbed = cast<atmosphere>(other);
  if(other == NULL) {
    std::cerr << "atmosphere::operator-=:  Invalid input type\n";
    throw atmosphere_error();
  }
  air_    -= subbed->air_;
  cloud_  -= subbed->cloud_;
  ground_ -= subbed->ground_;
  return ptr_t(this);
}

ptr_t atmosphere::operator*=(real_t value) {
  air_    *= value;
  cloud_  *= value;
  ground_ *= value;
  return ptr_t(this);
}
