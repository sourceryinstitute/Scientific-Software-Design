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
#include <iostream>
#include <exception>

#include "lorenz.h"

using namespace std;

struct LorenzError : public std::exception {
  virtual ~LorenzError() throw() {}
};

// default constructor
lorenz::lorenz ()
{}

// constructor using each element
lorenz::lorenz (const crd_t initial_state, real_t s, real_t r, real_t b)
	:  state_(initial_state), sigma_(s), rho_(r), beta_(b)
{}

const crd_t& lorenz::output() const {
  return state_;
}

integrand::ptr_t lorenz:: d_dt() const
{
  ptr_t result = ptr_t(new lorenz);
  result->state_.resize(3);
  result->state_.at(0) = sigma_*(state_.at(1) - state_.at(0));
  result->state_.at(1) = state_.at(0)*(rho_ - state_.at(2))
                       - state_.at(1);
  result->state_.at(2) = state_.at(0)*state_.at(1) - beta_*state_.at(2);
  return result;
}


void lorenz::operator+=(integrand::ptr_t rhs) {
  ptr_t other = cast<lorenz>(rhs);
  if(other == NULL) {
    std::cerr << "lorenz::operator+=:  Failed dynamic cast\n";
    throw LorenzError();
  }
  if(other->state_.size() != this->state_.size()) {
    std::cerr << "lorenz::operator+=:  Non-identical dimensions.\n";
    throw LorenzError();
  }
  
  for(size_t i = 0; i < state_.size(); ++i) {
    state_.at(i) += other->state_.at(i);
  }
}

integrand::ptr_t lorenz::operator*(real_t rhs) const
{
  ptr_t result = ptr_t(new lorenz(*this));
  for(size_t i = 0; i < result->state_.size(); ++i) {
    result->state_.at(i) *= rhs;
  }
  return result;
}

lorenz::~lorenz()
{}
