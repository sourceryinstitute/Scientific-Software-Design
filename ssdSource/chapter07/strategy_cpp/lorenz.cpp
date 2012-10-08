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
#include "lorenz.h"
#include "fmt.h"
#include <exception>

struct lorenz_error : public std::exception
{
  virtual ~lorenz_error() throw() {}
};

typedef lorenz::ptr_t ptr_t;

typedef lorenz::strategy_t strategy_t;

lorenz::lorenz(const crd_t& ste, real_t s, real_t r, real_t b, 
               strategy_t str) :
  integrand(str), state_(ste), sigma_(s), rho_(r), beta_(b)
{}
  
lorenz::~lorenz()
{}

ptr_t lorenz::clone() const
{
  return ptr_t(new lorenz(*this));
}

ptr_t lorenz::d_dt() const
{
  crd_t new_state(3);
  new_state.at(0) = sigma_ * (state_.at(1) - state_.at(0));
  new_state.at(1) = state_.at(0) * (rho_ - state_.at(2)) - state_.at(1);
  new_state.at(2) = state_.at(0) * state_.at(1) - beta_ * state_.at(2);
  return ptr_t(new lorenz(new_state, sigma_, rho_, beta_, 
               get_strategy()));
}

ptr_t lorenz::operator+=(ptr_t inval)
{
  Ref<lorenz> other = cast<lorenz>(inval);
  if((other == NULL) || (state_.size() != other->state_.size())) {
    std::cerr << "lorenz::operator+=:  Invalid input argument\n";
    throw lorenz_error();
  }
  size_t size = state_.size();
  for(size_t i = 0; i < size; ++i) {
    state_.at(i) += other->state_.at(i);
  }
  return ptr_t(this);
}

ptr_t lorenz::operator*=(real_t val)
{
  size_t size = state_.size();
  for(size_t i = 0; i < size; ++i) {
    state_.at(i) *= val;
  }
  return ptr_t(this);
}

void lorenz::set_coordinate(const crd_t& state)
{
  state_ = state;
}

const crd_t& lorenz::coordinate() const
{
  return state_;
}

real_t lorenz::sigma() const
{
  return sigma_;
}

real_t lorenz::rho() const
{
  return rho_;
}

real_t lorenz::beta() const
{
  return beta_;
}
