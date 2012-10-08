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
#include "timed_lorenz.h"
#include <exception>

struct timed_lorenz_error : public std::exception
{
  virtual ~timed_lorenz_error() throw() {}
};

typedef timed_lorenz::ptr_t ptr_t;
typedef timed_lorenz::strategy_t strategy_t;

timed_lorenz::timed_lorenz(const crd_t& ste, real_t s, real_t r, 
                           real_t b, strategy_t strat, double t_init) :
  lorenz(ste, s, r, b, strat), time_(t_init)
{}

timed_lorenz::~timed_lorenz()
{}

ptr_t timed_lorenz::clone() const
{
  return ptr_t(new timed_lorenz(*this));
}

ptr_t timed_lorenz::d_dt() const
{
  Ref<lorenz> parent = cast<lorenz>(lorenz::d_dt());
  return ptr_t(new timed_lorenz(parent->coordinate(), parent->sigma(),
               parent->rho(), parent->beta(),
               parent->get_strategy(), 1.0));
}

ptr_t timed_lorenz::operator+=(ptr_t inval)
{
  Ref<timed_lorenz> other = cast<timed_lorenz>(inval);
  if(other == NULL) {
    std::cerr << "timed_lorenz::operator+=:  Invalid input type\n";
    throw timed_lorenz_error();
  }
  lorenz::operator+=(other);
  time_ += other->time_;
  return ptr_t(this);
}

ptr_t timed_lorenz::operator*=(real_t val)
{
  lorenz::operator*=(val);
  time_ *= val;
  return ptr_t(this);
}

void timed_lorenz::set_time (real_t t)
{
  time_ = t;
}

real_t timed_lorenz::get_time() const
{
  return time_;
}
