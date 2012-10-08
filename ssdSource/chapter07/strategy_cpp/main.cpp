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
#include "explicit_euler.h"
#include "runge_kutta_2nd.h"
#include "fmt.h"
#include <iostream>

int main ()
{
  typedef lorenz::ptr_t ptr_t;

  static const int num_steps=2000;
  const real_t  sigma=10., rho=28., beta=8./3., dt=0.01;
  crd_t initial_condition(3, 1.0);

  Ref<lorenz> attractor = new lorenz(initial_condition, sigma, rho, beta,
                                     new explicit_euler);
  std::cout << " lorenz attractor:\n" 
       << fmt(attractor->coordinate(), 12, 9) << std::endl;
  for (int step = 0; step < 4*num_steps; ++step)
  {
    attractor->integrate(0.25*dt);
    std::cout << fmt(attractor->coordinate(), 12, 9) << std::endl;
  }

  Ref<timed_lorenz> timed_attractor
    = new timed_lorenz(initial_condition, sigma, rho, beta,
                                       new runge_kutta_2nd);
  std::cout << "\n timed_lorenz attractor:\n" 
            << fmt(timed_attractor->get_time(), 12, 9) << " "
            << fmt(timed_attractor->coordinate(), 12, 9) << std::endl;
  for (int i = 0; i < num_steps; ++i)
  {
    timed_attractor->integrate(dt);
    std::cout << fmt(timed_attractor->get_time(), 12, 9) << " "
        << fmt(timed_attractor->coordinate(), 12, 9) << std::endl;
  }
  
  return 0;
}
