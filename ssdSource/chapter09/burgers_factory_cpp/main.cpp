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
#include "field.h"
#include "field_factory.h"
#include "periodic_6th_factory.h"
#include "fmt.h"
#include "initializer.h"

int main ()
{
    typedef field::ptr_t field_ptr_t;
    typedef field_factory::ptr_t field_factory_ptr;

    const real_t t_final = 0.6;
    const real_t half = 0.5;
    const real_t nu = 1.0;
    const int grid_resolution = 128;

    real_t t = 0.0;

    field_factory_ptr field_creator =
        field_factory_ptr (new periodic_6th_factory());

    field_ptr_t u      = field_creator->create(u_initial,
                                grid_resolution);

    field_ptr_t u_half = field_creator->create(zero,
                                grid_resolution);

    real_t dt;

    while (t < t_final) {
        // use 2nd-order Runge-Kutta
        dt = u->runge_kutta_2nd_step(nu, grid_resolution);

        // first substep
        u_half = u + (u->xx()*nu - (u*u*half)->x())*dt*half;

        // second substep
        u += (u_half->xx()*nu - (u_half*u_half*half)->x())*dt;

        t += dt;
    }

    std::cout << " u at t = " << fmt(t) << std::endl;
    u->output();

    return 0;
}
