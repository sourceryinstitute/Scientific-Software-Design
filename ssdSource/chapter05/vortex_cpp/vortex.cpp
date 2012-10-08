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
#include <cmath>
#include "vortex.h"
#include "fmt.h"

vortex::vortex()
{}

vortex::~vortex()
{}

void vortex::destroy()
{
    // the following line start destruction of the vortex ring
    // breaks the ring and triggers the chain action
    this->next_ = NULL;
}

vortex::vortex(real_t radius, int num_points)
{
    const real_t pi = 3.1415926f;
    real_t delta;
    iterator current;

    delta = 2.0f*pi*radius/real_t(num_points);

    this->point_.push_back(radius);
    this->point_.push_back(0.0f);
    this->point_.push_back(0.0f);

    this->next_ = new vortex();

    current = this->next_;

    for (int i = 0; i < num_points-1; ++i)
    {
        crd_t temp_point;
        real_t theta;

        theta = real_t(i+1) * delta/radius;

        temp_point.push_back(radius*cos(theta));
        temp_point.push_back(radius*sin(theta));
        temp_point.push_back(0.0);

        (*current)->point_ = temp_point;

        (*current)->next_ = new vortex();

        current = (*current)->next_;
    }

    *current = this; // this will cause the last memory allocation
                     // returned by the new operation in the last
                     // iteration to be automatically deleted
}

void vortex::output() const
{
    iterator current;
    crd_t tmp;

    tmp = this->point_;

    std::cout << fmt(tmp) << std::endl;

    current = this->next_;

    do
    {
        tmp = (*current)->point_;
        std::cout << fmt(tmp) << std::endl;

        current = (*current)->next_;
    } while (*current != this);
}
