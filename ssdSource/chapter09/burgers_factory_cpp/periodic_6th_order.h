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
#ifndef _PERIODIC_6TH_ORDER_H_
#define _PERIODIC_6TH_ORDER_H_ 1

#include <vector>
#include "field.h"
#include "globals.h"

class periodic_6th_order : public field {
    public:
        typedef field::ptr_t ptr_t;

        periodic_6th_order(real_t (*initial) (real_t x),
                        int grid_resolution);

        periodic_6th_order(const crd_t& other);
        virtual ~periodic_6th_order();

        virtual ptr_t operator+=(ptr_t);
        virtual ptr_t operator-(ptr_t)  const;
        virtual ptr_t operator*(real_t) const;
        virtual ptr_t operator*(ptr_t)  const;

        virtual ptr_t x()                               const;
        virtual ptr_t xx()                              const;
        virtual real_t runge_kutta_2nd_step(real_t nu,
                                  int grid_resolution)  const;

        virtual ptr_t clone()   const;
        virtual void output()   const;

        static void set_grid(int num_grid_pts);
        static const crd_t & get_grid();
        static int get_grid_size();

    private:
        crd_t f_;
        static crd_t x_node_;
};

#endif
