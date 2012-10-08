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
#ifndef _FIELD_H_
#define _FIELD_H_ 1

#include "RefBase.h"
#include "globals.h"

class field : public virtual RefBase {
    public:
        typedef Ref<field> ptr_t;

        virtual ~field() {};

        virtual ptr_t operator+=(ptr_t)                   = 0;
        virtual ptr_t operator-(ptr_t)              const = 0;
        virtual ptr_t operator*(real_t)             const = 0;
        virtual ptr_t operator*(ptr_t)              const = 0;

        virtual ptr_t x()                           const = 0;
        virtual ptr_t xx()                          const = 0;
        virtual real_t runge_kutta_2nd_step(real_t,
                                            int)    const = 0;

        virtual ptr_t clone()                       const = 0;
        virtual void  output()                      const = 0;
    protected:
        field() {};
};


//
// define operators to be used for field ptr_t
//
inline field::ptr_t operator+= (field::ptr_t a, field::ptr_t b) {
    *a += b;
    return a;
}

inline field::ptr_t operator+ (field::ptr_t a, field::ptr_t b) {
    field::ptr_t tmp = a->clone();

    *tmp += b;
    return tmp;
}

inline field::ptr_t operator- (field::ptr_t a, field::ptr_t b) {
    field::ptr_t tmp = a->clone();

    tmp = *tmp - b;
    return tmp;
}


inline field::ptr_t operator*(field::ptr_t a, field::ptr_t b) {
    field::ptr_t tmp = a->clone();

    tmp = *tmp * b;
    return tmp;
}

inline field::ptr_t operator*(field::ptr_t a, real_t b) {
    field::ptr_t tmp = a->clone();

    tmp = *tmp * b;
    return tmp;
}

#endif
