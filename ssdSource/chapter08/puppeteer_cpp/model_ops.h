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
#ifndef _H_INTEGRAND_OPS
#define _H_INTEGRAND_OPS

#include "integrand.h"
#include <exception>

struct model_ops_exception : public std::exception {
  virtual ~model_ops_exception() throw() {}
};

typedef Ref<integrand> ptr_t;
 
inline ptr_t operator+(ptr_t a, ptr_t b) {
  if(a == NULL || b == NULL) {
    std::cerr << "ptr_t + ptr_t:  Neither pointer can be NULL\n";
    throw model_ops_exception();
  }
  ptr_t c = a->clone();
  *c += b;
  return c;
}

inline ptr_t operator-(ptr_t a, ptr_t b) {
  if(a == NULL || b == NULL) {
    std::cerr << "ptr_t - ptr_t:  Neither pointer can be NULL\n";
    throw model_ops_exception();
  }
  ptr_t c = a->clone();
  *c -= b;
  return c;
}

inline ptr_t operator*(ptr_t a, real_t b) {
  if(a == NULL) {
    std::cerr << "ptr_t * real_t:  Pointer must not be NULL\n";
    throw model_ops_exception();
  }
  ptr_t c = a->clone();
  *c *= b;
  return c;
}

inline ptr_t operator*(real_t b, ptr_t a) {
  return a*b;
}

#endif //!_H_INTEGRAND_OPS
