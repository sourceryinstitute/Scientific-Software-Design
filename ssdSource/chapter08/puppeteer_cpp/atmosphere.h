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
#ifndef _H_ATMOSPHERE_
#define _H_ATMOSPHERE_

#include "integrand.h"
#include "air.h"
#include "cloud.h"
#include "ground.h"

class atmosphere : public integrand {
public:
  typedef integrand::ptr_t ptr_t;
  typedef Ref<atmosphere> self_t;

  atmosphere(const air&, const cloud&, const ground&);
  virtual ~atmosphere();

  // The following methods do dynamic allocation.
  virtual ptr_t d_dt() const;
  virtual void  dRHS_dV(mat_t&) const;
  virtual ptr_t clone() const;
  virtual ptr_t  inverse_times(const mat_t&) const;
  virtual crd_t state_vector() const;
  
  // The following methods are destructive updates.
  virtual ptr_t operator+=(ptr_t);
  virtual ptr_t operator-=(ptr_t);
  virtual ptr_t operator*=(real_t);

private:
  air air_;
  cloud cloud_;
  ground ground_;
};

#endif //!_H_ATMOSPHERE_
