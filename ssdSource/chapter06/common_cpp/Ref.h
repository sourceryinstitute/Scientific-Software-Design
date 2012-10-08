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
#ifndef REF_H_
#define REF_H_

/**
 * A very simple invasive reference counted pointer.
 */

#include "pointer.h"

template <typename T>
class Ref : virtual public pointer<T> {
public:
  typedef Ref<T> Self_t;

  Ref() {}

  Ref(const Self_t &other) {
    this->assign(other.ptr());
    if(pointer<T>::ptr() != NULL) pointer<T>::ptr()->grab();
  }

  template <typename Other>
  Ref(Ref<Other> other){
    this->assign(other.ptr());
    if(pointer<T>::ptr() != NULL) pointer<T>::ptr()->grab();
  }

  template <typename Other>
  Ref(Other *other) {
    this->assign(other);
    if(pointer<T>::ptr() != NULL) pointer<T>::ptr()->grab();
  }

  ~Ref() {
    if(pointer<T>::ptr() != NULL) pointer<T>::ptr()->release();
  }

  Self_t& operator=(const Self_t &other) {
    if(other.ptr() != pointer<T>::ptr()) {
      if(pointer<T>::ptr() != NULL) pointer<T>::ptr()->release();
      this->assign( other.ptr());
      if(pointer<T>::ptr() != NULL) pointer<T>::ptr()->grab();
    }
    return *this;
  }

  Self_t& operator=(T *other) {
    if(pointer<T>::ptr() != other) {
      if(pointer<T>::ptr() != NULL) pointer<T>::ptr()->release();
      this->assign( other);
      if(pointer<T>::ptr() != NULL) pointer<T>::ptr()->grab();
    }
    return *this;
  }
};

template <typename Target, typename Source>
Ref<Target> cast(const Ref<Source> &pp) {
  return Ref<Target>(dynamic_cast<Target*>(pp.ptr()));
}
#endif
