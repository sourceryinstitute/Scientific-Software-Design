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
#ifndef POINTER_H_
#define POINTER_H_

/**
 * A simple wrapper for raw pointer.
 */

#include <iostream>
#include <unistd.h>

template <typename T>
class pointer
{
public:
  pointer() : ptr_(NULL){}
  pointer(T *other) : ptr_(other){}

  void assign(T *other) { ptr_ = other; }

  pointer& operator=(T * other) {
    assign(other);
    return *this;
  }
  const pointer& operator=(const T & other) {
    assign(const_cast<T*>(&other));
    return *this;
  }
  pointer& operator=(T & other) {
    assign(&other);
    return *this;
  }

  T* operator->() { return ptr_; }
  const T* operator->() const { return ptr_; }

  T& operator*() { return *ptr_; }
  const T& operator*() const { return *ptr_; }

  bool operator==(const T *other) const { return ptr_ == other; }
  bool operator!=(const T *other) const { return ptr_ != other; }

  T* ptr() const { return ptr_; }

private:
    T *ptr_;
};
#endif
