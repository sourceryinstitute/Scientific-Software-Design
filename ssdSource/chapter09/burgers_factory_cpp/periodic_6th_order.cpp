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
#include <exception>
#include "periodic_6th_order.h"
#include "fmt.h"
#include "mat.h"
#include "gaussian_elimination.h"

struct periodic_6th_order_error : public std::exception {
    virtual ~periodic_6th_order_error() throw () {};
};


const real_t pi = 3.14159265f;  

const real_t first_coeff_6th[5] =
      {1.0/3.0, 0.0, 14.0/9.0, 1.0/9.0, 0.0};

const real_t second_coeff_6th[5] =
      {2.0/11.0, 0.0, 12.0/11.0, 3.0/11.0, 0.0};


typedef field::ptr_t ptr_t;

crd_t periodic_6th_order::x_node_ = crd_t();

void periodic_6th_order::set_grid(int grid_resolution)
{
    if (x_node_.empty()) {
        for (int i=0; i< grid_resolution;++i) {
            x_node_.push_back (2.0*pi*(i)/real_t(grid_resolution));
        }

    }
}

int periodic_6th_order::get_grid_size() {
    return (x_node_.size());

}

const crd_t & periodic_6th_order::get_grid () {
    return x_node_;
}


periodic_6th_order::periodic_6th_order(real_t (*initial) (real_t x),
                int grid_resolution) {
    set_grid(grid_resolution);

    f_ = get_grid();

    for (int i = 0; i < f_.size(); ++i) {
        f_.at(i) = initial(f_.at(i));
    }
}

periodic_6th_order::periodic_6th_order(const crd_t & other)
    : f_(other)
{}


periodic_6th_order::~periodic_6th_order() {
}


ptr_t periodic_6th_order::operator+=(ptr_t rhs) {

    Ref<periodic_6th_order> other = cast<periodic_6th_order> (rhs);

    if ((other == NULL) || (f_.size() != other->f_.size())) {
        std::cerr << "periodic_6th_order::operator+= " <<
            "rhs is invalid" << std::endl;

        throw periodic_6th_order_error();
    }

    for (int i = 0; i < f_.size(); ++i) {
        f_.at(i) += other->f_.at(i);
    }

    return ptr_t(this);
}


ptr_t periodic_6th_order::operator-(ptr_t rhs)  const {
    Ref<periodic_6th_order> other = cast<periodic_6th_order> (rhs);

    if ((other == NULL) || (f_.size() != other->f_.size())) {
        std::cerr << "periodic_6th_order::operator- " <<
            "rhs is invalid" << std::endl;

        throw periodic_6th_order_error();
    }

    Ref<periodic_6th_order> result =
        Ref<periodic_6th_order>(new periodic_6th_order(*this));

    for (int i = 0; i < f_.size(); ++i) {
        result->f_.at(i) -= other->f_.at(i);
    }

    return result;
}


ptr_t periodic_6th_order::operator*(real_t rhs) const {
    Ref<periodic_6th_order> result =
        Ref<periodic_6th_order>(new periodic_6th_order(*this));

    for (int i = 0; i < f_.size(); ++i) {
        result->f_.at(i) *= rhs;
    }

    return result;
}


ptr_t periodic_6th_order::operator*(ptr_t rhs)  const {
    Ref<periodic_6th_order> other = cast<periodic_6th_order> (rhs);

    if ((other == NULL) || (f_.size() != other->f_.size())) {
        std::cerr << "periodic_6th_order::operator* " <<
            "rhs is invalid" << std::endl;

        throw periodic_6th_order_error();
    }

    Ref<periodic_6th_order> result =
        Ref<periodic_6th_order>(new periodic_6th_order(*this));

    for (int i = 0; i < f_.size(); ++i) {
        result->f_.at(i) *= other->f_.at(i);
    }

    return result;
}


ptr_t periodic_6th_order::clone()  const {
    return ptr_t (new periodic_6th_order(*this));
}


void periodic_6th_order::output() const {
    for (int i = 0; i < f_.size(); ++i) {
        std::cout << fmt(x_node_.at(i)) << " " << fmt(f_.at(i))
            << std::endl;
    }
}

ptr_t periodic_6th_order:: x()  const {
    int nx = get_grid_size();
    real_t dx = 2.0 * pi /real_t(nx);

    // __Initialize coeffecient matrix A and right handside b__
    mat_t A;
    A.clear_resize(nx, nx, 0.0);

    crd_t b = crd_t(nx, 0.0);

    crd_t coeff;

    for (int i = 0; i < 5; ++i)
        coeff.push_back(first_coeff_6th[i]);

    int x_east, x_west, x_east_plus1,
        x_east_plus2,x_west_minus1,x_west_minus2;

    for (int i = 0; i < nx; ++i) {
        x_east = (i+1)%nx;
        x_west = nx -1 - (nx-i)%nx;

        if (i == 1) {
            x_east_plus1=x_east+1;
            x_east_plus2=x_east+2;
            x_west_minus1=nx-1;
            x_west_minus2=nx-2;
        }
        else if (i == 2) {
            x_east_plus1=x_east+1;
            x_east_plus2=x_east+2;
            x_west_minus1=0;
            x_west_minus2=nx-1;
        }
        else if (i == nx-2) {
            x_east_plus1=0;
            x_east_plus2=1;
            x_west_minus1=x_west-1;
            x_west_minus2=x_west-2;
        }
        else if (i == nx-3) {
            x_east_plus1=nx-1;
            x_east_plus2=0;
            x_west_minus1=x_west-1;
            x_west_minus2=x_west-2;
        }
        else {
            x_east_plus1=x_east+1;
            x_east_plus2=x_east+2;
            x_west_minus1=x_west-1;
            x_west_minus2=x_west-2;
        }

        A(i,x_west_minus1) =coeff.at(1);
        A(i,x_west)        =coeff.at(0);
        A(i,i)             =1.0;
        A(i,x_east)        =coeff.at(0);
        A(i,x_east_plus1)  =coeff.at(1);

        b.at(i) = (0.25*coeff.at(3)*
              (f_.at(x_east_plus1)-f_.at(x_west_minus1))+
              0.5*coeff.at(2)*(f_.at(x_east) - f_.at(x_west)) +
              coeff.at(4)/6.0 * (f_.at(x_east_plus2)-
              f_.at(x_west_minus2)))/dx;
    }

    return ptr_t(new periodic_6th_order(gaussian_elimination(A,b)));
}


ptr_t periodic_6th_order:: xx()  const {
    int nx = get_grid_size();
    real_t dx = 2.0 * pi /real_t(nx);

    // __Initialize coeffecient matrix A and right handside b__
    mat_t A;
    A.clear_resize(nx, nx, 0.0);

    crd_t b = crd_t(nx, 0.0);

    crd_t coeff;

    for (int i = 0; i < 5; ++i)
        coeff.push_back(second_coeff_6th[i]);

    int x_east, x_west, x_east_plus1,
        x_east_plus2,x_west_minus1,x_west_minus2;

    for (int i = 0; i < nx; ++i) {
        x_east = (i+1)%nx;
        x_west = nx -1 - (nx-i)%nx;

        if (i == 1) {
            x_east_plus1=x_east+1;
            x_east_plus2=x_east+2;
            x_west_minus1=nx-1;
            x_west_minus2=nx-2;
        }
        else if (i == 2) {
            x_east_plus1=x_east+1;
            x_east_plus2=x_east+2;
            x_west_minus1=0;
            x_west_minus2=nx-1;
        }
        else if (i == nx-2) {
            x_east_plus1=0;
            x_east_plus2=1;
            x_west_minus1=x_west-1;
            x_west_minus2=x_west-2;
        }
        else if (i == nx-3) {
            x_east_plus1=nx-1;
            x_east_plus2=0;
            x_west_minus1=x_west-1;
            x_west_minus2=x_west-2;
        }
        else {
            x_east_plus1=x_east+1;
            x_east_plus2=x_east+2;
            x_west_minus1=x_west-1;
            x_west_minus2=x_west-2;
        }

        A(i,x_west_minus1) =coeff.at(1);
        A(i,x_west)        =coeff.at(0);
        A(i,i)             =1.0;
        A(i,x_east)        =coeff.at(0);
        A(i,x_east_plus1)  =coeff.at(1);

        b.at(i) = (0.25*coeff.at(3)*
              (f_.at(x_east_plus1)-2.0*f_.at(i)+f_.at(x_west_minus1))+
              coeff.at(2)*(f_.at(x_east)-2.0*f_.at(i)+f_.at(x_west)) +
              coeff.at(4)/9.0 * (f_.at(x_east_plus2)-
                  f_.at(i)+f_.at(x_west_minus2)))/(dx*dx);
    }

    return ptr_t(new periodic_6th_order(gaussian_elimination(A,b)));
}

real_t periodic_6th_order:: runge_kutta_2nd_step(real_t nu,
        int grid_resolution) const
{
    real_t dx, CFL, k_max;

    dx=2.0*pi/grid_resolution;

    k_max=grid_resolution*0.5;

    CFL=2.0/(24.0*(1-cos(k_max*dx))/11.0/(1.0+4.0/11.0*cos(k_max*dx))+ 
          3.0*(1.0-cos(2.0*k_max*dx))/22.0/(1.0+4.0/11.0*cos(k_max*dx)));
    return (CFL*dx*dx/nu);

}
