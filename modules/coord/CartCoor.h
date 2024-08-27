// clkcomb - Clock and phase bias products Combination
// Copyright (C) 2021 Yuanxin Pan
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
// 3. Neither the name of the copyright holder nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef CARTCOOR_H
#define CARTCOOR_H

#include <math.h>

struct CartCoor
{
    CartCoor(double _x=0, double _y=0, double _z=0) : x(_x), y(_y), z(_z){}
    CartCoor(const double *p) : x(p[0]), y(p[1]), z(p[2]){}
    ~CartCoor(){}

    inline double dist(CartCoor coor)const {
        return sqrt((x - coor.x)*(x - coor.x) +
                    (y - coor.y)*(y - coor.y) +
                    (z - coor.z)*(z - coor.z));
    }

    inline double dot(CartCoor c)const {
        return x*c.x + y*c.y + z*c.z;
    }

    inline double norm()const {
        return sqrt(x*x + y*y + z*z);
    }

    inline double *data() {
        return &x;
    }

    inline const double *data()const {
        return &x;
    }

    inline double &operator[](int i) {
        return *(&x+i);
    }

    inline const double &operator[](int i)const {
        return *(&x+i);
    }

    inline CartCoor &operator+=(CartCoor c) {
        x += c.x;
        y += c.y;
        z += c.z;
        return *this;
    }

    inline CartCoor &operator-=(CartCoor c) {
        x -= c.x;
        y -= c.y;
        z -= c.z;
        return *this;
    }

    inline CartCoor &operator*=(double d) {
        x *= d;
        y *= d;
        z *= d;
        return *this;
    }

    inline CartCoor &operator/=(double d) {
        x /= d;
        y /= d;
        z /= d;
        return *this;
    }

    double x, y, z;
};

inline CartCoor operator+(CartCoor lhs, CartCoor rhs) {
    return CartCoor(lhs.x+rhs.x, lhs.y+rhs.y, lhs.z+rhs.z);
}

inline CartCoor operator-(CartCoor lhs, CartCoor rhs) {
    return CartCoor(lhs.x-rhs.x, lhs.y-rhs.y, lhs.z-rhs.z);
}

inline CartCoor operator*(const CartCoor c, double d) {
    return CartCoor(c.x*d, c.y*d, c.z*d);
}

inline CartCoor operator*(double d, const CartCoor c) {
    return CartCoor(c.x*d, c.y*d, c.z*d);
}

inline CartCoor operator/(const CartCoor c, double d) {
    return CartCoor(c.x/d, c.y/d, c.z/d);
}

#endif
