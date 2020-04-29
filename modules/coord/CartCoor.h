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
