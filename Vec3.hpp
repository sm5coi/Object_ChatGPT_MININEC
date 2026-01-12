#ifndef VEC3_HPP
#define VEC3_HPP

#pragma once
#include <cmath>
#include <cassert>

struct Vec3
{
    double x;
    double y;
    double z;

    // --- Konstruktorer ---
    Vec3() : x(0.0), y(0.0), z(0.0) {}
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    // --- Operatorer ---
    Vec3 operator+(const Vec3& v) const
    {
        return Vec3(x + v.x, y + v.y, z + v.z);
    }

    Vec3 operator-(const Vec3& v) const
    {
        return Vec3(x - v.x, y - v.y, z - v.z);
    }

    Vec3 operator*(double s) const
    {
        return Vec3(x * s, y * s, z * s);
    }

    Vec3 operator/(double s) const
    {
        assert(s != 0.0);
        return Vec3(x / s, y / s, z / s);
    }

    Vec3& operator+=(const Vec3& v)
    {
        x += v.x; y += v.y; z += v.z;
        return *this;
    }

    Vec3& operator-=(const Vec3& v)
    {
        x -= v.x; y -= v.y; z -= v.z;
        return *this;
    }

    Vec3& operator*=(double s)
    {
        x *= s; y *= s; z *= s;
        return *this;
    }

    // --- Geometri ---
    double dot(const Vec3& v) const
    {
        return x * v.x + y * v.y + z * v.z;
    }

    double norm2() const
    {
        return dot(*this);
    }

    double norm() const
    {
        return std::sqrt(norm2());
    }

    Vec3 normalized() const
    {
        double n = norm();
        assert(n > 0.0);
        return *this / n;
    }
};

inline Vec3 operator*(double s, const Vec3& v)
{
    return v * s;
}



#endif // VEC3_HPP
