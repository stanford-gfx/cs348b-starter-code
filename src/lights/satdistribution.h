
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_LIGHTS_SATDISTRIBUTION_H
#define PBRT_LIGHTS_SATDISTRIBUTION_H

#include "sampling.h"

#include <functional>
#include <vector>

namespace pbrt {

class SummedAreaTable {
public:
    // |values| defines a piecewise-constant function over [0,1]^2
    SummedAreaTable(const Float *values, int nx, int ny)
        : nx(nx), ny(ny) {
        sum = integrate([&values,nx,ny](int x, int y) {
            return values[x + y * nx] / (nx * ny);
        });
    }

    Float Sum(const Bounds2f &extent) const {
        double s = ((lookup(extent.pMax.x, extent.pMax.y) -
                     lookup(extent.pMin.x, extent.pMax.y)) +
                    (lookup(extent.pMin.x, extent.pMin.y) -
                     lookup(extent.pMax.x, extent.pMin.y)));
        return std::max<Float>(s, 0);
    }

    Float Average(const Bounds2f &extent) const {
        return Sum(extent) / extent.Area();
    }

private:
    template <typename F> std::vector<double> integrate(F f) {
        std::vector<double> r(nx * ny);
        auto result = [&r,this](int x, int y) -> double & { return r[x + y * nx]; };

        result(0, 0) = f(0, 0);

        // sum across first scanline
        for (int x = 1; x < nx; ++x)
            result(x, 0) = f(x, 0) + result(x - 1, 0);

        // sum up first column
        for (int y = 1; y < ny; ++y)
            result(0, y) = f(0, y) + result(0, y - 1);

        // and all the rest of it
        for (int y = 1; y < ny; ++y)
            for (int x = 1; x < nx; ++x)
                result(x, y) = (f(x, y) + result(x - 1, y) +
                                result(x, y - 1) - result(x - 1, y - 1));

        return r;
    }

    double lookup(Float x, Float y) const {
        x *= nx + 1;
        y *= ny + 1;

        int x0 = (int)x, y0 = (int)y;

        Float v00 = lookup(x0, y0), v10 = lookup(x0 + 1, y0);
        Float v01 = lookup(x0, y0 + 1), v11 = lookup(x0 + 1, y0 + 1);

        // Bilinear interpolation
        Float dx = x - int(x), dy = y - int(y);
        return (1-dx)*(1-dy) * v00 + (1-dx)*dy * v01 + dx*(1-dy) * v10 + dx*dy * v11;
    }

    double lookup(int x, int y) const {
        // virtual zeros at lower boundaries
        if (x == 0 || y == 0)
            return 0;

        // reindex for actual stored values
        x = std::min(x - 1, nx - 1);
        y = std::min(y - 1, ny - 1);

        return sum[x + y * nx];
    }

    int nx, ny;
    std::vector<double> sum;
};

// Distribution2D based on a summed-area table.
class SATDistribution2D {
  public:
    SATDistribution2D(const Float *f, int nx, int ny)
        : sat(f, nx, ny), func(f, f + nx * ny), nx(nx), ny(ny) { }

    Point2f Sample(const Point2f &u, const Bounds2f &b, Float *pdf) const {
        if (sat.Sum(b) == 0) {
            *pdf = 0;
            return {};
        }

        // Marginal in first dimension
        auto Px = [&,this](Float x) -> Float {
            Bounds2f bx = b;
            bx.pMax.x = x;
            return sat.Sum(bx) / sat.Sum(b);
        };

        auto sample = [](std::function<Float(Float)> P,
                         Float u, Float min, Float max, int n) {
            while (std::ceil(n * max) - std::floor(n * min) > 1) {
                CHECK_LE(P(min), u);
                CHECK_GE(P(max), u);

                Float mid = (min + max) / 2;
                if (P(mid) > u)
                    max = mid;
                else
                    min = mid;
            }

            Float t = (u - P(min)) / (P(max) - P(min));
            //CHECK(t >= 0 && t <= 1) << t;
            return Clamp(Lerp(t, min, max), min, max);
        };

        Point2f p;
        p.x = sample(Px, u[0], b.pMin.x, b.pMax.x, nx);

        Bounds2f by(Point2f(std::floor(p.x * nx) / nx, b.pMin.y),
                    Point2f(std::ceil(p.x * nx) / nx, b.pMax.y));
        if (by.pMin.x == by.pMax.x)
            by.pMax.x += 1.f / nx;
        if (sat.Sum(by) <= 0) {
            // This can happen when we're provided a really narrow initial
            // bounding box, which happens in particular if the shading
            // point is in the plane of the portal.
            *pdf = 0;
            return {};
        }

        auto Py = [&, this](Float y) -> Float {
            Bounds2f byy = by;
            byy.pMax.y = y;
            return sat.Sum(byy) / sat.Sum(by);
        };
        p.y = sample(Py, u[1], b.pMin.y, b.pMax.y, ny);

        *pdf = Pdf(p, b);
        return p;
    }

    Float Pdf(const Point2f &p, const Bounds2f &b) const {
        if (sat.Sum(b) == 0)
            return 0;
        return Eval(p) / sat.Sum(b);
    }

  private:
    Float Eval(const Point2f &p) const {
        Point2i pi(std::min<int>(p[0] * nx, nx - 1),
                   std::min<int>(p[1] * ny, ny - 1));
        return func[pi[0] + pi[1] * nx];
    }

    SummedAreaTable sat;
    std::vector<Float> func;
    int nx, ny;
};

}

#endif
