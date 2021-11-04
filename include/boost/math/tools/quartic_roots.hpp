//  (C) Copyright Nick Thompson 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#ifndef BOOST_MATH_TOOLS_QUARTIC_ROOTS_HPP
#define BOOST_MATH_TOOLS_QUARTIC_ROOTS_HPP
#include <array>
#include <cmath>
#include <boost/math/tools/cubic_roots.hpp>

namespace boost::math::tools {

// Solves ax⁴ + bx³ + cx² + dx + e = 0.
// Only returns the real roots, as these are the only roots of interest in ray intersection problems.
// Follows Graphics Gems V: https://github.com/erich666/GraphicsGems/blob/master/gems/Roots3And4.c
template<typename Real>
std::array<Real, 4> quartic_roots(Real a, Real b, Real c, Real d, Real e) {
    using std::abs;
    using std::sqrt;
    std::array<Real, 4> roots(std::numeric_limits<Real>::quiet_NaN());
    if (std::abs(a) <= std::numeric_limits<Real>::min()) {
        auto cbrts = cubic_roots(b, c, d, e);
        roots[0] = cbrts[0];
        roots[1] = cbrts[1];
        roots[2] = cbrts[2];
        return roots;
    }
    if (std::abs(e) <= std::numeric_limits<Real>::min()) {
        auto v = cubic_roots(a, b, c, d);
        v.push_back(Real(0));
        return v;
    }
    // Now solve x⁴ + Ax³ + Bx² + Cx + D = 0.
    Real A = b/a;
    Real B = c/a;
    Real C = d/a;
    Real D = e/a;
    Real Asq = A*A;
    // Let x = y - A/4:
    // Mathematica: Expand[(y - A/4)^4 + A*(y - A/4)^3 + B*(y - A/4)^2 + C*(y - A/4) + D]
    // We now solve y⁴ + py² + qy + r = 0.
    Real p = B - 3*Asq/8;
    Real q = C - A*B/2 + Asq*A/8;
    Real r = D - A*C/4 + Asq*B/16 - 3*Asq*Asq/256;
    if (std::abs(r) <= std::numeric_limits<Real>::min()) {
        auto v = cubic_roots(Real(1), Real(0), p, q);
        v.push_back(Real(0));
        for (auto & y : v) {
            y -= A/4;
        }
        std::sort(v.begin(), v.end());
        return v;
    }
    // Biquadratic case:
    if (std::abs(q) <= std::numeric_limits<Real>::min()) {
        auto w = quadratic_roots(Real(1), p, r);
        std::vector<Real> v;
        for (auto r : w) {
            if (r >= 0) {
                Real rtr = sqrt(r);
                v.push_back(rtr - A/4);
                v.push_back(-rtr - A/4);
            }
        }
        std::sort(v.begin(), v.end());
        return v;
    }

    // Now split the depressed cubic into two quadratics:
    // y⁴ + py² + qy + r = (y² + sy + u)(y² - sy + v) = y⁴ + (v+u-s²)y² + s(v - u)y + uv
    // So p = v+u-s², q = s(v - u), r = uv.
    // Then (v+u)² - (v-u)² = 4uv = 4r = (p+s²)² - q²/s².
    // Multiply through by s² to get s²(p+s²)² - q² - 4rs² = 0, which is a cubic in s².
    // Then we let z = s², to get
    // z³ + 2pz² + (p² - 4r)z - q² = 0.
    auto z_roots = cubic_roots(Real(1), 2*p, p*p - 4*r, -q*q);
    // z = s², so s = sqrt(z).
    // No real roots:
    if (z_roots.back() <= 0) {
        std::vector<Real> v(0);
        return v;
    }
    Real s = std::sqrt(z_roots.back());

    // s is nonzero, because we took care of the biquadratic case.
    Real v = (p + s*s + q/s)/2;
    Real u = v - q/s;
    // Now solve y² + sy + u = 0:
    auto roots1 = quadratic_roots(Real(1), s, u);

    // Now solve y² - sy + v = 0:
    auto roots2 = quadratic_roots(Real(1), -s, v);
    for (auto root : roots2) {
        roots1.push_back(root);
    }

    for (auto& r : roots1) {
        r -= A/4;
    }

    // This is not super accurate. Clean up the roots with a Halley iterate.
    for (auto &r : roots1) {
        Real df = 4*a*r + 3*b;
        df = df*r + 2*c;
        df = df*r + d;
        Real d2f = 12*a*r + 6*b;
        d2f = d2f*r + 2*c;
        Real f = a*r + b;
        f = f*r + c;
        f = f*r + d;
        f = f*r + e;
        Real denom = 2*df*df - f*d2f;
        if (std::abs(denom) > std::numeric_limits<Real>::min())
        {
            r -= 2*f*df/denom;
        }
    }

    std::sort(roots1.begin(), roots1.end());
    return roots1;
}

}
#endif
