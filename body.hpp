// body.hpp
#pragma once

struct Vec {
    double x, y;
};

struct Body {
    double m;
    Vec    pos, vel, force;
    Body(double mass, Vec p, Vec v)
      : m(mass), pos(p), vel(v), force{0,0} {}
};
