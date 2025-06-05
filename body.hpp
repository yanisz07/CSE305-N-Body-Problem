#pragma once

struct Vec { double x, y; };

struct Body {
    double m;
    Vec    pos, vel, force;
    Body(double m_, Vec p_, Vec v_)
      : m(m_), pos(p_), vel(v_), force{0,0} {}
};
