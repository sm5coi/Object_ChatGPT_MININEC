#ifndef EXCITATION_HPP
#define EXCITATION_HPP

struct Excitation {
    int segment;        // index på feed-segment (0-baserat)
    double voltage;     // amplitud, t.ex. 1.0
    double phase;       // radianer (0 = 0°)
};

#endif // EXCITATION_HPP
