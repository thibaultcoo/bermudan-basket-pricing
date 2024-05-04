#include "pch.h"
#include "SingleTrajectory.h"
#include "iostream"

SingleTrajectory::SingleTrajectory() {}

// Creates a single trajectory following Black-Scholes theory
SingleTrajectory::SingleTrajectory(double _start, double _end, size_t _nbSteps)
    : start(_start), end(_end), steps(_nbSteps)
{
    step = (end - start) / steps;
    Values.reserve(steps + 1);
    Times.reserve(steps + 1);
}

// Helper function to append a value to a trajectory
void SingleTrajectory::AddValue(double val) {
    if (Values.empty()) {
        Times.push_back(start);
    }
    else {
        Times.push_back(Times.back() + step);
    }
    Values.push_back(val);
}

// Hepler function to reach out to a specific value from a trajectory
const double SingleTrajectory::GetValue(double time) {
    auto lower = std::lower_bound(Times.begin(), Times.end(), time);
    if (lower == Times.end()) {
        return Values.back();
    }

    size_t i = std::distance(Times.begin(), lower);
    if (*lower == time) {
        return Values[i];
    }
    else if (i > 0) {
        // Perform linear interpolation
        double t1 = Times[i - 1], t2 = *lower;
        double v1 = Values[i - 1], v2 = Values[i];
        return v1 + (v2 - v1) * (time - t1) / (t2 - t1);
    }

    return Values.front();
}