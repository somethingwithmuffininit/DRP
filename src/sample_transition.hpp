#pragma once
#include <Eigen/SparseCore>

namespace mix {
using TransitionRows = Eigen::SparseMatrix<double, Eigen::RowMajor>;

// Use the last row entry for any rounding remainder of the cumulative sum.
inline int sample_transition(const TransitionRows &rows, int from, double uniform) {
    int state = from;
    for (TransitionRows::InnerIterator entry(rows, from); entry; ++entry) {
        state = entry.col();
        uniform -= entry.value();
        if (uniform <= 0)
            break;
    }
    return state;
}
} // namespace mix
