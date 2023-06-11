#ifndef MAQ_H
#define MAQ_H

// Multi-Armed Qini (maq).
// https://github.com/grf-labs/maq
// Distributed under the MIT License.

#include "Data.hpp"
#include "SolverOptions.hpp"
#include "Solver.hpp"

namespace maq {

template <class T> // TODO: not necesarry to have x be const T&?
Solver<T> make_solver(const T& x, const SolverOptions& options) {return Solver<T>(x, options);}

} // namespace maq

#endif // MAQ_H
