#ifndef PSO_INPUT_VALIDATE_HPP
#define PSO_INPUT_VALIDATE_HPP

#include <src/input/parameters.hpp>

#include <string>
#include <vector>

namespace input{

/// @brief Checks every cross-field rule and returns all violations.
std::vector<std::string> validate(const Parameters& p);

} // namespace input

#endif // PSO_INPUT_VALIDATE_HPP
