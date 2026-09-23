#ifndef PSO_INPUT_INITIALIZE_HPP
#define PSO_INPUT_INITIALIZE_HPP

#include <src/input/parameters.hpp>
#include <toml++/toml.hpp>
#include <initializer_list>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace input{

// Reads a TOML input file and hands out typed parameters.
//
// Every get_or_add() call does two jobs at once: it marks (block, key) as a
// parameter this build knows about, and it returns the value. Registration is
// therefore a side effect of asking, exactly as in Parthenon's ParameterInput.
// Whatever remains in the file that nobody asked for is a typo -- see
// unknown_keys(), which reports it with the line number toml++ recorded.
class InputFile{
public:

    /// @brief Parses a TOML file from disk.
    explicit InputFile(const std::string& path);

    /// @brief Wraps an already-parsed TOML table.
    explicit InputFile(toml::table table);

    /// @brief Reads an optional parameter, registering it as known.
    template<typename T>
    T get_or_add(std::string_view block, std::string_view key, T fallback){
        const auto node = mark(block, key);

        if (!node){
            return fallback;
        }

        const auto value = node.template value<T>();

        if (!value){
            throw std::runtime_error(location(node) + ": '" + name(block, key)
                                   + "' has the wrong type");
        }

        return *value;
    }

    /// @brief Reads a required parameter; throws when it is absent.
    template<typename T>
    T get(std::string_view block, std::string_view key){
        const auto node = mark(block, key);

        if (!node){
            throw std::runtime_error("missing required parameter '"
                                   + name(block, key) + "'");
        }

        const auto value = node.template value<T>();

        if (!value){
            throw std::runtime_error(location(node) + ": '" + name(block, key)
                                   + "' has the wrong type");
        }

        return *value;
    }

    /// @brief Reads a fixed-length numeric array, checking its length.
    template<size_t N>
    std::array<Real, N> get_array(std::string_view block, std::string_view key){
        const auto node = mark(block, key);

        if (!node){
            throw std::runtime_error("missing required parameter '"
                                   + name(block, key) + "'");
        }

        const auto* array = node.as_array();

        if (!array || array->size() != N){
            throw std::runtime_error(location(node) + ": '" + name(block, key)
                                   + "' must be an array of " + std::to_string(N)
                                   + " numbers");
        }

        std::array<Real, N> out{};

        for (size_t i = 0; i < N; i++){
            const auto value = array->get(i)->value<Real>();

            if (!value){
                throw std::runtime_error(location(node) + ": '" + name(block, key)
                                       + "[" + std::to_string(i) + "]' is not a number");
            }

            out[i] = *value;
        }

        return out;
    }

    /// @brief Reads an enumeration written as a string label.
    template<typename E>
    E get_or_add_enum(std::string_view block, std::string_view key, E fallback,
                      std::initializer_list<std::pair<std::string_view, E>> labels){
        const auto node = mark(block, key);

        if (!node){
            return fallback;
        }

        const auto text = node.template value<std::string>();

        if (!text){
            throw std::runtime_error(location(node) + ": '" + name(block, key)
                                   + "' must be a quoted string");
        }

        for (const auto& [label, value] : labels){
            if (label == *text){
                return value;
            }
        }

        std::string allowed;

        for (const auto& [label, value] : labels){
            allowed += (allowed.empty() ? "" : ", ");
            allowed += '"';
            allowed += label;
            allowed += '"';
        }

        throw std::runtime_error(location(node) + ": '" + name(block, key)
                               + "' is \"" + *text + "\"; expected one of " + allowed);
    }

    /// @brief Reads a high-precision scalar, keeping every digit of a quoted literal.
    LongReal get_or_add_long(std::string_view block, std::string_view key, LongReal fallback);

    /// @brief Lists file keys no accessor asked for, with their line numbers.
    std::vector<std::string> unknown_keys() const;

    /// @brief Returns every key this build understands.
    const std::set<std::string>& known_keys() const{
        return touched_;
    }

private:

    /// @brief Joins a block and key into a single dotted name.
    static std::string name(std::string_view block, std::string_view key){
        return std::string(block) + '.' + std::string(key);
    }

    /// @brief Renders a node's source line for an error message.
    static std::string location(const toml::node_view<toml::node>& node){
        return "line " + std::to_string(node.node()->source().begin.line);
    }

    /// @brief Records a key as known and returns whatever the file holds for it.
    toml::node_view<toml::node> mark(std::string_view block, std::string_view key){
        touched_.emplace(name(block, key));

        return tbl_[block][key];
    }

    toml::table tbl_;
    std::set<std::string> touched_;
};


/// @brief Binds an input file onto Parameters, field by field.
Parameters read_parameters(InputFile& in);

/// @brief Reads a file and rejects it if it names anything unknown.
Parameters read_parameters(const std::string& path);

} // namespace input

#endif // PSO_INPUT_INITIALIZE_HPP
