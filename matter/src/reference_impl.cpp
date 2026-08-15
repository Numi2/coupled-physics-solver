#include "numi/matter/matter.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>

namespace numi::matter {
namespace {

using Matrix = std::array<double, 9>;
using Vector = std::array<double, 3>;

[[nodiscard]] bool finite(const double value) noexcept {
    return std::isfinite(value);
}

template <std::size_t Count>
[[nodiscard]] bool finite(const std::array<double, Count>& values) noexcept {
    return std::ranges::all_of(values, [](const double value) {
        return finite(value);
    });
}

void fail(std::string* error, const std::string& message) {
    if (error != nullptr) *error = message;
}

[[nodiscard]] Matrix transposed(const Matrix& value) noexcept {
    Matrix result{};
    for (std::size_t row = 0u; row < 3u; ++row) {
        for (std::size_t column = 0u; column < 3u; ++column) {
            result[3u * row + column] = value[3u * column + row];
        }
    }
    return result;
}

[[nodiscard]] Matrix multiplied(
    const Matrix& left,
    const Matrix& right
) noexcept {
    Matrix result{};
    for (std::size_t row = 0u; row < 3u; ++row) {
        for (std::size_t column = 0u; column < 3u; ++column) {
            for (std::size_t inner = 0u; inner < 3u; ++inner) {
                result[3u * row + column] +=
                    left[3u * row + inner] * right[3u * inner + column];
            }
        }
    }
    return result;
}

[[nodiscard]] Matrix added(
    const Matrix& left,
    const Matrix& right
) noexcept {
    Matrix result{};
    for (std::size_t index = 0u; index < result.size(); ++index) {
        result[index] = left[index] + right[index];
    }
    return result;
}

[[nodiscard]] Matrix scaled(
    const Matrix& value,
    const double scale
) noexcept {
    Matrix result{};
    for (std::size_t index = 0u; index < result.size(); ++index) {
        result[index] = scale * value[index];
    }
    return result;
}

[[nodiscard]] double determinant(const Matrix& value) noexcept {
    return value[0] * (value[4] * value[8] - value[5] * value[7]) -
        value[1] * (value[3] * value[8] - value[5] * value[6]) +
        value[2] * (value[3] * value[7] - value[4] * value[6]);
}

[[nodiscard]] Matrix cofactor(const Matrix& value) noexcept {
    return {
        value[4] * value[8] - value[5] * value[7],
        value[5] * value[6] - value[3] * value[8],
        value[3] * value[7] - value[4] * value[6],
        value[2] * value[7] - value[1] * value[8],
        value[0] * value[8] - value[2] * value[6],
        value[1] * value[6] - value[0] * value[7],
        value[1] * value[5] - value[2] * value[4],
        value[2] * value[3] - value[0] * value[5],
        value[0] * value[4] - value[1] * value[3],
    };
}

[[nodiscard]] Matrix cofactorDirectional(
    const Matrix& value,
    const Matrix& direction
) noexcept {
    const auto product = [&](const std::size_t a, const std::size_t b) {
        return direction[a] * value[b] + value[a] * direction[b];
    };
    return {
        product(4u, 8u) - product(5u, 7u),
        product(5u, 6u) - product(3u, 8u),
        product(3u, 7u) - product(4u, 6u),
        product(2u, 7u) - product(1u, 8u),
        product(0u, 8u) - product(2u, 6u),
        product(1u, 6u) - product(0u, 7u),
        product(1u, 5u) - product(2u, 4u),
        product(2u, 3u) - product(0u, 5u),
        product(0u, 4u) - product(1u, 3u),
    };
}

[[nodiscard]] double contraction(
    const Matrix& left,
    const Matrix& right
) noexcept {
    double result = 0.0;
    for (std::size_t index = 0u; index < left.size(); ++index) {
        result += left[index] * right[index];
    }
    return result;
}

[[nodiscard]] Vector multiplied(
    const Matrix& value,
    const Vector& vector
) noexcept {
    Vector result{};
    for (std::size_t row = 0u; row < 3u; ++row) {
        for (std::size_t column = 0u; column < 3u; ++column) {
            result[row] += value[3u * row + column] * vector[column];
        }
    }
    return result;
}

[[nodiscard]] double dot(const Vector& left, const Vector& right) noexcept {
    return left[0] * right[0] + left[1] * right[1] + left[2] * right[2];
}

[[nodiscard]] Matrix outer(const Vector& left, const Vector& right) noexcept {
    Matrix result{};
    for (std::size_t row = 0u; row < 3u; ++row) {
        for (std::size_t column = 0u; column < 3u; ++column) {
            result[3u * row + column] = left[row] * right[column];
        }
    }
    return result;
}

[[nodiscard]] bool evaluateBytecode(
    const ScalarBytecode& bytecode,
    const MaterialProgram& material,
    const ConstitutiveDifferentialInput& input,
    const std::array<double, 9>& direction,
    double& output
) noexcept {
    if (bytecode.maximumStack > NM_EXPRESSION_STACK_CAPACITY) return false;
    std::array<double, NM_EXPRESSION_STACK_CAPACITY> stack{};
    std::size_t stackSize = 0u;
    const auto push = [&](const double value) {
        if (stackSize >= stack.size() || !finite(value)) return false;
        stack[stackSize++] = value;
        return true;
    };
    const auto pop = [&](double& value) {
        if (stackSize == 0u) return false;
        value = stack[--stackSize];
        return true;
    };
    for (const NMExpressionInstructionGPU& instruction :
         bytecode.instructions) {
        double left = 0.0;
        double right = 0.0;
        switch (instruction.opcode) {
        case NM_EXPR_CONSTANT:
            if (!push(static_cast<double>(instruction.immediate.x))) return false;
            break;
        case NM_EXPR_PARAMETER:
            if (instruction.index >= material.parameters.size() ||
                !push(material.parameters[instruction.index].defaultValue)) {
                return false;
            }
            break;
        case NM_EXPR_F:
            if (instruction.index >= input.deformation.size() ||
                !push(input.deformation[instruction.index])) return false;
            break;
        case NM_EXPR_DF:
            if (instruction.index >= direction.size() ||
                !push(direction[instruction.index])) return false;
            break;
        case NM_EXPR_RATE:
            if (instruction.index >= input.deformationRate.size() ||
                !push(input.deformationRate[instruction.index])) return false;
            break;
        case NM_EXPR_DT:
            if (!push(input.timestep)) return false;
            break;
        case NM_EXPR_TEMPERATURE:
            if (!push(input.temperature)) return false;
            break;
        case NM_EXPR_STATE:
        case NM_EXPR_NEXT_STATE:
            return false;
        case NM_EXPR_ADD:
        case NM_EXPR_SUBTRACT:
        case NM_EXPR_MULTIPLY:
        case NM_EXPR_DIVIDE:
        case NM_EXPR_MIN:
        case NM_EXPR_MAX:
            if (!pop(right) || !pop(left)) return false;
            if (instruction.opcode == NM_EXPR_ADD && !push(left + right)) return false;
            if (instruction.opcode == NM_EXPR_SUBTRACT && !push(left - right)) return false;
            if (instruction.opcode == NM_EXPR_MULTIPLY && !push(left * right)) return false;
            if (instruction.opcode == NM_EXPR_DIVIDE) {
                if (!(std::abs(right) > 1.0e-12) || !push(left / right)) return false;
            }
            if (instruction.opcode == NM_EXPR_MIN && !push(std::min(left, right))) return false;
            if (instruction.opcode == NM_EXPR_MAX && !push(std::max(left, right))) return false;
            break;
        case NM_EXPR_NEGATE:
        case NM_EXPR_LOG:
        case NM_EXPR_EXP:
        case NM_EXPR_SQRT:
        case NM_EXPR_ABS:
        case NM_EXPR_POW_INTEGER:
            if (!pop(left)) return false;
            if (instruction.opcode == NM_EXPR_NEGATE && !push(-left)) return false;
            if (instruction.opcode == NM_EXPR_LOG) {
                if (!(left > 0.0) || !push(std::log(left))) return false;
            }
            if (instruction.opcode == NM_EXPR_EXP && !push(std::exp(left))) return false;
            if (instruction.opcode == NM_EXPR_SQRT) {
                if (left < 0.0 || !push(std::sqrt(left))) return false;
            }
            if (instruction.opcode == NM_EXPR_ABS && !push(std::abs(left))) return false;
            if (instruction.opcode == NM_EXPR_POW_INTEGER) {
                if (instruction.integer < 0 && !(std::abs(left) > 1.0e-12)) return false;
                if (!push(std::pow(left, instruction.integer))) return false;
            }
            break;
        case NM_EXPR_CLAMP: {
            double upper = 0.0;
            double lower = 0.0;
            if (!pop(upper) || !pop(lower) || !pop(left) || upper < lower ||
                !push(std::clamp(left, lower, upper))) return false;
            break;
        }
        default:
            return false;
        }
    }
    if (stackSize != 1u || !finite(stack[0])) return false;
    output = stack[0];
    return true;
}

[[nodiscard]] bool passiveDifferential(
    const ConstitutiveProgram& program,
    const ConstitutiveDifferentialInput& input,
    Matrix& stress,
    Matrix& tangent
) noexcept {
    const std::array<double, 9> zero{};
    for (std::size_t component = 0u; component < stress.size(); ++component) {
        if (!evaluateBytecode(
                program.stress[component], program.material, input, zero,
                stress[component]) ||
            !evaluateBytecode(
                program.tangentVector[component], program.material, input,
                input.deformationDirection, tangent[component])) {
            return false;
        }
        if (program.dissipation.has_value()) {
            double viscousStress = 0.0;
            double viscousTangent = 0.0;
            if (!evaluateBytecode(
                    program.viscousStress[component], program.material, input,
                    zero, viscousStress) ||
                !evaluateBytecode(
                    program.viscousTangentVector[component], program.material,
                    input, input.rateDirection, viscousTangent)) {
                return false;
            }
            stress[component] += viscousStress;
            tangent[component] += viscousTangent;
        }
    }
    return finite(stress) && finite(tangent);
}

void mixedDifferential(
    const MixedMaterialSource& mixed,
    const ConstitutiveDifferentialInput& input,
    const Matrix& passiveStress,
    const Matrix& passiveTangent,
    Matrix& stress,
    Matrix& tangent
) {
    const Matrix& deformation = input.deformation;
    const Matrix& direction = input.deformationDirection;
    const double j = determinant(deformation);
    const Matrix cof = cofactor(deformation);
    const Matrix dcof = cofactorDirectional(deformation, direction);
    const double dj = contraction(cof, direction);
    const Matrix kirchhoff = multiplied(passiveStress, transposed(deformation));
    const Matrix dkirchhoff = added(
        multiplied(passiveTangent, transposed(deformation)),
        multiplied(passiveStress, transposed(direction))
    );
    const double hydrostatic =
        (kirchhoff[0] + kirchhoff[4] + kirchhoff[8]) / 3.0;
    const double dhydrostatic =
        (dkirchhoff[0] + dkirchhoff[4] + dkirchhoff[8]) / 3.0;
    Matrix deviatoric = kirchhoff;
    Matrix ddeviatoric = dkirchhoff;
    deviatoric[0] -= hydrostatic;
    deviatoric[4] -= hydrostatic;
    deviatoric[8] -= hydrostatic;
    ddeviatoric[0] -= dhydrostatic;
    ddeviatoric[4] -= dhydrostatic;
    ddeviatoric[8] -= dhydrostatic;
    const Matrix inverseTranspose = scaled(cof, 1.0 / j);
    const Matrix dinverseTranspose = added(
        scaled(dcof, 1.0 / j),
        scaled(cof, -dj / (j * j))
    );
    const double pressure =
        input.pressure + mixed.biotCoefficient * input.porePressure;
    stress = added(
        multiplied(deviatoric, inverseTranspose),
        scaled(cof, -pressure)
    );
    tangent = added(
        added(
            multiplied(ddeviatoric, inverseTranspose),
            multiplied(deviatoric, dinverseTranspose)
        ),
        scaled(dcof, -pressure)
    );

    Vector fibre = mixed.fibreDirection;
    const double fibreNorm = std::sqrt(dot(fibre, fibre));
    for (double& value : fibre) value /= fibreNorm;
    const Vector deformedFibre = multiplied(deformation, fibre);
    const Vector directionFibre = multiplied(direction, fibre);
    const double fibreLength = std::sqrt(dot(deformedFibre, deformedFibre));
    if (fibreLength > 1.0e-12 && mixed.maximumActiveTension > 0.0) {
        Vector currentFibre{};
        for (std::size_t index = 0u; index < currentFibre.size(); ++index) {
            currentFibre[index] = deformedFibre[index] / fibreLength;
        }
        const double projectedDirection = dot(currentFibre, directionFibre);
        Vector currentFibreDirection{};
        for (std::size_t index = 0u;
             index < currentFibreDirection.size(); ++index) {
            currentFibreDirection[index] =
                (directionFibre[index] -
                 currentFibre[index] * projectedDirection) /
                fibreLength;
        }
        const double tension = mixed.maximumActiveTension *
            std::clamp(input.activation, 0.0, 1.0);
        const Matrix cauchy = scaled(
            outer(currentFibre, currentFibre), tension
        );
        const Matrix dcauchy = scaled(
            added(
                outer(currentFibreDirection, currentFibre),
                outer(currentFibre, currentFibreDirection)
            ),
            tension
        );
        stress = added(stress, multiplied(cauchy, cof));
        tangent = added(
            tangent,
            added(multiplied(dcauchy, cof), multiplied(cauchy, dcof))
        );
    }
}

} // namespace

bool evaluateStateFreeConstitutiveDifferentialFP64(
    const ConstitutiveProgram& program,
    const ConstitutiveDifferentialInput& input,
    ConstitutiveDifferentialFP64& output,
    std::string* error
) {
    if (!program.material.internalState.empty() ||
        program.gpu.stateCount != 0u ||
        program.gpu.projectionKind != NM_MATERIAL_PROJECTION_GENERIC ||
        program.material.hint == ConstitutiveHint::polyconvexICNN) {
        fail(error,
            "FP64 differential supports state-free, non-learned generic projection only");
        return false;
    }
    if (!finite(input.deformation) ||
        !finite(input.deformationDirection) ||
        !finite(input.deformationRate) ||
        !finite(input.rateDirection) ||
        !finite(input.timestep) || input.timestep < 0.0 ||
        !finite(input.temperature) || !finite(input.pressure) ||
        !finite(input.porePressure) || !finite(input.activation)) {
        fail(error, "FP64 differential input is non-finite or has a negative timestep");
        return false;
    }
    const double j = determinant(input.deformation);
    if (!(j > 1.0e-12) || !finite(j)) {
        fail(error, "FP64 differential deformation is singular or inverted");
        return false;
    }
    const Vector fibre = program.material.mixed.fibreDirection;
    const double fibreNorm = std::sqrt(dot(fibre, fibre));
    if (!(fibreNorm > 1.0e-12) || !finite(fibreNorm)) {
        fail(error, "FP64 differential mixed fibre is invalid");
        return false;
    }

    ConstitutiveDifferentialFP64 candidate;
    candidate.determinant = j;
    if (!passiveDifferential(
            program, input, candidate.passiveFirstPiola,
            candidate.passiveTangent)) {
        fail(error, "FP64 compiled constitutive bytecode evaluation failed");
        return false;
    }
    mixedDifferential(
        program.material.mixed,
        input,
        candidate.passiveFirstPiola,
        candidate.passiveTangent,
        candidate.coupledFirstPiola,
        candidate.coupledTangent
    );
    if (!finite(candidate.coupledFirstPiola) ||
        !finite(candidate.coupledTangent)) {
        fail(error, "FP64 mixed constitutive result is non-finite");
        return false;
    }
    output = candidate;
    if (error != nullptr) error->clear();
    return true;
}

} // namespace numi::matter
