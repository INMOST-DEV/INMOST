//
// Created by bvdmitri on 10.02.19.
//

#include "ttsp_alternating.h"

namespace TTSP {

    AlternatingParameterHandler::AlternatingParameterHandler(const OptimizationParameter &parameter) :
            parameter(parameter), direction(AlternatingDirection::RIGHT), current_index(0) {}

    AlternatingParameterHandler::AlternatingParameterHandler(const AlternatingParameterHandler &other) :
            parameter(other.parameter), direction(other.direction), current_index(other.current_index) {}

    const OptimizationParameter &AlternatingParameterHandler::GetParameter() const {
        return parameter;
    }

    std::size_t AlternatingParameterHandler::GetDirection() const {
        return direction;
    }

    std::size_t AlternatingParameterHandler::GetCurrentIndex() const {
        return current_index;
    }

    std::size_t AlternatingParameterHandler::NextIndex() {
        std::size_t count = parameter.GetValues().size() - 1;
        std::size_t index = current_index;

        switch_direction:
        switch (direction) {
            case RIGHT:
                if (index == count) {
                    direction = LEFT;
                    goto switch_direction;
                } else {
                    index += 1;
                }
                break;
            case STAY1:
                break;
            case LEFT:
                if (index == 0) {
                    direction = RIGHT;
                    goto switch_direction;
                } else {
                    index -= 1;
                }
                break;
            case STAY2:
                break;
            default:
                break;
        }
        return index;
    }

    void AlternatingParameterHandler::NextDirection() {
        direction = (direction + 1) % 4;
    }

    void AlternatingParameterHandler::UpdateIndex(std::size_t index) {
        current_index = index;
    }

    AlternatingOptimizer::AlternatingOptimizer(const OptimizationParametersSpace &space, const OptimizerProperties &properties, std::size_t buffer_capacity) :
            OptimizerInterface(space, properties, buffer_capacity), current_handler_index(0) {
        const OptimizationParameters &parameters = space.GetParameters();
        handlers.reserve(parameters.size());
        std::for_each(parameters.begin(), parameters.end(), [this](const OptimizationParametersEntry &entry) {
            handlers.emplace_back(AlternatingParameterHandler(entry.first));
        });
    }

    OptimizationParametersSuggestion AlternatingOptimizer::Suggest(const std::function<OptimizationFunctionInvokeResult(const OptimizationParameterPoints &,
                                                                                                                        const OptimizationParameterPoints &,
                                                                                                                        void *)> &invoke, void *data) {

        OptimizationParameterPoints points(space.GetParameters().size());

        if (results.size() < 2) {
            if (results.size() == 0) {
                std::transform(handlers.begin(), handlers.end(), points.begin(), [](const AlternatingParameterHandler &h) {
                    return std::make_pair(h.GetParameter().GetName(), h.GetParameter().GetValues().at(h.GetCurrentIndex()));
                });
            } else if (results.size() == 1) {
                int i = 0;
                std::transform(handlers.begin(), handlers.end(), points.begin(), [this, &i](AlternatingParameterHandler &h) {
                    return std::make_pair(
                            h.GetParameter().GetName(),
                            h.GetParameter().GetValues().at(i++ == current_handler_index ? h.NextIndex() : h.GetCurrentIndex())
                    );
                });
            }
        } else {
            AlternatingParameterHandler       &current     = handlers.at(current_handler_index);
            const OptimizationParameterResult &last        = results.at(0);
            const OptimizationParameterResult &before_last = results.at(1);

            if (last.IsGood() && (last.GetMetricsAfter() < before_last.GetMetricsAfter())) {
                current.UpdateIndex(current.NextIndex());
            } else {
                current.NextDirection();
            }

            current_handler_index = (current_handler_index + 1) % (handlers.size());

            int i = 0;
            std::transform(handlers.begin(), handlers.end(), points.begin(), [this, &i](AlternatingParameterHandler &h) {
                return std::make_pair(
                        h.GetParameter().GetName(),
                        h.GetParameter().GetValues().at(i++ == current_handler_index ? h.NextIndex() : h.GetCurrentIndex())
                );
            });
        }

        return OptimizationParametersSuggestion(handlers.at(current_handler_index).GetParameter(), GetCurrentPoints(), points);
    }

    const OptimizationParameterPoints AlternatingOptimizer::GetCurrentPoints() const noexcept {
        OptimizationParameterPoints points(space.GetParameters().size());
        std::transform(handlers.cbegin(), handlers.cend(), points.begin(), [](const AlternatingParameterHandler &h) {
            return std::make_pair(h.GetParameter().GetName(), h.GetParameter().GetValues().at(h.GetCurrentIndex()));
        });
        return points;
    }

    AlternatingOptimizer::~AlternatingOptimizer() {}

}
