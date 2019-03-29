//
// Created by bvdmitri on 10.02.19.
//

#include "ttsp_alternating.h"

namespace TTSP {

    AlternatingParameterHandler::AlternatingParameterHandler(const OptimizationParameter &parameter) :
            parameter(parameter), direction(AlternatingDirection::RIGHT), current_index(parameter.GetClosestIndexTo(parameter.GetDefaultValue())) {}

    AlternatingParameterHandler::AlternatingParameterHandler(const AlternatingParameterHandler &other) :
            parameter(other.parameter), direction(other.direction), current_index(other.current_index) {}

    const OptimizationParameter &AlternatingParameterHandler::GetParameter() const {
        return parameter;
    }

    AlternatingDirection AlternatingParameterHandler::GetDirection() const {
        return direction;
    }

    std::size_t AlternatingParameterHandler::GetCurrentIndex() const {
        return current_index;
    }

    std::size_t AlternatingParameterHandler::NextIndex() const {
        std::size_t count = parameter.GetValues().size() - 1;
        std::size_t index = current_index;

        if (count < 2) {
            return index;
        }

        switch_direction:
        switch (direction) {
            case AlternatingDirection::RIGHT:
                if (index == count) {
                    direction = AlternatingDirection::LEFT;
                    goto switch_direction;
                } else {
                    index += 1;
                }
                break;
            case AlternatingDirection::STAY:
                break;
            case AlternatingDirection::LEFT:
                if (index == 0) {
                    direction = AlternatingDirection::RIGHT;
                    goto switch_direction;
                } else {
                    index -= 1;
                }
                break;
            default:
                break;
        }
        return index;
    }

    void AlternatingParameterHandler::NextDirection() {
        switch (direction) {
            case AlternatingDirection::RIGHT:
                direction = AlternatingDirection::STAY;
                break;
            case AlternatingDirection::STAY:
                direction = AlternatingDirection::LEFT;
                break;
            case AlternatingDirection::LEFT:
                direction = AlternatingDirection::RIGHT;
                break;
            default:
                break;
        }
    }

    void AlternatingParameterHandler::UpdateIndex(std::size_t index) {
        current_index = index;
    }

    AlternatingOptimizer::AlternatingOptimizer(const std::string &name, const OptimizationParameters &space, const OptimizerProperties &properties, std::size_t buffer_capacity) :
            OptimizerInterface(name, space, properties, buffer_capacity), current_handler_index(0) {
        const OptimizationParameterEntries &parameters = space.GetParameterEntries();
        handlers.reserve(parameters.size());
        std::for_each(parameters.cbegin(), parameters.cend(), [this](const OptimizationParametersEntry &entry) {
            handlers.emplace_back(AlternatingParameterHandler(entry.first));
        });
    }

    SuggestionChangedParameters AlternatingOptimizer::AlgorithmMakeSuggestion(const std::function<OptimizationFunctionInvokeResult(const OptimizationParameterPoints &,
                                                                                                                                   const OptimizationParameterPoints &,
                                                                                                                                   void *)> &invoke, void *data) const {
        const AlternatingParameterHandler &handler = handlers.at(current_handler_index);
        return std::vector<SuggestionChangedParameter>{SuggestionChangedParameter(current_handler_index, handler.GetParameter().GetValues().at(handler.NextIndex()))};
    }

    bool AlternatingOptimizer::UpdateSpaceWithLatestResults() {
        AlternatingParameterHandler       &current = handlers.at(current_handler_index);
        const OptimizationParameterResult &last    = results.at(0);

        bool is_updated        = false;
        bool is_metrics_better = (last.GetMetricsBefore() < 0.0 || (last.GetMetricsAfter() < last.GetMetricsBefore()));
        bool is_direction_stay = (current.GetDirection() == AlternatingDirection::STAY);

        if (last.IsGood() && (is_metrics_better || is_direction_stay)) {
            current.UpdateIndex(current.NextIndex());
            parameters.Update(current_handler_index, parameters.GetParameter(current_handler_index).GetValues().at(current.GetCurrentIndex()), last.GetMetricsAfter());
            is_updated = true;
        }

        if (!is_metrics_better) {
            current.NextDirection();
        }

        current_handler_index = (current_handler_index + 1) % (handlers.size());

        return is_updated;
    }

    AlternatingOptimizer::~AlternatingOptimizer() {}


}
