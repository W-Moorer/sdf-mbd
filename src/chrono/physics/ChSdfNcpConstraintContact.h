// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2024 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
//
// SDF-NCP contact constraints for Chrono's descriptor path.
//
// This module is intentionally independent from ChContactContainerNSC/SMC. It
// lets an SDF frontend provide gap/normal/contact-point samples, then injects
// frictionless unilateral normal constraints into the Chrono multibody solver.
//
// =============================================================================

#ifndef CH_SDF_NCP_CONSTRAINT_CONTACT_H
#define CH_SDF_NCP_CONSTRAINT_CONTACT_H

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "chrono/collision/ChSdfNcpContact.h"
#include "chrono/core/ChMatrix.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChPhysicsItem.h"
#include "chrono/solver/ChConstraintTwoBodies.h"

namespace chrono {
namespace sdfncp {

class ChSdfNcpPatchProjectedConstraint : public ChConstraintTwoBodies {
  public:
    ChSdfNcpPatchProjectedConstraint() = default;

    virtual ChSdfNcpPatchProjectedConstraint* Clone() const override {
        return new ChSdfNcpPatchProjectedConstraint(*this);
    }

    void ConfigurePatchProjection(bool enabled,
                                  int patch_id,
                                  const ChVector3d& point_abs,
                                  double area,
                                  double strength,
                                  double radius,
                                  const ChVector3d& force_per_pressure_abs = ChVector3d(0, 0, 0),
                                  const ChVector3d& torque_per_pressure_abs = ChVector3d(0, 0, 0)) {
        m_patch_projection_enabled = enabled;
        m_patch_id = patch_id;
        m_point_abs = point_abs;
        m_area = std::max(0.0, area);
        m_projection_strength = std::max(0.0, std::min(1.0, strength));
        m_projection_radius = std::max(0.0, radius);
        m_force_per_pressure_abs = force_per_pressure_abs;
        m_torque_per_pressure_abs = torque_per_pressure_abs;
        m_patch_neighbors.clear();
    }

    void SetPatchNeighbors(const std::vector<ChSdfNcpPatchProjectedConstraint*>& neighbors) {
        m_patch_neighbors = neighbors;
    }

    virtual void Project() override {
        ChConstraintTwoBodies::Project();
        if (!m_patch_projection_enabled || m_block_projection_enabled || m_projection_strength <= 0.0 ||
            m_patch_neighbors.empty()) {
            return;
        }

        const double projected = ComputeProjectedMultiplier(GetLagrangeMultiplier());
        SetLagrangeMultiplier(projected);
    }

    virtual bool ProjectGroupAndIncrementState(double old_lambda_current,
                                               double sharpness,
                                               bool record_violation_history,
                                               double& max_delta_lambda) override {
        if (!m_patch_projection_enabled || !m_block_projection_enabled || m_projection_strength <= 0.0 ||
            m_patch_neighbors.empty()) {
            return false;
        }

        std::vector<ChSdfNcpPatchProjectedConstraint*> group;
        group.reserve(m_patch_neighbors.size() + 1);
        group.push_back(this);
        for (auto* neighbor : m_patch_neighbors) {
            if (neighbor && neighbor != this && neighbor->IsActive() && neighbor->m_patch_id == m_patch_id) {
                group.push_back(neighbor);
            }
        }

        if (group.size() <= 1) {
            return false;
        }

        auto* leader = group.front();
        for (auto* constraint : group) {
            if (constraint < leader) {
                leader = constraint;
            }
        }

        // Only the first row encountered for a patch performs the block solve.
        // Other rows discard their scalar tentative update; the patch solution
        // already changed their multiplier and variable state consistently.
        if (this != leader) {
            SetLagrangeMultiplier(old_lambda_current);
            return true;
        }

        SolvePatchBlockLcpAndIncrementState(group, old_lambda_current, sharpness, record_violation_history,
                                            max_delta_lambda);
        return true;
    }

    void SetPatchBlockProjection(bool enabled,
                                 double laplacian_compliance = 0.0,
                                 double wrench_closure_strength = 0.0) {
        m_block_projection_enabled = enabled;
        m_block_laplacian_compliance = std::max(0.0, laplacian_compliance);
        m_block_wrench_closure_strength = std::max(0.0, wrench_closure_strength);
    }

  private:
    static void AddPatchGraphLaplacian(const std::vector<ChSdfNcpPatchProjectedConstraint*>& group,
                                       double compliance,
                                       std::vector<std::vector<double>>& W) {
        const size_t n = group.size();
        if (n <= 1 || compliance <= 0.0) {
            return;
        }

        double mean_area = 0.0;
        for (const auto* constraint : group) {
            mean_area += std::max(constraint ? constraint->m_area : 0.0, 1.0e-12);
        }
        mean_area /= static_cast<double>(n);
        double h = 2.0 * std::sqrt(std::max(mean_area, 1.0e-12));
        for (const auto* constraint : group) {
            if (constraint && constraint->m_projection_radius > 0.0) {
                h = constraint->m_projection_radius;
                break;
            }
        }
        h = std::max(h, 1.0e-9);
        const double radius2 = std::max(4.0 * h * h, 1.0e-24);

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                if (!group[i] || !group[j]) {
                    continue;
                }
                const double distance = (group[i]->m_point_abs - group[j]->m_point_abs).Length();
                if (!(distance > 1.0e-14) || distance > 3.0 * h) {
                    continue;
                }
                const double area_weight =
                    std::sqrt(std::max(group[i]->m_area * group[j]->m_area, 1.0e-24));
                const double spatial_weight = std::exp(-(distance * distance) / radius2);
                const double w = compliance * area_weight * spatial_weight;
                W[i][i] += w;
                W[j][j] += w;
                W[i][j] -= w;
                W[j][i] -= w;
            }
        }
    }

    static void AddPatchWrenchClosureRegularization(const std::vector<ChSdfNcpPatchProjectedConstraint*>& group,
                                                    double strength,
                                                    std::vector<std::vector<double>>& W) {
        const size_t n = group.size();
        if (n <= 1 || strength <= 0.0) {
            return;
        }

        double mean_diag = 0.0;
        double max_torque2 = 0.0;
        for (size_t i = 0; i < n; ++i) {
            mean_diag += std::abs(W[i][i]);
            if (group[i]) {
                max_torque2 = std::max(max_torque2, group[i]->m_torque_per_pressure_abs.Length2());
            }
        }
        mean_diag /= static_cast<double>(n);
        if (!(mean_diag > 1.0e-30) || !(max_torque2 > 1.0e-30)) {
            return;
        }

        const double scale = strength * mean_diag / max_torque2;
        for (size_t i = 0; i < n; ++i) {
            if (!group[i]) {
                continue;
            }
            for (size_t j = 0; j < n; ++j) {
                if (!group[j]) {
                    continue;
                }
                W[i][j] += scale * group[i]->m_torque_per_pressure_abs.Dot(group[j]->m_torque_per_pressure_abs);
            }
        }
    }

    void SolvePatchBlockLcpAndIncrementState(const std::vector<ChSdfNcpPatchProjectedConstraint*>& group,
                                             double old_lambda_current,
                                             double sharpness,
                                             bool record_violation_history,
                                             double& max_delta_lambda) {
        const size_t n = group.size();
        std::vector<double> old_lambda(n, 0.0);
        std::vector<double> base_cv(n, 0.0);
        std::vector<double> residual(n, 0.0);
        std::vector<std::vector<double>> W(n, std::vector<double>(n, 0.0));

        for (size_t i = 0; i < n; ++i) {
            old_lambda[i] = group[i] == this ? old_lambda_current : group[i]->GetLagrangeMultiplier();
            group[i]->SetLagrangeMultiplier(old_lambda[i]);
        }

        for (size_t i = 0; i < n; ++i) {
            base_cv[i] = group[i]->ComputeJacobianTimesState();
            residual[i] = base_cv[i] + group[i]->GetRightHandSide() +
                          group[i]->GetComplianceTerm() * old_lambda[i];
        }

        // Build the local Delassus block by perturbing descriptor variables with
        // each row's Eq_j = M^-1 J_j^T.  The perturbation is reverted immediately,
        // so the global solver state is unchanged until the final block update.
        for (size_t j = 0; j < n; ++j) {
            group[j]->IncrementState(1.0);
            for (size_t i = 0; i < n; ++i) {
                W[i][j] = group[i]->ComputeJacobianTimesState() - base_cv[i];
            }
            group[j]->IncrementState(-1.0);
        }
        for (size_t i = 0; i < n; ++i) {
            W[i][i] += group[i]->GetComplianceTerm();
        }
        AddPatchGraphLaplacian(group, m_block_laplacian_compliance, W);
        AddPatchWrenchClosureRegularization(group, m_block_wrench_closure_strength, W);

        std::vector<double> lambda = old_lambda;
        for (auto& value : lambda) {
            value = std::max(0.0, value);
        }

        const int local_iterations = std::max(12, static_cast<int>(4 * n));
        for (int iter = 0; iter < local_iterations; ++iter) {
            double max_step = 0.0;
            for (size_t i = 0; i < n; ++i) {
                double wi = residual[i];
                for (size_t j = 0; j < n; ++j) {
                    wi += W[i][j] * (lambda[j] - old_lambda[j]);
                }
                const double diag = std::max(std::abs(W[i][i]), 1.0e-30);
                const double updated = std::max(0.0, lambda[i] - wi / diag);
                max_step = std::max(max_step, std::abs(updated - lambda[i]));
                lambda[i] = updated;
            }
            if (max_step < 1.0e-10) {
                break;
            }
        }

        const double relaxation = std::max(0.0, std::min(1.0, m_projection_strength * sharpness));
        for (size_t i = 0; i < n; ++i) {
            const double new_lambda = std::max(0.0, old_lambda[i] + relaxation * (lambda[i] - old_lambda[i]));
            group[i]->SetLagrangeMultiplier(new_lambda);
            const double true_delta = new_lambda - old_lambda[i];
            group[i]->IncrementState(true_delta);
            if (record_violation_history) {
                max_delta_lambda = std::max(max_delta_lambda, std::abs(true_delta));
            }
        }
    }

    double ComputeProjectedMultiplier(double trial_lambda) const {
        trial_lambda = std::max(0.0, trial_lambda);
        if (!m_patch_projection_enabled || m_projection_strength <= 0.0 || m_patch_neighbors.empty()) {
            return trial_lambda;
        }

        double weighted_sum = 0.0;
        double weight_total = 0.0;
        const double radius2 = std::max(4.0 * m_projection_radius * m_projection_radius, 1.0e-24);
        for (const auto* neighbor : m_patch_neighbors) {
            if (!neighbor || neighbor == this || !neighbor->IsActive() || neighbor->m_patch_id != m_patch_id) {
                continue;
            }
            const double distance = (m_point_abs - neighbor->m_point_abs).Length();
            if (m_projection_radius > 0.0 && distance > 3.0 * m_projection_radius) {
                continue;
            }
            const double area_weight = std::sqrt(std::max(m_area * neighbor->m_area, 1.0e-24));
            const double spatial_weight = std::exp(-(distance * distance) / radius2);
            const double weight = area_weight * spatial_weight;
            weighted_sum += weight * std::max(0.0, neighbor->GetLagrangeMultiplier());
            weight_total += weight;
        }

        if (weight_total <= 1.0e-30) {
            return trial_lambda;
        }

        const double neighbor_average = weighted_sum / weight_total;
        return std::max(0.0, (1.0 - m_projection_strength) * trial_lambda + m_projection_strength * neighbor_average);
    }

    bool m_patch_projection_enabled = false;
    bool m_block_projection_enabled = false;
    int m_patch_id = -1;
    ChVector3d m_point_abs = ChVector3d(0, 0, 0);
    ChVector3d m_force_per_pressure_abs = ChVector3d(0, 0, 0);
    ChVector3d m_torque_per_pressure_abs = ChVector3d(0, 0, 0);
    double m_area = 0.0;
    double m_projection_strength = 0.0;
    double m_projection_radius = 0.0;
    double m_block_laplacian_compliance = 0.0;
    double m_block_wrench_closure_strength = 0.0;
    std::vector<ChSdfNcpPatchProjectedConstraint*> m_patch_neighbors;
};

/// One frictionless SDF-NCP contact sample between two Chrono rigid bodies.
///
/// The gap is defined as g = phi(x_A), where x_A is a world-space point attached
/// to body_a and phi is the SDF carried by body_b. The normal points in the
/// world direction of increasing gap. The resulting velocity Jacobian is
///
///   g_dot = n^T(v_A + w_A x r_A - v_B - w_B x r_B)
///
/// with body angular speeds represented in Chrono body-local coordinates.
///
/// The descriptor backend interprets lambda_n as a local normal contact
/// intensity. A positive force scale converts that intensity into the actual
/// sample force. When a physical sample area is available,
/// force_i = lambda_n_i * area_i * normal_i; otherwise weight provides a
/// backward-compatible unit-sample scale for older generators.
struct ChSdfNcpDescriptorContact {
    std::shared_ptr<ChBody> body_a;
    std::shared_ptr<ChBody> body_b;
    ChVector3d point_abs = ChVector3d(0, 0, 0);
    ChVector3d normal_abs = ChVector3d(0, 1, 0);
    double gap = 0.0;
    double weight = 1.0;
    double area = 0.0;
    double patch_area = 0.0;
    double gap_rate = 0.0;
    double predicted_gap = 0.0;
    bool use_custom_jacobian = false;
    ChVector3d custom_Cq_a_trans_abs = ChVector3d(0, 0, 0);
    ChVector3d custom_Cq_a_rot_local = ChVector3d(0, 0, 0);
    ChVector3d custom_Cq_b_trans_abs = ChVector3d(0, 0, 0);
    ChVector3d custom_Cq_b_rot_local = ChVector3d(0, 0, 0);
    bool use_custom_wrench = false;
    ChVector3d force_on_body_a_per_lambda_abs = ChVector3d(0, 0, 0);
    ChVector3d torque_on_body_a_per_lambda_abs = ChVector3d(0, 0, 0);
    ChVector3d force_on_body_b_per_lambda_abs = ChVector3d(0, 0, 0);
    ChVector3d torque_on_body_b_per_lambda_abs = ChVector3d(0, 0, 0);
    int contact_id = -1;
    int patch_id = -1;
};

struct ChSdfNcpDescriptorContactState {
    ChSdfNcpDescriptorContact sample;
    double lambda_n = 0.0;
    double lambda_force = 0.0;
    double penetration = 0.0;
    double ncp_residual = 0.0;
    double complementarity_error = 0.0;
    bool active = false;
};

/// Persistent SDF contact manifold manager.
///
/// This manager lives between an SDF geometry query frontend and an SDF-NCP
/// backend.  It is intentionally backend-neutral: it preserves sample identity,
/// maps transient connected components to persistent patch ids, applies
/// hysteretic active-set logic, and clears opening open-gap warm-start
/// multipliers before samples reach the descriptor.
class ChSdfContactManifoldManager {
  public:
    struct Settings {
        double dt = 0.001;
        double gap_on = 1.0e-5;
        double gap_off = 2.0e-3;
        double lambda_hold_force = 1.0e-10;
        double patch_overlap_threshold = 0.25;
        int release_steps = 1;
        bool use_prediction = true;
        bool cleanup_opening_gap = true;
    };

    struct Diagnostics {
        int candidate_count = 0;
        int active_count = 0;
        int patch_count = 0;
        int reused_patch_count = 0;
        int new_patch_count = 0;
        int released_sample_count = 0;
        double active_area = 0.0;
    };

    void SetSettings(const Settings& settings) {
        m_settings = settings;
        m_settings.dt = std::max(0.0, m_settings.dt);
        m_settings.gap_off = std::max(m_settings.gap_on, m_settings.gap_off);
        m_settings.lambda_hold_force = std::max(0.0, m_settings.lambda_hold_force);
        m_settings.patch_overlap_threshold = std::max(0.0, std::min(1.0, m_settings.patch_overlap_threshold));
        m_settings.release_steps = std::max(0, m_settings.release_steps);
    }

    const Diagnostics& GetDiagnostics() const {
        return m_diagnostics;
    }

    std::vector<ChSdfNcpDescriptorContact> Update(const std::vector<ChSdfNcpDescriptorContact>& candidates) {
        m_diagnostics = Diagnostics();
        m_diagnostics.candidate_count = static_cast<int>(candidates.size());

        std::vector<ChSdfNcpDescriptorContact> active;
        active.reserve(candidates.size());
        std::set<int> seen_ids;
        for (auto sample : candidates) {
            if (!PrepareSample(sample)) {
                continue;
            }

            const int id = sample.contact_id;
            if (id >= 0) {
                seen_ids.insert(id);
            }
            auto& memory = m_sample_memory[id];
            const bool was_active = memory.active;
            const bool near_closed = sample.gap <= m_settings.gap_on;
            const bool predicted_closed =
                m_settings.use_prediction && sample.gap_rate < 0.0 && sample.predicted_gap <= m_settings.gap_on;
            const bool within_hysteresis = sample.gap <= m_settings.gap_off;
            const bool held_by_force = memory.lambda_force > m_settings.lambda_hold_force;
            const bool opening_open_gap =
                m_settings.cleanup_opening_gap && sample.gap > m_settings.gap_on && sample.gap_rate >= 0.0;

            bool activate = near_closed || predicted_closed;
            if (!activate && was_active && within_hysteresis && !opening_open_gap) {
                activate = sample.gap_rate < 0.0 || held_by_force;
            }
            if (!activate && was_active && within_hysteresis && !opening_open_gap &&
                memory.release_count < m_settings.release_steps) {
                activate = true;
                memory.release_count++;
            }

            const double scale = ContactForceScale(sample);
            memory.last_gap = sample.gap;
            memory.last_gap_rate = sample.gap_rate;
            memory.area = scale;
            if (activate) {
                memory.active = true;
                memory.release_count = 0;
                m_diagnostics.active_area += scale;
                active.push_back(std::move(sample));
            } else {
                if (was_active) {
                    m_diagnostics.released_sample_count++;
                }
                memory.active = false;
                memory.release_count = 0;
                memory.lambda_n = 0.0;
                memory.lambda_force = 0.0;
            }
        }

        for (auto& item : m_sample_memory) {
            if (item.first >= 0 && seen_ids.count(item.first) == 0) {
                item.second.active = false;
                item.second.lambda_n = 0.0;
                item.second.lambda_force = 0.0;
            }
        }

        std::map<int, std::vector<size_t>> transient_patches;
        int fallback_patch_id = 0;
        for (size_t i = 0; i < active.size(); ++i) {
            const int patch_id = active[i].patch_id >= 0 ? active[i].patch_id : 100000000 + fallback_patch_id++;
            transient_patches[patch_id].push_back(i);
        }

        std::map<int, PatchMemory> next_patches;
        std::set<int> used_old_patch_ids;
        for (const auto& item : transient_patches) {
            const auto& indices = item.second;
            std::set<int> sample_ids;
            double patch_area = 0.0;
            for (size_t index : indices) {
                if (active[index].contact_id >= 0) {
                    sample_ids.insert(active[index].contact_id);
                }
                patch_area += ContactForceScale(active[index]);
            }

            int persistent_patch_id = -1;
            double best_overlap = 0.0;
            for (const auto& old_item : m_patch_memory) {
                if (used_old_patch_ids.count(old_item.first) > 0) {
                    continue;
                }
                const double overlap = PatchOverlap(sample_ids, old_item.second.sample_ids);
                if (overlap > best_overlap) {
                    best_overlap = overlap;
                    persistent_patch_id = old_item.first;
                }
            }
            if (persistent_patch_id >= 0 && best_overlap >= m_settings.patch_overlap_threshold) {
                used_old_patch_ids.insert(persistent_patch_id);
                m_diagnostics.reused_patch_count++;
            } else {
                persistent_patch_id = m_next_patch_id++;
                m_diagnostics.new_patch_count++;
            }

            PatchMemory patch;
            patch.patch_id = persistent_patch_id;
            patch.sample_ids = sample_ids;
            patch.area = patch_area;
            next_patches[persistent_patch_id] = patch;

            for (size_t index : indices) {
                active[index].patch_id = persistent_patch_id;
                active[index].patch_area = patch_area;
                if (!(active[index].area > 1.0e-30) && active[index].weight > 1.0e-30) {
                    active[index].area = active[index].weight;
                }
                if (active[index].contact_id >= 0) {
                    m_sample_memory[active[index].contact_id].patch_id = persistent_patch_id;
                }
            }
        }

        m_patch_memory = std::move(next_patches);
        m_diagnostics.active_count = static_cast<int>(active.size());
        m_diagnostics.patch_count = static_cast<int>(m_patch_memory.size());

        std::sort(active.begin(), active.end(), [](const ChSdfNcpDescriptorContact& a,
                                                   const ChSdfNcpDescriptorContact& b) {
            if (a.patch_id != b.patch_id) {
                return a.patch_id < b.patch_id;
            }
            return a.contact_id < b.contact_id;
        });
        return active;
    }

    void ObserveSolvedStates(const std::vector<ChSdfNcpDescriptorContactState>& states) {
        std::set<int> active_ids;
        for (const auto& state : states) {
            const int id = state.sample.contact_id;
            if (id < 0) {
                continue;
            }
            auto& memory = m_sample_memory[id];
            memory.active = state.active;
            memory.patch_id = state.sample.patch_id;
            memory.lambda_n = state.lambda_n;
            memory.lambda_force = state.lambda_force;
            memory.last_gap = state.sample.gap;
            memory.last_gap_rate = state.sample.gap_rate;
            memory.area = ContactForceScale(state.sample);
            active_ids.insert(id);
        }
        for (auto& item : m_sample_memory) {
            if (item.first >= 0 && active_ids.count(item.first) == 0) {
                item.second.active = false;
                item.second.lambda_n = 0.0;
                item.second.lambda_force = 0.0;
            }
        }
    }

  private:
    struct SampleMemory {
        bool active = false;
        int release_count = 0;
        int patch_id = -1;
        double lambda_n = 0.0;
        double lambda_force = 0.0;
        double last_gap = 0.0;
        double last_gap_rate = 0.0;
        double area = 0.0;
    };

    struct PatchMemory {
        int patch_id = -1;
        std::set<int> sample_ids;
        double area = 0.0;
    };

    static double ContactForceScale(const ChSdfNcpDescriptorContact& sample) {
        if (sample.area > 1.0e-30 && std::isfinite(sample.area)) {
            return sample.area;
        }
        if (sample.weight > 1.0e-30 && std::isfinite(sample.weight)) {
            return sample.weight;
        }
        return 1.0;
    }

    static ChVector3d PointVelocityAbs(const ChBody& body, const ChVector3d& point_abs) {
        return body.GetPosDt() + body.GetAngVelParent().Cross(point_abs - body.GetPos());
    }

    bool PrepareSample(ChSdfNcpDescriptorContact& sample) const {
        if (!sample.body_a || !sample.body_b) {
            return false;
        }
        ChVector3d normal = sample.normal_abs;
        const double length = normal.Length();
        if (!(length > 1.0e-14) || !std::isfinite(length)) {
            return false;
        }
        sample.normal_abs = normal / length;
        sample.weight = std::max(0.0, sample.weight);
        sample.area = std::max(0.0, sample.area);
        sample.patch_area = std::max(0.0, sample.patch_area);
        const ChVector3d va = PointVelocityAbs(*sample.body_a, sample.point_abs);
        const ChVector3d vb = PointVelocityAbs(*sample.body_b, sample.point_abs);
        sample.gap_rate = sample.normal_abs.Dot(va - vb);
        sample.predicted_gap = sample.gap + m_settings.dt * sample.gap_rate;
        return true;
    }

    static double PatchOverlap(const std::set<int>& a, const std::set<int>& b) {
        if (a.empty() || b.empty()) {
            return 0.0;
        }
        size_t intersection = 0;
        auto ia = a.begin();
        auto ib = b.begin();
        while (ia != a.end() && ib != b.end()) {
            if (*ia == *ib) {
                intersection++;
                ++ia;
                ++ib;
            } else if (*ia < *ib) {
                ++ia;
            } else {
                ++ib;
            }
        }
        const size_t union_size = a.size() + b.size() - intersection;
        return union_size > 0 ? static_cast<double>(intersection) / static_cast<double>(union_size) : 0.0;
    }

    Settings m_settings;
    Diagnostics m_diagnostics;
    int m_next_patch_id = 0;
    std::map<int, SampleMemory> m_sample_memory;
    std::map<int, PatchMemory> m_patch_memory;
};

/// Fixed-capacity SDF-NCP contact item for Chrono multibody systems.
///
/// The contact generator is called during Update/descriptor preparation and
/// should return the active patch samples for the current body state. Inactive
/// slots are disabled, so the item can be used in systems where the number of
/// active patch samples changes over time without changing the item's maximum
/// constraint capacity.
class ChSdfNcpConstraintContactSet : public ChPhysicsItem {
  public:
    using ContactGenerator = std::function<std::vector<ChSdfNcpDescriptorContact>()>;

    explicit ChSdfNcpConstraintContactSet(size_t max_contacts = 64) {
        SetMaxContacts(max_contacts);
    }

    virtual ChSdfNcpConstraintContactSet* Clone() const override {
        return new ChSdfNcpConstraintContactSet(*this);
    }

    void SetMaxContacts(size_t max_contacts) {
        if (max_contacts == 0) {
            throw std::invalid_argument("ChSdfNcpConstraintContactSet max_contacts must be positive.");
        }
        m_slots.resize(max_contacts);
        for (auto& slot : m_slots) {
            slot.constraint.SetMode(ChConstraint::Mode::UNILATERAL);
            slot.constraint.SetDisabled(true);
            slot.constraint.SetRightHandSide(0.0);
        }
    }

    size_t GetMaxContacts() const {
        return m_slots.size();
    }

    void SetContactGenerator(ContactGenerator generator) {
        m_generator = std::move(generator);
    }

    void SetSmoothingEps(double eps) {
        m_eps = RegularizedSmoothingEps(eps);
    }

    double GetSmoothingEps() const {
        return m_eps;
    }

    void SetRecoveryClamp(double clamp) {
        m_recovery_clamp = std::max(0.0, clamp);
    }

    double GetRecoveryClamp() const {
        return m_recovery_clamp;
    }

    /// Enable or disable the predictive active-set filter.
    ///
    /// Candidate samples can be generated over a broad band for patch continuity,
    /// but only samples satisfying the policy enter the descriptor as unilateral
    /// constraints.  The filter is intentionally generic: it uses the contact gap,
    /// normal relative velocity, a one-step predicted gap, and persistent
    /// contact_id hysteresis.
    void SetActiveSetPolicy(bool enabled,
                            double step_size,
                            double gap_on,
                            double gap_off,
                            bool use_prediction = true) {
        m_active_set_enabled = enabled;
        m_active_step_size = std::max(0.0, step_size);
        m_active_gap_on = gap_on;
        m_active_gap_off = std::max(gap_on, gap_off);
        m_active_use_prediction = use_prediction;
    }

    bool IsActiveSetPolicyEnabled() const {
        return m_active_set_enabled;
    }

    /// Set a diagonal pressure compliance for local sample intensities.
    ///
    /// The descriptor row is scaled by the sample force scale, so the actual CFM
    /// term is scale * pressure_compliance.  This regularizes p_i in
    /// force_i = p_i * area_i * n_i without changing the geometric contact
    /// frontend.
    void SetPressureCompliance(double compliance) {
        m_pressure_compliance = std::max(0.0, compliance);
    }

    double GetPressureCompliance() const {
        return m_pressure_compliance;
    }

    /// Enable patch-level pressure aggregation.
    ///
    /// When enabled, generated samples are first filtered by the active-set
    /// policy and then grouped by patch_id. Each patch contributes one pressure
    /// unknown with a preassembled Jacobian row sum_i area_i J_i. This is a
    /// generic local patch pressure mode; it does not depend on a particular
    /// benchmark geometry.
    void SetPatchPressureAggregation(bool enabled) {
        m_patch_pressure_aggregation = enabled;
    }

    bool IsPatchPressureAggregationEnabled() const {
        return m_patch_pressure_aggregation;
    }

    /// Enable row-wise patch-aware multiplier projection inside the descriptor solver.
    ///
    /// This keeps every sample as a standard unilateral descriptor row, so joint
    /// and driver constraints remain coupled through Chrono's global solver.  The
    /// projection step only regularizes pressure distribution inside connected
    /// SDF patches.
    void SetPatchPressureProjection(bool enabled, double strength, double radius, bool block_projection = true) {
        m_patch_projection_enabled = enabled;
        m_patch_projection_strength = std::max(0.0, std::min(1.0, strength));
        m_patch_projection_radius = std::max(0.0, radius);
        m_patch_block_projection_enabled = block_projection;
    }

    /// Add patch graph Laplacian regularization to descriptor block solves.
    ///
    /// This keeps sample-level nonnegative pressure unknowns in Chrono's global
    /// descriptor solve while mildly coupling neighboring pressures inside the
    /// same persistent patch.  Setting the value to zero recovers the pure
    /// Delassus block LCP W = J M^{-1} J^T + CFM.
    void SetPatchBlockLaplacianCompliance(double compliance) {
        m_patch_block_laplacian_compliance = std::max(0.0, compliance);
    }

    /// Add a generic patch-level torque regularization to the descriptor block.
    ///
    /// The regularization is optional and acts only inside one persistent patch.
    /// It keeps all sample pressures nonnegative while discouraging pressure
    /// distributions that create a large free torque on the contacted body.  A
    /// value of zero disables this term.
    void SetPatchBlockWrenchClosureStrength(double strength) {
        m_patch_block_wrench_closure_strength = std::max(0.0, strength);
    }

    /// Use a velocity-level Signorini residual in Chrono's descriptor RHS.
    ///
    /// When enabled, each active sample row represents
    ///
    ///   0 <= p_i perpendicular v_n,i + beta min(g_i,0) / dt + W p >= 0
    ///
    /// where p_i is the local sample pressure/intensity and W is supplied by the
    /// descriptor solve or by the patch Delassus block projection.
    void SetDescriptorVelocityLevelNcp(bool enabled,
                                       double step_size,
                                       double baumgarte = 0.2,
                                       double rhs_scale = 1.0) {
        m_descriptor_velocity_ncp_enabled = enabled;
        m_descriptor_velocity_step = std::max(step_size, 1.0e-12);
        m_descriptor_velocity_baumgarte = std::max(0.0, baumgarte);
        m_descriptor_velocity_rhs_scale = std::max(0.0, rhs_scale);
    }

    size_t GetNumActiveContacts() const {
        size_t count = 0;
        for (const auto& slot : m_slots) {
            if (slot.state.active) {
                count++;
            }
        }
        return count;
    }

    const std::vector<ChSdfNcpDescriptorContactState>& GetContactStates() const {
        return m_public_states;
    }

    double MaxPenetration() const {
        double value = 0.0;
        for (const auto& state : m_public_states) {
            value = std::max(value, state.penetration);
        }
        return value;
    }

    double MaxLambdaForce() const {
        double value = 0.0;
        for (const auto& state : m_public_states) {
            value = std::max(value, state.lambda_force);
        }
        return value;
    }

    double NcpResidualNorm() const {
        double sum = 0.0;
        for (const auto& state : m_public_states) {
            sum += state.ncp_residual * state.ncp_residual;
        }
        return std::sqrt(sum);
    }

    double MaxComplementarityError() const {
        double value = 0.0;
        for (const auto& state : m_public_states) {
            value = std::max(value, state.complementarity_error);
        }
        return value;
    }

    virtual unsigned int GetNumConstraintsUnilateral() override {
        return static_cast<unsigned int>(m_slots.size());
    }

    virtual void Update(double time, UpdateFlags update_flags) override {
        ChPhysicsItem::Update(time, update_flags);
        RefreshContacts();
    }

    virtual void IntStateGatherReactions(const unsigned int off_L, ChVectorDynamic<>& L) override {
        unsigned int cnt = 0;
        for (auto& slot : m_slots) {
            if (slot.constraint.IsActive()) {
                L(off_L + cnt) = slot.constraint.GetLagrangeMultiplier();
                cnt++;
            }
        }
    }

    virtual void IntStateScatterReactions(const unsigned int off_L, const ChVectorDynamic<>& L) override {
        unsigned int cnt = 0;
        for (auto& slot : m_slots) {
            if (slot.constraint.IsActive()) {
                slot.constraint.SetLagrangeMultiplier(L(off_L + cnt));
                cnt++;
            }
        }
        UpdatePublicStatesFromMultipliers(1.0);
    }

    virtual void IntLoadResidual_CqL(const unsigned int off_L,
                                     ChVectorDynamic<>& R,
                                     const ChVectorDynamic<>& L,
                                     const double c) override {
        unsigned int cnt = 0;
        for (auto& slot : m_slots) {
            if (slot.constraint.IsActive()) {
                slot.constraint.AddJacobianTransposedTimesScalarInto(R, L(off_L + cnt) * c);
                cnt++;
            }
        }
    }

    virtual void IntLoadConstraint_C(const unsigned int off_L,
                                     ChVectorDynamic<>& Qc,
                                     const double c,
                                     bool do_clamp,
                                     double recovery_clamp) override {
        unsigned int cnt = 0;
        for (auto& slot : m_slots) {
            if (slot.constraint.IsActive()) {
                const double scale = ContactForceScale(slot.state.sample);
                const double rhs = ConstraintCorrectionRhs(slot.state.sample, scale, c, recovery_clamp, do_clamp);
                Qc(off_L + cnt) += rhs;
                cnt++;
            }
        }
    }

    virtual void IntToDescriptor(const unsigned int off_v,
                                 const ChStateDelta& v,
                                 const ChVectorDynamic<>& R,
                                 const unsigned int off_L,
                                 const ChVectorDynamic<>& L,
                                 const ChVectorDynamic<>& Qc) override {
        unsigned int cnt = 0;
        for (auto& slot : m_slots) {
            if (slot.constraint.IsActive()) {
                slot.constraint.SetLagrangeMultiplier(L(off_L + cnt));
                slot.constraint.SetRightHandSide(Qc(off_L + cnt));
                cnt++;
            }
        }
    }

    virtual void IntFromDescriptor(const unsigned int off_v,
                                   ChStateDelta& v,
                                   const unsigned int off_L,
                                   ChVectorDynamic<>& L) override {
        unsigned int cnt = 0;
        for (auto& slot : m_slots) {
            if (slot.constraint.IsActive()) {
                L(off_L + cnt) = slot.constraint.GetLagrangeMultiplier();
                cnt++;
            }
        }
        UpdatePublicStatesFromMultipliers(1.0);
    }

    virtual void InjectConstraints(ChSystemDescriptor& descriptor) override {
        RefreshContacts();
        for (auto& slot : m_slots) {
            if (slot.constraint.IsActive()) {
                descriptor.InsertConstraint(&slot.constraint);
            }
        }
    }

    virtual void ConstraintsBiReset() override {
        for (auto& slot : m_slots) {
            slot.constraint.SetRightHandSide(0.0);
        }
    }

    virtual void ConstraintsBiLoad_C(double factor = 1, double recovery_clamp = 0.1, bool do_clamp = false) override {
        for (auto& slot : m_slots) {
            if (slot.constraint.IsActive()) {
                const double scale = ContactForceScale(slot.state.sample);
                const double rhs =
                    ConstraintCorrectionRhs(slot.state.sample, scale, factor, recovery_clamp, do_clamp);
                slot.constraint.SetRightHandSide(slot.constraint.GetRightHandSide() + rhs);
            }
        }
    }

    virtual void LoadConstraintJacobians() override {
        RefreshContacts();
    }

    virtual void ConstraintsFetch_react(double factor = 1) override {
        m_last_reaction_factor = factor;
        UpdatePublicStatesFromMultipliers(factor);
    }

  private:
    struct Slot {
        ChSdfNcpPatchProjectedConstraint constraint;
        ChSdfNcpDescriptorContactState state;
    };

    static ChVector3d BodyLocalTorqueRow(const ChBody& body, const ChVector3d& point_abs, const ChVector3d& force_dir) {
        const ChVector3d r_abs = point_abs - body.GetPos();
        const ChVector3d torque_abs = r_abs.Cross(force_dir);
        return body.GetRotMat().transpose() * torque_abs;
    }

    static double ContactForceScale(const ChSdfNcpDescriptorContact& sample) {
        if (sample.area > 1.0e-30 && std::isfinite(sample.area)) {
            return sample.area;
        }
        if (sample.weight > 1.0e-30 && std::isfinite(sample.weight)) {
            return sample.weight;
        }
        return 1.0;
    }

    double ConstraintCorrectionRhs(const ChSdfNcpDescriptorContact& sample,
                                   double scale,
                                   double factor,
                                   double recovery_clamp,
                                   bool do_clamp) const {
        if (m_descriptor_velocity_ncp_enabled) {
            double rhs = scale * m_descriptor_velocity_rhs_scale * m_descriptor_velocity_baumgarte *
                         std::min(0.0, sample.gap) / std::max(m_descriptor_velocity_step, 1.0e-12);
            if (do_clamp) {
                rhs = std::max(rhs, -scale * recovery_clamp);
            }
            return rhs;
        }

        double rhs = factor * scale * sample.gap;
        if (do_clamp) {
            rhs = std::max(rhs, -scale * recovery_clamp);
        }
        return rhs;
    }

    static std::pair<ChVector3d, ChVector3d> ContactForceTorquePerPressure(const ChSdfNcpDescriptorContact& sample,
                                                                            double scale) {
        if (sample.use_custom_wrench) {
            return {sample.force_on_body_a_per_lambda_abs, sample.torque_on_body_a_per_lambda_abs};
        }
        const ChVector3d force = sample.normal_abs * scale;
        const ChVector3d torque =
            sample.body_a ? (sample.point_abs - sample.body_a->GetPos()).Cross(force) : ChVector3d(0, 0, 0);
        return {force, torque};
    }

    static ChVector3d PointVelocityAbs(const ChBody& body, const ChVector3d& point_abs) {
        return body.GetPosDt() + body.GetAngVelParent().Cross(point_abs - body.GetPos());
    }

    static double NormalGapRate(const ChSdfNcpDescriptorContact& sample, const ChVector3d& normal) {
        if (!sample.body_a || !sample.body_b) {
            return 0.0;
        }
        const ChVector3d va = PointVelocityAbs(*sample.body_a, sample.point_abs);
        const ChVector3d vb = PointVelocityAbs(*sample.body_b, sample.point_abs);
        return normal.Dot(va - vb);
    }

    bool PrepareSample(ChSdfNcpDescriptorContact& sample) const {
        if (!sample.body_a || !sample.body_b) {
            return false;
        }
        ChVector3d normal = sample.normal_abs;
        const double normal_length = normal.Length();
        if (!(normal_length > 1.0e-14) || !std::isfinite(normal_length)) {
            return false;
        }
        normal /= normal_length;
        sample.normal_abs = normal;
        sample.weight = std::max(0.0, sample.weight);
        sample.area = std::max(0.0, sample.area);
        sample.patch_area = std::max(0.0, sample.patch_area);
        sample.gap_rate = NormalGapRate(sample, normal);
        sample.predicted_gap = sample.gap + m_active_step_size * sample.gap_rate;
        return true;
    }

    bool ShouldActivate(const ChSdfNcpDescriptorContact& sample) const {
        if (!m_active_set_enabled) {
            return true;
        }
        const bool persistent = sample.contact_id >= 0 && m_persistent_active_ids.count(sample.contact_id) > 0;
        const bool near_closed = sample.gap <= m_active_gap_on;
        const bool predicted_closed =
            m_active_use_prediction && m_active_step_size > 0.0 && sample.gap_rate < 0.0 &&
            sample.predicted_gap <= m_active_gap_on;
        const bool retained = persistent && sample.gap <= m_active_gap_off && sample.gap_rate < 0.0;
        return near_closed || predicted_closed || retained;
    }

    std::vector<ChSdfNcpDescriptorContact> AggregatePatchContacts(std::vector<ChSdfNcpDescriptorContact> generated,
                                                                  std::set<int>& next_active_ids) const {
        struct PatchAccumulator {
            bool initialized = false;
            std::shared_ptr<ChBody> body_a;
            std::shared_ptr<ChBody> body_b;
            int patch_id = -1;
            int min_gap_contact_id = -1;
            double min_gap = std::numeric_limits<double>::max();
            double min_gap_rate = 0.0;
            double min_predicted_gap = 0.0;
            double total_area = 0.0;
            double total_weight = 0.0;
            double center_weight = 0.0;
            double gap_sum = 0.0;
            double gap_rate_sum = 0.0;
            double predicted_gap_sum = 0.0;
            ChVector3d center = ChVector3d(0, 0, 0);
            ChVector3d normal_sum = ChVector3d(0, 0, 0);
            ChVector3d Cq_a_trans = ChVector3d(0, 0, 0);
            ChVector3d Cq_a_rot_local = ChVector3d(0, 0, 0);
            ChVector3d Cq_b_trans = ChVector3d(0, 0, 0);
            ChVector3d Cq_b_rot_local = ChVector3d(0, 0, 0);
            ChVector3d torque_a_abs = ChVector3d(0, 0, 0);
            ChVector3d torque_b_abs = ChVector3d(0, 0, 0);
        };

        std::map<int, PatchAccumulator> patches;
        int fallback_patch_id = 0;
        for (auto& sample : generated) {
            if (!PrepareSample(sample) || !ShouldActivate(sample)) {
                continue;
            }
            if (sample.contact_id >= 0) {
                next_active_ids.insert(sample.contact_id);
            }
            const int patch_id = sample.patch_id >= 0 ? sample.patch_id : 100000000 + fallback_patch_id++;
            const double scale = ContactForceScale(sample);
            const ChVector3d force_row_a = sample.normal_abs * scale;
            const ChVector3d force_row_b = -force_row_a;
            const ChVector3d torque_a_local = BodyLocalTorqueRow(*sample.body_a, sample.point_abs, force_row_a);
            const ChVector3d torque_b_local = BodyLocalTorqueRow(*sample.body_b, sample.point_abs, force_row_b);

            auto& patch = patches[patch_id];
            if (!patch.initialized) {
                patch.initialized = true;
                patch.body_a = sample.body_a;
                patch.body_b = sample.body_b;
                patch.patch_id = patch_id;
                patch.min_gap_contact_id = sample.contact_id;
            }
            if (sample.gap < patch.min_gap) {
                patch.min_gap = sample.gap;
                patch.min_gap_rate = sample.gap_rate;
                patch.min_predicted_gap = sample.predicted_gap;
                patch.min_gap_contact_id = sample.contact_id;
            }
            patch.total_area += std::max(0.0, sample.area);
            patch.total_weight += std::max(0.0, sample.weight);
            patch.center += sample.point_abs * scale;
            patch.center_weight += scale;
            patch.gap_sum += sample.gap * scale;
            patch.gap_rate_sum += sample.gap_rate * scale;
            patch.predicted_gap_sum += sample.predicted_gap * scale;
            patch.normal_sum += force_row_a;
            patch.Cq_a_trans += force_row_a;
            patch.Cq_a_rot_local += torque_a_local;
            patch.Cq_b_trans += force_row_b;
            patch.Cq_b_rot_local += torque_b_local;
            patch.torque_a_abs += sample.body_a->GetRotMat() * torque_a_local;
            patch.torque_b_abs += sample.body_b->GetRotMat() * torque_b_local;
        }

        std::vector<ChSdfNcpDescriptorContact> aggregated;
        aggregated.reserve(patches.size());
        for (const auto& item : patches) {
            const auto& patch = item.second;
            if (!patch.initialized || !patch.body_a || !patch.body_b) {
                continue;
            }
            const double force_norm = patch.Cq_a_trans.Length();
            ChSdfNcpDescriptorContact contact;
            contact.body_a = patch.body_a;
            contact.body_b = patch.body_b;
            contact.point_abs = patch.center_weight > 1.0e-30 ? patch.center / patch.center_weight : ChVector3d(0, 0, 0);
            contact.normal_abs = force_norm > 1.0e-30 ? patch.Cq_a_trans / force_norm : ChVector3d(0, 1, 0);
            contact.gap = patch.center_weight > 1.0e-30 ? patch.gap_sum / patch.center_weight : patch.min_gap;
            contact.weight = patch.total_weight > 1.0e-30 ? patch.total_weight : 1.0;
            contact.area = patch.total_area > 1.0e-30 ? patch.total_area : patch.total_weight;
            contact.patch_area = contact.area;
            contact.gap_rate = patch.center_weight > 1.0e-30 ? patch.gap_rate_sum / patch.center_weight : patch.min_gap_rate;
            contact.predicted_gap =
                patch.center_weight > 1.0e-30 ? patch.predicted_gap_sum / patch.center_weight : patch.min_predicted_gap;
            contact.contact_id = patch.min_gap_contact_id;
            contact.patch_id = patch.patch_id;
            contact.use_custom_jacobian = true;
            contact.custom_Cq_a_trans_abs = patch.Cq_a_trans;
            contact.custom_Cq_a_rot_local = patch.Cq_a_rot_local;
            contact.custom_Cq_b_trans_abs = patch.Cq_b_trans;
            contact.custom_Cq_b_rot_local = patch.Cq_b_rot_local;
            contact.use_custom_wrench = true;
            contact.force_on_body_a_per_lambda_abs = patch.Cq_a_trans;
            contact.torque_on_body_a_per_lambda_abs = patch.torque_a_abs;
            contact.force_on_body_b_per_lambda_abs = patch.Cq_b_trans;
            contact.torque_on_body_b_per_lambda_abs = patch.torque_b_abs;
            aggregated.push_back(contact);
        }
        return aggregated;
    }

    void LoadJacobian(Slot& slot) {
        auto& sample = slot.state.sample;
        if (!PrepareSample(sample)) {
            slot.constraint.SetDisabled(true);
            slot.state.active = false;
            return;
        }

        if (!ShouldActivate(sample)) {
            slot.constraint.SetDisabled(true);
            slot.constraint.SetRightHandSide(0.0);
            slot.state.active = false;
            slot.state.penetration = std::max(0.0, -sample.gap);
            slot.state.lambda_n = 0.0;
            slot.state.lambda_force = 0.0;
            slot.state.ncp_residual = 0.0;
            slot.state.complementarity_error = ComplementarityError(sample.gap, 0.0);
            return;
        }
        const double force_scale = ContactForceScale(sample);

        slot.constraint.SetVariables(&sample.body_a->Variables(), &sample.body_b->Variables());
        slot.constraint.SetMode(ChConstraint::Mode::UNILATERAL);
        slot.constraint.SetDisabled(false);
        slot.constraint.SetComplianceTerm(force_scale * m_pressure_compliance);

        auto Cq_a = slot.constraint.Get_Cq_a();
        auto Cq_b = slot.constraint.Get_Cq_b();
        Cq_a.setZero();
        Cq_b.setZero();

        if (sample.use_custom_jacobian) {
            Cq_a(0) = sample.custom_Cq_a_trans_abs.x();
            Cq_a(1) = sample.custom_Cq_a_trans_abs.y();
            Cq_a(2) = sample.custom_Cq_a_trans_abs.z();
            Cq_a(3) = sample.custom_Cq_a_rot_local.x();
            Cq_a(4) = sample.custom_Cq_a_rot_local.y();
            Cq_a(5) = sample.custom_Cq_a_rot_local.z();

            Cq_b(0) = sample.custom_Cq_b_trans_abs.x();
            Cq_b(1) = sample.custom_Cq_b_trans_abs.y();
            Cq_b(2) = sample.custom_Cq_b_trans_abs.z();
            Cq_b(3) = sample.custom_Cq_b_rot_local.x();
            Cq_b(4) = sample.custom_Cq_b_rot_local.y();
            Cq_b(5) = sample.custom_Cq_b_rot_local.z();
        } else {
            const ChVector3d scaled_normal = sample.normal_abs * force_scale;
            const ChVector3d torque_a = BodyLocalTorqueRow(*sample.body_a, sample.point_abs, scaled_normal);
            const ChVector3d torque_b = BodyLocalTorqueRow(*sample.body_b, sample.point_abs, -scaled_normal);

            Cq_a(0) = scaled_normal.x();
            Cq_a(1) = scaled_normal.y();
            Cq_a(2) = scaled_normal.z();
            Cq_a(3) = torque_a.x();
            Cq_a(4) = torque_a.y();
            Cq_a(5) = torque_a.z();

            Cq_b(0) = -scaled_normal.x();
            Cq_b(1) = -scaled_normal.y();
            Cq_b(2) = -scaled_normal.z();
            Cq_b(3) = torque_b.x();
            Cq_b(4) = torque_b.y();
            Cq_b(5) = torque_b.z();
        }

        slot.state.active = true;
        slot.state.penetration = std::max(0.0, -sample.gap);
        slot.state.lambda_n = std::max(0.0, slot.constraint.GetLagrangeMultiplier());
        slot.state.lambda_force = slot.state.lambda_n * force_scale;
        slot.state.ncp_residual = std::abs(SmoothFischerBurmeister(sample.gap, slot.state.lambda_n, 1.0e-12));
        slot.state.complementarity_error = ComplementarityError(sample.gap, slot.state.lambda_n);
    }

    void RefreshContacts() {
        if (!m_generator) {
            DisableAllSlots();
            return;
        }

        std::map<int, double> previous_lambdas;
        for (const auto& slot : m_slots) {
            if (slot.state.sample.contact_id >= 0) {
                previous_lambdas[slot.state.sample.contact_id] = slot.constraint.GetLagrangeMultiplier();
            }
        }

        std::vector<ChSdfNcpDescriptorContact> generated = m_generator();
        std::set<int> next_active_ids;
        if (m_patch_pressure_aggregation) {
            generated = AggregatePatchContacts(std::move(generated), next_active_ids);
        }
        if (generated.size() > m_slots.size()) {
            generated.resize(m_slots.size());
        }

        for (size_t i = 0; i < m_slots.size(); i++) {
            auto& slot = m_slots[i];
            if (i < generated.size()) {
                slot.state.sample = std::move(generated[i]);
                double warm_lambda = 0.0;
                const int contact_id = slot.state.sample.contact_id;
                if (contact_id >= 0) {
                    const auto found = previous_lambdas.find(contact_id);
                    if (found != previous_lambdas.end()) {
                        warm_lambda = std::max(0.0, found->second);
                    }
                }
                slot.constraint.SetLagrangeMultiplier(warm_lambda);
                LoadJacobian(slot);
                if (!m_patch_pressure_aggregation && slot.state.active && slot.state.sample.contact_id >= 0) {
                    next_active_ids.insert(slot.state.sample.contact_id);
                }
            } else {
                slot.constraint.SetDisabled(true);
                slot.constraint.SetRightHandSide(0.0);
                slot.state = ChSdfNcpDescriptorContactState();
            }
        }
        ConfigurePatchProjectionGroups();
        m_persistent_active_ids = std::move(next_active_ids);
        UpdatePublicStatesFromMultipliers(m_last_reaction_factor);
    }

    void ConfigurePatchProjectionGroups() {
        std::map<int, std::vector<ChSdfNcpPatchProjectedConstraint*>> patch_groups;
        double area_sum = 0.0;
        int area_count = 0;
        for (auto& slot : m_slots) {
            slot.constraint.ConfigurePatchProjection(false, -1, ChVector3d(0, 0, 0), 0.0, 0.0, 0.0);
            if (!m_patch_projection_enabled || !slot.state.active || slot.state.sample.patch_id < 0) {
                continue;
            }
            const double area = ContactForceScale(slot.state.sample);
            const auto force_torque = ContactForceTorquePerPressure(slot.state.sample, area);
            area_sum += area;
            area_count++;
            slot.constraint.ConfigurePatchProjection(true,
                                                      slot.state.sample.patch_id,
                                                      slot.state.sample.point_abs,
                                                      area,
                                                      m_patch_projection_strength,
                                                      m_patch_projection_radius,
                                                      force_torque.first,
                                                      force_torque.second);
            slot.constraint.SetPatchBlockProjection(m_patch_block_projection_enabled,
                                                    m_patch_block_laplacian_compliance,
                                                    m_patch_block_wrench_closure_strength);
            patch_groups[slot.state.sample.patch_id].push_back(&slot.constraint);
        }

        if (m_patch_projection_enabled && m_patch_projection_radius <= 0.0 && area_count > 0) {
            const double radius = 2.0 * std::sqrt(std::max(area_sum / static_cast<double>(area_count), 1.0e-12));
            for (auto& slot : m_slots) {
                if (!slot.state.active || slot.state.sample.patch_id < 0) {
                    continue;
                }
                const double area = ContactForceScale(slot.state.sample);
                const auto force_torque = ContactForceTorquePerPressure(slot.state.sample, area);
                slot.constraint.ConfigurePatchProjection(true,
                                                          slot.state.sample.patch_id,
                                                          slot.state.sample.point_abs,
                                                          area,
                                                          m_patch_projection_strength,
                                                          radius,
                                                          force_torque.first,
                                                          force_torque.second);
                slot.constraint.SetPatchBlockProjection(m_patch_block_projection_enabled,
                                                        m_patch_block_laplacian_compliance,
                                                        m_patch_block_wrench_closure_strength);
            }
        }

        for (auto& item : patch_groups) {
            for (auto* constraint : item.second) {
                if (constraint) {
                    constraint->SetPatchNeighbors(item.second);
                }
            }
        }
    }

    void DisableAllSlots() {
        for (auto& slot : m_slots) {
            slot.constraint.SetDisabled(true);
            slot.constraint.SetRightHandSide(0.0);
            slot.state = ChSdfNcpDescriptorContactState();
        }
        m_public_states.clear();
    }

    void UpdatePublicStatesFromMultipliers(double factor) {
        m_public_states.clear();
        for (auto& slot : m_slots) {
            if (!slot.constraint.IsActive()) {
                continue;
            }
            slot.state.lambda_n = std::max(0.0, slot.constraint.GetLagrangeMultiplier() * factor);
            slot.state.lambda_force = slot.state.lambda_n * ContactForceScale(slot.state.sample);
            slot.state.penetration = std::max(0.0, -slot.state.sample.gap);
            slot.state.ncp_residual =
                std::abs(SmoothFischerBurmeister(slot.state.sample.gap, slot.state.lambda_n, m_eps));
            slot.state.complementarity_error = ComplementarityError(slot.state.sample.gap, slot.state.lambda_n);
            m_public_states.push_back(slot.state);
        }
    }

    std::vector<Slot> m_slots;
    std::vector<ChSdfNcpDescriptorContactState> m_public_states;
    ContactGenerator m_generator;
    double m_eps = 1.0e-6;
    double m_recovery_clamp = 0.1;
    double m_last_reaction_factor = 1.0;
    bool m_active_set_enabled = false;
    bool m_active_use_prediction = true;
    double m_active_step_size = 0.0;
    double m_active_gap_on = 0.0;
    double m_active_gap_off = 0.0;
    double m_pressure_compliance = 0.0;
    bool m_patch_pressure_aggregation = false;
    bool m_patch_projection_enabled = false;
    bool m_patch_block_projection_enabled = true;
    double m_patch_projection_strength = 0.0;
    double m_patch_projection_radius = 0.0;
    double m_patch_block_laplacian_compliance = 0.0;
    double m_patch_block_wrench_closure_strength = 0.0;
    bool m_descriptor_velocity_ncp_enabled = false;
    double m_descriptor_velocity_step = 1.0e-3;
    double m_descriptor_velocity_baumgarte = 0.2;
    double m_descriptor_velocity_rhs_scale = 1.0;
    std::set<int> m_persistent_active_ids;
};

}  // namespace sdfncp
}  // namespace chrono

#endif
