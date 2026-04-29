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
#include <memory>
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

/// One frictionless SDF-NCP contact sample between two Chrono rigid bodies.
///
/// The gap is defined as g = phi(x_A), where x_A is a world-space point attached
/// to body_a and phi is the SDF carried by body_b. The normal points in the
/// world direction of increasing gap. The resulting velocity Jacobian is
///
///   g_dot = n^T(v_A + w_A x r_A - v_B - w_B x r_B)
///
/// with body angular speeds represented in Chrono body-local coordinates.
struct ChSdfNcpDescriptorContact {
    std::shared_ptr<ChBody> body_a;
    std::shared_ptr<ChBody> body_b;
    ChVector3d point_abs = ChVector3d(0, 0, 0);
    ChVector3d normal_abs = ChVector3d(0, 1, 0);
    double gap = 0.0;
    double weight = 1.0;
    int contact_id = -1;
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
                double rhs = c * slot.state.sample.gap;
                if (do_clamp) {
                    rhs = std::max(rhs, -recovery_clamp);
                }
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
                double rhs = factor * slot.state.sample.gap;
                if (do_clamp) {
                    rhs = std::max(rhs, -recovery_clamp);
                }
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
        ChConstraintTwoBodies constraint;
        ChSdfNcpDescriptorContactState state;
    };

    static ChVector3d BodyLocalTorqueRow(const ChBody& body, const ChVector3d& point_abs, const ChVector3d& force_dir) {
        const ChVector3d r_abs = point_abs - body.GetPos();
        const ChVector3d torque_abs = r_abs.Cross(force_dir);
        return body.GetRotMat().transpose() * torque_abs;
    }

    static void LoadJacobian(Slot& slot) {
        auto& sample = slot.state.sample;
        if (!sample.body_a || !sample.body_b) {
            slot.constraint.SetDisabled(true);
            slot.state.active = false;
            return;
        }

        ChVector3d normal = sample.normal_abs;
        const double normal_length = normal.Length();
        if (!(normal_length > 1.0e-14) || !std::isfinite(normal_length)) {
            slot.constraint.SetDisabled(true);
            slot.state.active = false;
            return;
        }
        normal /= normal_length;
        sample.normal_abs = normal;
        sample.weight = std::max(0.0, sample.weight);

        slot.constraint.SetVariables(&sample.body_a->Variables(), &sample.body_b->Variables());
        slot.constraint.SetMode(ChConstraint::Mode::UNILATERAL);
        slot.constraint.SetDisabled(false);
        slot.constraint.SetRightHandSide(0.0);

        auto Cq_a = slot.constraint.Get_Cq_a();
        auto Cq_b = slot.constraint.Get_Cq_b();
        Cq_a.setZero();
        Cq_b.setZero();

        const ChVector3d torque_a = BodyLocalTorqueRow(*sample.body_a, sample.point_abs, normal);
        const ChVector3d torque_b = BodyLocalTorqueRow(*sample.body_b, sample.point_abs, -normal);

        Cq_a(0) = normal.x();
        Cq_a(1) = normal.y();
        Cq_a(2) = normal.z();
        Cq_a(3) = torque_a.x();
        Cq_a(4) = torque_a.y();
        Cq_a(5) = torque_a.z();

        Cq_b(0) = -normal.x();
        Cq_b(1) = -normal.y();
        Cq_b(2) = -normal.z();
        Cq_b(3) = torque_b.x();
        Cq_b(4) = torque_b.y();
        Cq_b(5) = torque_b.z();

        slot.state.active = true;
        slot.state.penetration = std::max(0.0, -sample.gap);
        slot.state.lambda_n = std::max(0.0, slot.constraint.GetLagrangeMultiplier());
        slot.state.lambda_force = slot.state.lambda_n * sample.weight;
        slot.state.ncp_residual = std::abs(SmoothFischerBurmeister(sample.gap, slot.state.lambda_force, 1.0e-12));
        slot.state.complementarity_error = ComplementarityError(sample.gap, slot.state.lambda_force);
    }

    void RefreshContacts() {
        if (!m_generator) {
            DisableAllSlots();
            return;
        }

        std::vector<ChSdfNcpDescriptorContact> generated = m_generator();
        if (generated.size() > m_slots.size()) {
            generated.resize(m_slots.size());
        }

        for (size_t i = 0; i < m_slots.size(); i++) {
            auto& slot = m_slots[i];
            if (i < generated.size()) {
                slot.state.sample = std::move(generated[i]);
                LoadJacobian(slot);
            } else {
                slot.constraint.SetDisabled(true);
                slot.constraint.SetRightHandSide(0.0);
                slot.state = ChSdfNcpDescriptorContactState();
            }
        }
        UpdatePublicStatesFromMultipliers(m_last_reaction_factor);
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
            slot.state.lambda_force = slot.state.lambda_n * slot.state.sample.weight;
            slot.state.penetration = std::max(0.0, -slot.state.sample.gap);
            slot.state.ncp_residual =
                std::abs(SmoothFischerBurmeister(slot.state.sample.gap, slot.state.lambda_force, m_eps));
            slot.state.complementarity_error = ComplementarityError(slot.state.sample.gap, slot.state.lambda_force);
            m_public_states.push_back(slot.state);
        }
    }

    std::vector<Slot> m_slots;
    std::vector<ChSdfNcpDescriptorContactState> m_public_states;
    ContactGenerator m_generator;
    double m_eps = 1.0e-6;
    double m_recovery_clamp = 0.1;
    double m_last_reaction_factor = 1.0;
};

}  // namespace sdfncp
}  // namespace chrono

#endif
