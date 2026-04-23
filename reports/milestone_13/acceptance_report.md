# Milestone 13 Acceptance Report

## Executive Summary

Milestone 13 implemented adaptive carry-over v1 to replace fixed carry-over (alpha=constant) with context-aware alpha. The implementation is **functional** but the success criteria are **not fully met**.

## Implementation Status

| Task | Status | Notes |
|------|--------|-------|
| Implement adaptive alpha formula | COMPLETE | alpha = f(theta, match, topology) |
| Add topology event detection | COMPLETE | Born/Drifted/Stable |
| Integrate with existing transport system | COMPLETE | Preserves all switches |
| Create experimental matrix | COMPLETE | 5 configs across 2 cases |
| Verify A-only baseline | PASS | Matches expected |
| Verify C-adaptive improves Case A | FAIL | y_error worsened 49% |
| Verify C-adaptive not worse on Case B | PASS | y_error improved 51% |

## Bug Fix

During development, a critical bug was discovered and fixed:
- **Bug**: `topology_birth_event` was copied from previous track without reset, causing alpha to always be 0 for matched tracks
- **Fix**: Reset `topology_birth_event = 0` when track status is `TrackStatus::Alive`
- **Impact**: Without this fix, adaptive alpha was effectively disabled

## Performance Summary

### Case A (Rotational Contact)
- A-only: y_error = 0.0632m (baseline)
- C-adaptive: y_error = 0.0941m (**49% worse**)
- C-fixed: y_error = 0.0750m (reference)

### Case B (Multi-Patch)
- A-only: y_error = 0.0639m (baseline)
- C-adaptive: y_error = 0.0314m (**51% better**)
- C-fixed: y_error = 0.0740m (reference)

## Parameter Analysis

| Parameter | Value | Effect |
|-----------|-------|--------|
| alpha_base | 0.8 | Too aggressive, results in avg alpha ~0.78 |
| theta_scale | 0.3 rad | Too large, theta ~0.001 causes negligible decay |
| match_scale | 0.05 m | Too large, match_dist ~0.001 causes negligible decay |

The exponential decay functions are too gentle for the small theta and match distances observed, resulting in adaptive alpha being nearly constant (~0.78) instead of varying meaningfully.

## Recommendations for Next Iteration

1. **Reduce alpha_base** to 0.4-0.5 for more conservative behavior
2. **Decrease theta_scale** to 0.05-0.1 for meaningful sensitivity
3. **Consider match distance thresholds** for binary carry-over disable
4. **Case-aware parameters** could resolve the Case A/B conflict

## Verdict

**Partial Success**: Adaptive carry-over works correctly and significantly improves Case B, but degrades Case A. The fundamental tradeoff between rotational contact and multi-patch stability remains unresolved. The adaptive mechanism itself is sound, but the parameter values need tuning for Case A scenarios.
