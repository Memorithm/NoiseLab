//! Complete outcome-blind manifest for the preregistered U2 surrogate panel.
//!
//! Materialization expands only frozen pair/null/repetition indices and their
//! deterministic seeds. It never reads source series, surrogate outputs,
//! scores, p-values, decisions, or protected holdouts.

use crate::u2_plan::{U2ExecutionPlan, U2SurrogateJob, U2_FROZEN_NULLS};

/// Materialize every preregistered surrogate descriptor in stable order.
///
/// Ordering is pair index, then null family in [`U2_FROZEN_NULLS`] order, then
/// repetition. The returned descriptors remain outcome-blind.
#[must_use]
pub fn materialize_u2_manifest(plan: U2ExecutionPlan) -> Vec<U2SurrogateJob> {
    let mut jobs = Vec::with_capacity(plan.surrogate_job_count());
    for pair_index in 0..plan.pairs.len() {
        for null_family in U2_FROZEN_NULLS {
            for repetition in 0..plan.surrogates_per_null {
                let job = plan
                    .surrogate_job(pair_index, null_family, repetition)
                    .expect("loop bounds are inside the preregistered contract");
                jobs.push(job);
            }
        }
    }
    jobs
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use super::*;
    use crate::u2_plan::U2NullFamily;

    #[test]
    fn manifest_is_complete_and_deterministic() {
        let plan = U2ExecutionPlan::preregistered(0x1234_5678).unwrap();
        let a = materialize_u2_manifest(plan);
        let b = materialize_u2_manifest(plan);
        assert_eq!(a, b);
        assert_eq!(a.len(), plan.surrogate_job_count());
    }

    #[test]
    fn every_descriptor_identity_is_unique() {
        let plan = U2ExecutionPlan::preregistered(7).unwrap();
        let jobs = materialize_u2_manifest(plan);
        let identities = jobs
            .iter()
            .map(|job| {
                (
                    job.pair_index,
                    job.null_family as u8,
                    job.repetition,
                    job.seed,
                )
            })
            .collect::<BTreeSet<_>>();
        assert_eq!(identities.len(), jobs.len());
    }

    #[test]
    fn stable_order_starts_with_first_pair_and_shuffled_null() {
        let plan = U2ExecutionPlan::preregistered(1).unwrap();
        let jobs = materialize_u2_manifest(plan);
        let first = jobs.first().unwrap();
        assert_eq!(first.pair_index, 0);
        assert_eq!(first.null_family, U2NullFamily::ShuffledMarginal);
        assert_eq!(first.repetition, 0);
    }
}
