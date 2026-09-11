use noiselab::{run_stage_u1_panel, StageU1PanelConfig};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let config = StageU1PanelConfig::default();
    let result = run_stage_u1_panel(&config)?;

    println!(
        "left\tright\tpurpose\tfine_distance\tterminal_distance\tconvergence\tscaling_gap\tempirical_p\tevidence"
    );
    for comparison in result.comparisons {
        println!(
            "{:?}\t{:?}\t{:?}\t{:.12}\t{:.12}\t{:.12}\t{:.12}\t{:.12}\t{:?}",
            comparison.left,
            comparison.right,
            comparison.purpose,
            comparison.fine_scale_distance,
            comparison.terminal_scale_distance,
            comparison.convergence_score,
            comparison.variance_scaling_exponent_gap,
            comparison.empirical_p_value,
            comparison.evidence,
        );
    }

    Ok(())
}
