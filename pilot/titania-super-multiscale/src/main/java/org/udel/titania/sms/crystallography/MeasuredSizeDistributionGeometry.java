package org.udel.titania.sms.crystallography;

import java.util.List;

/** Future input: a measured size distribution; results are mass-weighted over its components. */
public record MeasuredSizeDistributionGeometry(List<Component> components) {
    public record Component(double diameterNm, double massFraction) {}
    public MeasuredSizeDistributionGeometry {
        double s = components.stream().mapToDouble(Component::massFraction).sum();
        if (Math.abs(s - 1) > 1e-9) throw new IllegalArgumentException("mass fractions must sum to 1");
    }
    public static String status() { return "MISSING_INPUT: no measured size distribution for A600-R1100"; }
}
