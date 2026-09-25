package org.udel.titania.sms.crystallography;

import org.udel.titania.sms.core.Phys;
import org.udel.titania.sms.core.SiteType;

/**
 * Question 7, the core of the model. A site class in trilayer k holds (multiplicity per 1x1 cell) x
 * (cells in the surface parallel to the particle at that depth). Nothing is given a capacity by default:
 * bridging, in-plane and sub-bridging each carry their own count, and the audit checks the sum against
 * the particle's total oxygen.
 */
public final class CapacityBuilder {
    private final ParticleGeometry g;
    private final RutileStructure s;
    private final Phys phys;
    private final double molarMass;

    public CapacityBuilder(ParticleGeometry g, RutileStructure s, Phys phys, double molarMass) {
        this.g = g; this.s = s; this.phys = phys; this.molarMass = molarMass;
    }
    public ParticleGeometry geometry() { return g; }
    public RutileStructure structure() { return s; }

    /** umol O sites per gram of oxide. */
    public double capacity(int k, SiteType t) {
        double depth = s.depthNm(k, t);
        if (depth >= g.radiusNm()) return 0;
        double sites = s.multiplicity(t) / s.cellArea110Nm2() * g.areaAtDepthNm2(depth);
        return sites / g.massG() / phys.NA * 1e6;
    }
    public double surfaceBridgingCapacity() { return capacity(1, SiteType.BRI); }
    public double totalOxygen() { return 2.0 / molarMass * 1e6; }
    /** Shell capacity from volume and bulk O density, the internal form of Question 7. */
    public double shellCapacity(double innerDepthNm, double outerDepthNm) {
        return g.shellVolumeNm3(innerDepthNm, outerDepthNm) * s.oxygenSitesPerNm3() / g.massG() / phys.NA * 1e6;
    }

    public record Audit(double trilayerDensityPerNm3, double bulkDensityFromMassPerNm3, double densityRelDiff,
                        double discreteSumUmolG, double totalOxygenUmolG, double sumRelDiff, int trilayers) {}
    public Audit audit(double densityGcm3) {
        double rhoO = 2 * densityGcm3 / molarMass * phys.NA * 1e-21;
        double sum = 0; int k = 1;
        for (; s.depthNm(k, SiteType.BRI) < g.radiusNm(); k++)
            for (SiteType t : SiteType.values()) sum += capacity(k, t);
        return new Audit(s.oxygenSitesPerNm3(), rhoO, s.oxygenSitesPerNm3() / rhoO - 1,
            sum, totalOxygen(), sum / totalOxygen() - 1, k - 1);
    }
}
