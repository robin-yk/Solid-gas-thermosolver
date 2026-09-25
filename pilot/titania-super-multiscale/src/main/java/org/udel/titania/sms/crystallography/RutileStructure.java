package org.udel.titania.sms.crystallography;

import java.util.EnumMap;
import java.util.Map;
import org.udel.titania.sms.core.SiteType;

/**
 * Questions 6-7 start here: the rutile cell projected on [110]. Oxygen heights, which O belongs to which
 * site class, and how many of each there are per (110) 1x1 cell are all counted from atomic positions.
 */
public final class RutileStructure {
    public final double a, c, u;
    private final Map<SiteType, Integer> mult = new EnumMap<>(SiteType.class);
    private final Map<SiteType, Double> offset = new EnumMap<>(SiteType.class);

    public RutileStructure(double aNm, double cNm, double u) {
        this.a = aNm; this.c = cNm; this.u = u;
        double[][] ti = {{0, 0, 0}, {0.5, 0.5, 0.5}};
        double[][] o = {{u, u, 0}, {1 - u, 1 - u, 0}, {0.5 + u, 0.5 - u, 0.5}, {0.5 - u, 0.5 + u, 0.5}};
        double d = d110();
        // every Ti must sit on a (110) Ti plane, or the trilayer picture is wrong
        for (double[] t : ti) {
            double h = height(t), off = h - Math.round(h / d) * d;
            if (Math.abs(off) > 1e-12) throw new IllegalStateException("Ti off its (110) plane");
        }
        for (SiteType t : SiteType.values()) mult.put(t, 0);
        for (double[] p : o) {
            double h = height(p), off = h - Math.round(h / d) * d;
            SiteType t = off > 1e-9 ? SiteType.BRI : off < -1e-9 ? SiteType.SBR : SiteType.IPL;
            mult.merge(t, 1, Integer::sum);
            offset.put(t, off);
        }
    }
    private double height(double[] f) { return (f[0] + f[1]) * a / Math.sqrt(2); }

    public double d110() { return a / Math.sqrt(2); }
    /** (110) 1x1 surface cell: a*sqrt(2) by c. One trilayer of this area has the volume of one conventional cell. */
    public double cellArea110Nm2() { return a * Math.sqrt(2) * c; }
    public int multiplicity(SiteType t) { return mult.get(t); }
    public int oxygenPerCellPerTrilayer() { return mult.values().stream().mapToInt(Integer::intValue).sum(); }
    /** Height of a site class above (+) or below (-) its own Ti plane. */
    public double offsetFromTiPlane(SiteType t) { return offset.get(t); }
    /** Depth below the outermost bridging O of site class t in trilayer k (k = 1 at the surface). */
    public double depthNm(int k, SiteType t) {
        return (k - 1) * d110() + offset.get(SiteType.BRI) - offset.get(t);
    }
    public double oxygenSitesPerNm3() { return oxygenPerCellPerTrilayer() / (cellArea110Nm2() * d110()); }
    public Map<SiteType, Integer> multiplicities() { return mult; }
}
