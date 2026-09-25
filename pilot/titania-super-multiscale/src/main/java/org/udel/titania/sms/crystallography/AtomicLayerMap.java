package org.udel.titania.sms.crystallography;

import java.util.ArrayList;
import java.util.List;
import org.udel.titania.sms.core.SiteType;

/**
 * Question 6. d110 is only the repeat. A paper's "layer k" is read as the k-th O-(Ti2O2)-O trilayer
 * (bridging, in-plane, sub-bridging); that reading is an assumption and every row carries it.
 * A paper's "atomic layer n" is read as the n-th distinct oxygen height from the surface.
 */
public final class AtomicLayerMap {
    public static final String TRILAYER = "ASSUMED_TRILAYER_CONVENTION";
    public static final String ORDERED_O = "ASSUMED_ORDERED_O_HEIGHT";
    private final RutileStructure s;
    public AtomicLayerMap(RutileStructure s) { this.s = s; }

    public AtomicLayer trilayer(int k) {
        return new AtomicLayer("Layer " + k, k, s.depthNm(k, SiteType.BRI), s.depthNm(k, SiteType.IPL),
            s.depthNm(k, SiteType.SBR), TRILAYER);
    }

    /** One entry per oxygen height, ordered by depth, down to maxDepthNm. */
    public record OHeight(int atomicOLayer, int trilayer, SiteType site, double depthNm) {}
    public List<OHeight> oxygenHeights(double maxDepthNm) {
        List<OHeight> out = new ArrayList<>();
        for (int k = 1; ; k++) {
            if (s.depthNm(k, SiteType.BRI) > maxDepthNm) break;
            for (SiteType t : SiteType.values()) {
                double d = s.depthNm(k, t);
                if (d <= maxDepthNm) out.add(new OHeight(0, k, t, d));
            }
        }
        out.sort((x, y) -> Double.compare(x.depthNm(), y.depthNm()));
        List<OHeight> numbered = new ArrayList<>();
        for (int i = 0; i < out.size(); i++) { OHeight h = out.get(i); numbered.add(new OHeight(i + 1, h.trilayer(), h.site(), h.depthNm())); }
        return numbered;
    }
}
