package org.udel.titania.sms.crystallography;

import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import org.udel.titania.sms.core.SiteType;

/** Question 3. Only rutile (110) has layer-resolved energies, so it is the one facet, and every output says so. */
public final class FacetModel {
    public record Facet(String millerIndex, double areaFraction, Map<SiteType, Double> siteDensityPerNm2) {}
    public final List<Facet> facets;
    public final String label;
    private FacetModel(List<Facet> f, String label) { this.facets = f; this.label = label; }

    public static FacetModel rutile110(RutileStructure s) {
        Map<SiteType, Double> d = new EnumMap<>(SiteType.class);
        for (SiteType t : SiteType.values()) d.put(t, s.multiplicity(t) / s.cellArea110Nm2());
        return new FacetModel(List.of(new Facet("110", 1.0, d)), "rutile-110-equivalent");
    }
    public boolean facetWeighted() { return facets.size() > 1; }
}
