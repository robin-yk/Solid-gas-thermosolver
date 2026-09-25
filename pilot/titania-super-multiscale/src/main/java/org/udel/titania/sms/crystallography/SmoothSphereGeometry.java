package org.udel.titania.sms.crystallography;

/** Baseline: a smooth rutile sphere. Area and volume come from the same radius, so they agree. */
public final class SmoothSphereGeometry implements ParticleGeometry {
    private final double diameterNm, densityGcm3;
    public SmoothSphereGeometry(double diameterNm, double densityGcm3) {
        if (!(diameterNm > 0) || !(densityGcm3 > 0)) throw new IllegalArgumentException("diameter and density must be positive");
        this.diameterNm = diameterNm; this.densityGcm3 = densityGcm3;
    }
    public String id() { return "smooth-sphere-" + (int) Math.round(diameterNm) + "nm"; }
    public double diameterNm() { return diameterNm; }
    public double radiusNm() { return diameterNm / 2; }
    public double volumeNm3() { double r = radiusNm(); return 4.0 / 3.0 * Math.PI * r * r * r; }
    public double areaNm2() { double r = radiusNm(); return 4 * Math.PI * r * r; }
    public double massG() { return densityGcm3 * volumeNm3() * 1e-21; }
    public double areaM2PerGram() { return areaNm2() * 1e-18 / massG(); }
    public double shellVolumeNm3(double inner, double outer) {
        double r = radiusNm(), a = Math.max(0, r - Math.max(inner, outer)), b = Math.max(0, r - Math.min(inner, outer));
        return 4.0 / 3.0 * Math.PI * (b * b * b - a * a * a);
    }
    public double areaAtDepthNm2(double d) { double x = Math.max(0, radiusNm() - d); return 4 * Math.PI * x * x; }
}
