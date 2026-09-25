package org.udel.titania.sms.crystallography;

/**
 * A measured BET area may enter only for the sample it was measured on, and it defines an area, not a
 * particle: multiplying it into spherical shell volumes is refused (the geometry validator).
 */
public final class SampleMatchedBetMorphology implements ParticleGeometry {
    private final String sampleId, betSampleId;
    private final double betM2g;
    public SampleMatchedBetMorphology(String sampleId, String betSampleId, double betM2g) {
        this.sampleId = sampleId; this.betSampleId = betSampleId; this.betM2g = betM2g;
    }
    public boolean matched() { return sampleId.equals(betSampleId); }
    private void require() {
        if (!matched()) throw new IllegalStateException("BET value of " + betSampleId + " does not belong to " + sampleId);
    }
    public String id() { return "bet-" + betSampleId; }
    public double areaM2PerGram() { require(); return betM2g; }
    public double radiusNm() { throw refuse(); }
    public double volumeNm3() { throw refuse(); }
    public double areaNm2() { throw refuse(); }
    public double massG() { throw refuse(); }
    public double shellVolumeNm3(double i, double o) { throw refuse(); }
    public double areaAtDepthNm2(double d) { throw refuse(); }
    private static IllegalStateException refuse() {
        return new IllegalStateException("a BET area does not define a particle; it cannot be multiplied into spherical shells");
    }
}
