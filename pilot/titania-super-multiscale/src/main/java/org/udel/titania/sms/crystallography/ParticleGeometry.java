package org.udel.titania.sms.crystallography;

/** Question 1-2: what a particle is, in the units every other module needs. */
public interface ParticleGeometry {
    String id();
    double radiusNm();
    double volumeNm3();
    double areaNm2();
    double massG();
    double areaM2PerGram();
    /** Volume between two depths measured inward from the surface. */
    double shellVolumeNm3(double innerDepthNm, double outerDepthNm);
    /** Area of the surface parallel to the outer one, at a given depth. */
    double areaAtDepthNm2(double depthNm);
}
