package org.udel.titania.sms.core;

/** Physical constants, read once from the registry so they carry a source. */
public final class Phys {
    public final double kB, NA, h;
    public Phys(Spec s) { kB = s.d("boltzmann_eV_K"); NA = s.d("avogadro"); h = s.d("planck_eV_s"); }
}
