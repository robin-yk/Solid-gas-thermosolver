package org.udel.titania.sms.core;

/** Oxygen site classes of one rutile (110) O-(Ti2O2)-O trilayer. */
public enum SiteType {
    BRI("bridging, outermost O of the trilayer"),
    IPL("in-plane O; at the surface this is the basal site"),
    SBR("sub-bridging, innermost O of the trilayer");
    public final String meaning;
    SiteType(String m) { meaning = m; }
}
