package org.udel.titania.sms.crystallography;

/** A literature layer label next to real coordinates. depthStatus says how the label was placed. */
public record AtomicLayer(String paperLayerLabel, int orderedIndex, double minDepthNm,
                          double representativeDepthNm, double maxDepthNm, String depthStatus) {}
