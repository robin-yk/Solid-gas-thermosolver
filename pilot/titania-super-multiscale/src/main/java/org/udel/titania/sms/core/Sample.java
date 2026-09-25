package org.udel.titania.sms.core;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/** Measured inputs of one sample. Null means not recorded; nothing is filled in. */
public record Sample(String id, Double inventory, Double inventoryAlt, String inventoryStatus, double rate,
                     String atmosphere, Double treatmentTC, Double treatmentTimeS, Double coolingRateKs,
                     double reactionTC, Double observationTimeS, String note) {

    public static List<Sample> load(Spec spec) {
        List<Sample> out = new ArrayList<>();
        for (Map<String, String> s : spec.samples) {
            Map<String, String> h = spec.histories.stream().filter(x -> x.get("sample").equals(s.get("sample")))
                .findFirst().orElseThrow(() -> new IllegalStateException("no history row for " + s.get("sample")));
            out.add(new Sample(s.get("sample"), num(s.get("inventory_umol_g")), num(s.get("inventory_alt_umol_g")),
                s.get("inventory_status"), Double.parseDouble(s.get("rate_co_umol_g_s")),
                h.get("atmosphere"), num(h.get("treatment_T_C")), num(h.get("treatment_time_s")),
                num(h.get("cooling_rate_K_s")), Double.parseDouble(h.get("reaction_T_C")),
                num(h.get("initial_rate_observation_time_s")), s.get("note")));
        }
        return out;
    }
    static Double num(String s) { return s == null || s.isBlank() ? null : Double.parseDouble(s); }

    public List<String> missingInputs() {
        List<String> m = new ArrayList<>();
        if (inventory == null) m.add("INVENTORY");
        if (atmosphere == null || atmosphere.isBlank()) m.add("ATMOSPHERE");
        if (treatmentTimeS == null) m.add("TREATMENT_TIME");
        if (coolingRateKs == null) m.add("COOLING_RATE");
        if (observationTimeS == null) m.add("INITIAL_RATE_OBSERVATION_TIME");
        m.add("INVENTORY_AND_RATE_ERROR_BARS");
        return m;
    }
}
