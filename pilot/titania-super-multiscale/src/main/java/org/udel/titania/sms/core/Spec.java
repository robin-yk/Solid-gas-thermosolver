package org.udel.titania.sms.core;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/** Every number the model uses comes from specification/ or experimental-data/, never from code. */
public final class Spec {
    public record Param(String id, String raw, String unit, String role, String source) {}

    private final Map<String, Param> params = new LinkedHashMap<>();
    public final Path root;
    public final List<Map<String, String>> energySets, transportEdges, polaronStates, reconstructionLimits, samples, histories;

    public Spec(Path root) throws IOException {
        this.root = root;
        for (Map<String, String> r : Csv.read(root.resolve("specification/parameter_registry.csv")))
            params.put(r.get("id"), new Param(r.get("id"), r.get("value"), r.get("unit"), r.get("role"), r.get("source")));
        energySets = Csv.read(root.resolve("specification/energy_sets.csv"));
        transportEdges = Csv.read(root.resolve("specification/transport_edges.csv"));
        polaronStates = Csv.read(root.resolve("specification/polaron_states.csv"));
        reconstructionLimits = Csv.read(root.resolve("specification/reconstruction_limits.csv"));
        samples = Csv.read(root.resolve("experimental-data/samples.csv"));
        histories = Csv.read(root.resolve("experimental-data/treatment_histories.csv"));
    }

    public Param param(String id) {
        Param p = params.get(id);
        if (p == null) throw new IllegalArgumentException("parameter not in registry: " + id);
        return p;
    }
    public double d(String id) { return Double.parseDouble(param(id).raw()); }
    public List<Double> list(String id) {
        List<Double> out = new ArrayList<>();
        for (String s : param(id).raw().split(";")) out.add(Double.parseDouble(s.trim()));
        return out;
    }
    public String source(String id) { return param(id).source(); }
    public Map<String, Param> all() { return params; }
}
