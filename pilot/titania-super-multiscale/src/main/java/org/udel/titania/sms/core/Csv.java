package org.udel.titania.sms.core;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/** Minimal RFC-4180 reader: quoted fields, commas inside quotes, a header row. */
public final class Csv {
    private Csv() {}

    public static List<Map<String, String>> read(Path path) throws IOException {
        List<String> lines = Files.readAllLines(path);
        List<Map<String, String>> rows = new ArrayList<>();
        if (lines.isEmpty()) return rows;
        List<String> header = split(lines.get(0));
        for (int i = 1; i < lines.size(); i++) {
            if (lines.get(i).isBlank()) continue;
            List<String> f = split(lines.get(i));
            Map<String, String> row = new LinkedHashMap<>();
            for (int j = 0; j < header.size(); j++) row.put(header.get(j), j < f.size() ? f.get(j) : "");
            rows.add(row);
        }
        return rows;
    }

    static List<String> split(String line) {
        List<String> out = new ArrayList<>();
        StringBuilder cur = new StringBuilder();
        boolean q = false;
        for (int i = 0; i < line.length(); i++) {
            char c = line.charAt(i);
            if (q) {
                if (c == '"' && i + 1 < line.length() && line.charAt(i + 1) == '"') { cur.append('"'); i++; }
                else if (c == '"') q = false;
                else cur.append(c);
            } else if (c == '"') q = true;
            else if (c == ',') { out.add(cur.toString()); cur.setLength(0); }
            else cur.append(c);
        }
        out.add(cur.toString());
        return out;
    }

    public static String esc(Object o) {
        String s = o == null ? "" : String.valueOf(o);
        if (s.contains(",") || s.contains("\"") || s.contains("\n")) return "\"" + s.replace("\"", "\"\"") + "\"";
        return s;
    }

    public static void write(Path path, List<String> header, List<List<Object>> rows) throws IOException {
        StringBuilder sb = new StringBuilder();
        sb.append(String.join(",", header)).append('\n');
        for (List<Object> r : rows) {
            for (int i = 0; i < r.size(); i++) { if (i > 0) sb.append(','); sb.append(esc(r.get(i))); }
            sb.append('\n');
        }
        Files.createDirectories(path.getParent());
        Files.writeString(path, sb.toString());
    }
}
