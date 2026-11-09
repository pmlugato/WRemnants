def get_labels_colors_procs_sorted(procs):
    preferred = ["flatBkg", "bkgConst", "bkgExp", "background", "signal"]
    ordered = [p for p in preferred if p in procs]
    ordered.extend([p for p in procs if p not in ordered])

    labels_map = {
        "flatBkg": "Flat bkg",
        "bkgConst": "Flat bkg",
        "bkgExp": "Exp bkg",
        "background": "Background",
        "signal": "Signal",
    }
    colors_map = {
        "flatBkg": "lightgray",
        "bkgConst": "lightgray",
        "bkgExp": "tab:orange",
        "background": "tab:orange",
        "signal": "tab:blue",
    }
    fallback = ["tab:green", "tab:red", "tab:purple"]

    labels = [labels_map.get(p, p) for p in ordered]
    colors = []
    for i, proc in enumerate(ordered):
        colors.append(colors_map.get(proc, fallback[i % len(fallback)]))
    return labels, colors, ordered


def eta_label(lo, hi):
    return f"eta [{lo:.2f}, {hi:.2f}]"


def charge_label(c):
    return f"charge {int(c)}"


def pt_label(center):
    return f"pT bin {int(center)}"


translate_selection = {
    "bkmm_jpsimc_kaon1eta": eta_label,  # takes (low_edge, high_edge)
    "bkmm_kaon_eta": eta_label,  # takes (low_edge, high_edge)
    "bkmm_kaon_pt": pt_label,  # takes (low_edge, high_edge)
    "bkmm_kaon_charge": charge_label,  # takes bin center (e.g. -1 or 1)
}
