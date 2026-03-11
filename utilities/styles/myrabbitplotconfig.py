def eta_label(lo, hi):
    return f"Kaon eta [{lo:.1f}, {hi:.1f}]"


def charge_label(c):
    return f"Kaon charge {int(c)}"


translate_selection = {
    "bkmm_jpsimc_kaon1eta": eta_label,  # takes (low_edge, high_edge)
    "bkmm_kaon_charge": charge_label,  # takes bin center (e.g. -1 or 1)
}
