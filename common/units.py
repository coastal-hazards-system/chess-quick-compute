"""Unit conversion for the calculator edge.

Compute applications always work in SI. The GUI shows values in the selected unit
system (SI or US) and converts at the edge. Each unit symbol maps to a factor that
converts a value *in that unit* to SI:  value_SI = value_in_unit * TO_SI[unit].
"""

TO_SI = {
    "": 1.0,            # dimensionless
    # length
    "m": 1.0, "ft": 0.3048, "km": 1000.0, "mi": 1609.344, "nm": 1852.0,
    # time
    "s": 1.0, "hr": 3600.0,
    # velocity
    "m/s": 1.0, "km/h": 1000.0 / 3600.0, "ft/s": 0.3048, "mph": 0.44704, "kt": 0.514444,
    # temperature difference (no scale change between C and K)
    "C": 1.0,
    # atmospheric pressure (canonical Pa)
    "hPa": 100.0, "mb": 100.0, "inHg": 3386.389,
    # acceleration
    "m/s^2": 1.0, "ft/s^2": 0.3048,
    # pressure
    "Pa": 1.0, "psf": 47.880259,
    # energy density (per unit area -> per length here): N/m vs lb/ft
    "N/m": 1.0, "lb/ft": 14.593903,
    # energy flux: N/s vs lb/s
    "N/s": 1.0, "lb/s": 4.4482216,
    # angle
    "deg": 1.0, "rad": 1.0,
    # percent
    "%": 1.0,
    # unit weight (specific weight): canonical N/m^3
    "N/m^3": 1.0, "kN/m^3": 1000.0, "lb/ft^3": 157.08746,
    # weight / force: canonical N
    "N": 1.0, "kN": 1000.0, "lb": 4.4482216, "tons": 8896.4432,   # tons = US short ton-force
    # force per unit length (N/m vs lb/ft already above); moment per unit length
    "N-m/m": 1.0, "lb-ft/ft": 4.4482216,   # (47.880259 * 0.3048^2 = 4.4482)
    # volume transport rate: canonical m^3/yr
    "m^3/yr": 1.0, "yd^3/yr": 0.764554858,   # 1 yd^3 = 0.9144^3 m^3
    # sea-level trend rate: canonical m/yr (year identical in both systems)
    "yr": 1.0, "m/yr": 1.0, "ft/yr": 0.3048, "mm/yr": 0.001, "in/yr": 0.0254,
    # volume: canonical m^3
    "m^3": 1.0, "yd^3": 0.764554858,
    # overtopping rate per unit length (= area flux m^2/s): canonical m^3/s/m
    "m^3/s/m": 1.0, "ft^3/s/ft": 0.09290304,   # 1 ft^3/s/ft = 0.3048^2 m^3/s/m
    # phi grain-size scale (dimensionless; identical in both systems)
    "phi": 1.0,
    # areas, flows, impulse per area, inverse length, nautical mile
    "m^2": 1.0, "ft^2": 0.09290304,
    "m^3/s": 1.0, "ft^3/s": 0.028316846592,
    "m^2/s": 1.0, "ft^2/s": 0.09290304,
    "m^2/s^2": 1.0, "ft^2/s^2": 0.09290304,
    "N*s/m^2": 1.0, "lb*s/ft^2": 47.880259,
    "1/m": 1.0, "1/ft": 3.280839895,
    "nmi": 1852.0,
    # unit strings carried through unchanged (identical in SI and US, or already
    # expressed in their own unit): declared so both front-ends agree explicitly
    # rather than relying on the unknown-unit fallback
    "kg/m^3": 1.0, "deg C": 1.0, "1/s": 1.0, "1/yr": 1.0, "h": 1.0, "min": 1.0,
    "kg": 1.0, "tonne": 1.0, "x1000": 1.0, "m^(1/3)": 1.0, "m^B": 1.0,
    "s or 1/s": 1.0, "g": 1.0, "mm": 1.0, "log10 yr": 1.0,
}


def factor(unit: str) -> float:
    return TO_SI.get(unit, 1.0)


def to_si(value: float, unit: str) -> float:
    """Value given in `unit` -> SI."""
    return value * factor(unit)


def from_si(value_si: float, unit: str) -> float:
    """SI value -> displayed in `unit`."""
    f = factor(unit)
    return value_si / f if f else value_si


def convert(value: float, from_unit: str, to_unit: str) -> float:
    """Convert a value expressed in `from_unit` to `to_unit` (same physical dim)."""
    return from_si(to_si(value, from_unit), to_unit)
