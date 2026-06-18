"""Persisted user preferences for the CHESS-QC desktop GUI (QSettings).

Stores the output-decimal count and the active theme so the user's choices
survive across sessions. Kept tiny and dependency-free on purpose.
"""
from __future__ import annotations

from .qt import QtCore
from .theme import THEMES, DEFAULT_THEME, VIBES, DEFAULT_VIBE, BADGE_STYLES, DEFAULT_BADGE

_ORG, _APP = "CHESS-QC", "QuickCompute"
DEFAULT_DECIMALS = 2
UNIT_CHOICES = ("SI", "US")
TCWIND_SPEED_CHOICES = ("km/h", "m/s")   # SI display unit for wind / TC-translation speeds
DEFAULT_TCWIND_SPEED = "km/h"


def _store() -> "QtCore.QSettings":
    return QtCore.QSettings(_ORG, _APP)


def get_decimals(default: int = DEFAULT_DECIMALS) -> int:
    try:
        n = int(_store().value("decimals", default))
    except (TypeError, ValueError):
        return default
    return min(max(n, 0), 8)


def set_decimals(n: int) -> None:
    _store().setValue("decimals", int(n))


def get_theme(default: str = DEFAULT_THEME) -> str:
    t = _store().value("theme", default)
    return t if t in THEMES else default


def set_theme(t: str) -> None:
    if t in THEMES:
        _store().setValue("theme", t)


def get_vibe(default: str = DEFAULT_VIBE) -> str:
    v = _store().value("vibe", default)
    return v if v in VIBES else default


def set_vibe(v: str) -> None:
    if v in VIBES:
        _store().setValue("vibe", v)


def get_badge(default: str = DEFAULT_BADGE) -> str:
    b = _store().value("badge", default)
    return b if b in BADGE_STYLES else default


def set_badge(b: str) -> None:
    if b in BADGE_STYLES:
        _store().setValue("badge", b)


def get_units(default: str = "SI") -> str:
    """Launcher-set unit system the calculators open in. Falls back to `default`
    (usually the app's own AppMeta.default_system) when unset/invalid."""
    u = _store().value("units", default)
    return u if u in UNIT_CHOICES else default


def set_units(u: str) -> None:
    if u in UNIT_CHOICES:
        _store().setValue("units", u)


def get_tcwind_speed(default: str = DEFAULT_TCWIND_SPEED) -> str:
    """Launcher-set SI display unit for wind / TC-translation speeds (km/h or m/s)."""
    s = _store().value("tcwind_speed", default)
    return s if s in TCWIND_SPEED_CHOICES else default


def set_tcwind_speed(s: str) -> None:
    if s in TCWIND_SPEED_CHOICES:
        _store().setValue("tcwind_speed", s)
