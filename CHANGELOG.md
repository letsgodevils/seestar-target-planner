# Changelog

## 1.2.0 (2026-11-10)

- New **Min altitude** setting (10°–50° in 5° steps, default 25°): the lowest altitude counted as imaging time.
  Seasons, peak nights, dark windows, chart threshold lines and hours in darkness all follow it. The default can be
  changed in the `ss-config` block in `index.html`.
- **Moon distance rule:** a bright Moon (25% lit or more, above the horizon) no longer rules out a night if it stays
  at least 60° from a dual-band target or 90° from a UV/IR target. Dark windows count those nights; each chart adds
  a line saying how close the Moon comes tonight; new **Moon distance** sort lists the safest targets first.
- **Red light** button: switches the page to red-only light (no green or blue anywhere, charts included) with
  outlined buttons instead of solid fills, to protect night vision at the telescope. Remembered per browser.
- **Back to top** button appears after scrolling down.
- On phones, the filter/sort bar now scrolls away instead of staying pinned.
- Low-target thresholds never exceed the chosen setting (previously they could be 1–2° above 25° for targets whose
  highest point is 33–35°).

## 1.1.0 (2026-11-10)

- Tonight's altitude chart for every target and lookup: altitude curve, Moon path, twilight shading, the imaging
  threshold and a "now" marker, with hover/drag read-out of altitude, compass direction and time. Times use the
  device's time zone.
- Max altitude now counts only full darkness. The column shows tonight's highest point in darkness, with the best
  for the next 12 months underneath.
- New sort option: Max altitude (tonight, in darkness).
- Fixed: the no-JavaScript message showed a placeholder instead of text.
- Fixed: group sizes such as "9x4 (each)" no longer get a stray arcminute mark.

## 1.0.0 (2026-11-08)

First public release of the stand-alone web version.

- Top 50 galaxy and top 50 nebula targets, chosen for Seestar smart telescopes.
- For any location: season, peak night, next two moon-dark windows, maximum altitude, and low / near-zenith warnings.
- Telescope selector for the Seestar S30, S30 Pro, S50 and S50 Pro, with framing tags (fits, fits lengthwise,
  small, or mosaic size) and integration times rescaled by focal ratio.
- Lookup of any Messier, NGC or IC object (OpenNGC catalog), with estimated integration time.
- Town search via OpenStreetMap Nominatim, plus a built-in city list and manual coordinates.
- Filter by month, sort by rank, next dark window, peak date, integration or framing; pages of 10–50 or all.
- Wikipedia links for every target.
- Runs entirely in the browser; strict Content Security Policy with no inline scripts and no third-party
  scripts or fonts.
