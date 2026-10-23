# The bundled LAT CALDB

Part of the [irfs architecture docs](README.md).

The Fermitools do not use a separately installed HEASoft CALDB. The LAT
calibration database lives in this repository under `caldb/CALDB/` (about
48 MB, about 590 files). `caldb/CMakeLists.txt` installs it unchanged to
`${FERMI_INSTALL_DATADIR}/caldb`, which is `share/fermitools/data/caldb`.
Conda activation points `CALDB` at that directory (`$FERMI_DIR/data/caldb`).

Lookups go through `irfUtil`'s private copy of the HEASoft `HDgtcalf`
routine (`irfUtil/src/HDgtcalf*.c`, written by Ziqin Pan at GSFC in 2003,
built "standalone" without the rest of HEASoft).

## Layout

```
$CALDB/
  software/tools/
    caldb.config          # mission/instrument → index location
    alias_config.fits     # only checked for existence
  data/glast/lat/
    caldb.indx            # CIF: one row per (IRF, event type, component) HDU
    crcif.par             # HEASoft parameter file for building a CIF; unused at runtime
    bcf/
      irf_index.fits      # event-class and event-type bit mappings
      valid_evclass_selections.txt
      valid_evtype_selections.txt
      v2update.out        # log left from a maintenance script
      ea/     aeff_<IRF>_<FB|PSF|EDISP>.fits     (and older aeff_<irf>_<front|back>.fits)
      psf/    psf_<IRF>_<FB|PSF|EDISP>.fits
      edisp/  edisp_<IRF>_<FB|PSF|EDISP>.fits
```

`ea/` also contains `README` (about the DC1 file), `aeff_all.root` and
other legacy files. `psf/` and `edisp/` contain ROOT parameter files
(`*_thin_parameters.root`, `*_thick_parameters.root`). Nothing in the built
code reads the ROOT files.

`caldb.config` has a `GBM` line as well, but no GBM data is shipped.

## `caldb.config`

Whitespace-separated columns: mission, instrument, device, index
directory, index file, data device, data directory. The `GLAST LAT` line
is

```
GLAST LAT CALDB data/glast/lat caldb.indx CALDB data/glast/lat
```

`latResponse::IrfLoader::find_cif` takes columns 4 and 5 (`data/glast/lat`
+ `caldb.indx`) and joins them to `$CALDB`. `HDgtcalf` parses the same
file on its own.

## `caldb.indx`

A standard Calibration Index File (`CIF` extension, about 3,270 rows).
The columns that matter:

| Column | Use |
|--------|-----|
| `TELESCOP`, `INSTRUME` | Always `GLAST`, `LAT`. |
| `DETNAM` | The **event type**: `FRONT`, `BACK`, `PSF0`–`PSF3`, `EDISP0`–`EDISP3`. |
| `FILTER` | `NONE`/`NULL`; not used for selection. |
| `CAL_DIR`, `CAL_FILE`, `CAL_XNO` | File and **HDU number** of the component. |
| `CAL_CNAM` | Component: `EFF_AREA`, `PHI_DEP`, `EFFICIENCY_PARS`, `RPSF`, `PSF_SCALING`, `FISHEYE_CORR`, `EDISP`, `EDISP_SCALING`. |
| `CAL_CBD` | Boundary strings. The first is `VERSION(<IRF name>)`, which is how rows are matched to an IRF name. Others hold `CLASS(...)`, energy, `CTHETA` and `PHI` ranges. |
| `CAL_VSD`, `CAL_VST` | Validity start date and time (all `2007-01-17 00:00:00` in the shipped CALDB). |
| `CAL_QUAL` | Must be 0 for `HDgtcalf` to return the row. |
| `CAL_DATE` | Date the row was added. `latResponse` uses it as a filter (below). |

For one IRF (for example `P8R3_SOURCE_V3`) there is one row per event
type per component. For Pass 8 that is 10 event types × 8 components
(`EDISP` and `EDISP_SCALING` for edisp; `RPSF`, `PSF_SCALING`,
`FISHEYE_CORR` for psf; `EFF_AREA`, `PHI_DEP`, `EFFICIENCY_PARS` for
aeff). The event types share three files per component, one for each
partition: `_FB`, `_PSF` and `_EDISP`. Each file holds the event types of
its partition as consecutive HDUs.

### How the index is used

1. **Listing IRF names** (`latResponse::IrfLoader::read_caldb_indx`, at
   library load). For every row with `CAL_DATE` after `2007-01-01`, it
   tokenizes the first `CAL_CBD` string on `(` and `)` and takes the second
   token as the IRF name. The result is the list of names this loader will
   try. It includes Pass 6 names (`P6_V3_DIFFUSE`, …) and
   `P7SOURCE_V6MC`, which later fail silently (see
   [loading-and-names.md](loading-and-names.md)).
2. **Finding files** (`irfUtil::HdCaldb::getFiles`, on first use of a
   component). Calls `HDgtcalf("GLAST", "LAT", detnam=<event type>,
   "NONE", codename=<CAL_CNAM>, "-", "-", "-", "-",
   "VERSION.eq.<IRF>")`. With all dates set to `"-"` the validity-time
   filter is off, so **every epoch's row is returned, in index order**.
   The match is case-insensitive on `TELESCOP`, `INSTRUME`, `DETNAM`,
   `FILTER` and `CAL_CNAM` (each also matches a row value of `NONE`), and
   requires `CAL_QUAL == 0` and a `CAL_CBD` match. At most 100 rows are
   returned (`HdCaldb::s_maxret`). `IrfHdus` keeps them as (full path,
   HDU number as a string) pairs, which `tip` accepts as file + extension.
3. **Single-date lookups** (`HdCaldb::operator()` and
   `irfUtil::Util::getCaldbFile`) pass a date and time and return the
   first match with the latest validity start. Only the legacy
   `dc1Response`/`dc2Response` and `test_irfUtil` use these.

## `irf_index.fits`

Read by `irfUtil::EventTypeMapper` (and directly by
`dataSubselector::Cuts`). The path is hard-coded as
`$CALDB/data/glast/lat/bcf/irf_index.fits`; it does not come from
`caldb.config`.

| Extension | Columns | Meaning |
|-----------|---------|---------|
| `BITMASK_MAPPING` | `EVENT_CLASS`, `BITPOSITION`, `EVENT_TYPES` | One row per IRF name (50 rows: P7, P7REP, P8R2, P8R3 V2 and V3). `BITPOSITION` is the `EVENT_CLASS` bit in FT1 (for example `P8R3_SOURCE_V3` → 7, i.e. `evclass=128`). `EVENT_TYPES` is a mask of allowed event-type bits: 3 (FRONT/BACK) for Pass 7, 1023 (all 10) for Pass 8. |
| `EVENT_TYPE_MAPPING` | `EVENT_TYPE`, `BITPOSITION`, `EVENT_TYPE_PARTITION` | `FRONT` 0 and `BACK` 1 (partition `none`), `PSF0`–`PSF3` 2–5 (`PSF`), `EDISP0`–`EDISP3` 6–9 (`EDISP`). |

Pass 7 and Pass 8 reuse event-class bit positions, so `BITPOSITION` alone
does not identify an IRF; `dataSubselector` combines it with the pass
version from the FT1 header.

`caldb/python/irf_index.txt` is the text source for this file;
`make_irf_index.py` writes the FITS version.

## `valid_ev*_selections.txt`

Comma-separated `selection_mask,allowed_mask` pairs, generated by
`build_caldb.py` from `caldb/python/valid_selections_p8r3.yaml`.
`dataSubselector::BitMaskCut` reads them to decide which further event-class
and event-type selections are valid after an earlier one. The irfs
libraries do not read them.

## IRF FITS files

Each component HDU has `ENERG_LO`, `ENERG_HI` (MeV), `CTHETA_LO`,
`CTHETA_HI` as its first four columns, followed by the parameter columns.
Each cell is a vector over energy × cos θ. Headers carry the standard
CALDB keywords (`CCNM0001`, `CBD10001`, `CVSD0001`, `CVST0001`) and, for
PSF and energy dispersion, a version keyword:

| HDU (per event type) | Columns after the bounds | Version keyword |
|----------------------|--------------------------|-----------------|
| `EFFECTIVE AREA_<type>` | `EFFAREA` (m²) | — |
| `PHI_DEPENDENCE_<type>` | `PHIDEP0`, `PHIDEP1` | — |
| `EFFICIENCY_PARAMS_<type>` | `EFFICIENCY_PARS` (two rows: p0 and p1 coefficient sets; no bounds columns) | — |
| `RPSF_<type>` | `NCORE`, `NTAIL`, `SCORE`, `STAIL`, `GCORE`, `GTAIL` | `PSFVER` (3 for Pass 8 and P7REP) |
| `PSF_SCALING_PARAMS_<type>` | `PSFSCALE` (`c0`, `c1`, `beta`) | — |
| `FISHEYE_CORRECTION_<type>` | `MEAN`, `MEDIAN`, `PEAK` | — (never read) |
| `ENERGY DISPERSION_<type>` | `F`, `S1`, `K1`, `BIAS`, `BIAS2`, `S2`, `K2`, `PINDEX1`, `PINDEX2` | `EDISPVER` (3 for Pass 8, 2 for Pass 7) |
| `EDISP_SCALING_PARAMS_<type>` | `EDISPSCALE` (6 coefficients) | — |

How these are evaluated is in [lat-response.md](lat-response.md).

## Maintenance scripts

`caldb/python/` holds the scripts used to build and update the CALDB.
They are **not installed** and need `astropy`, `yaml` and, for
`build_caldb.py`, the old `GtApp` wrapper and an `INST_DIR` that contains
`irfs/caldb/python` (an SCons-era layout). Expect to adapt them before use.

| Script | Purpose |
|--------|---------|
| `build_caldb.py` | Builds a CALDB tree from IRF files, writes the valid-selection tables. |
| `IrfConsolidator.py` | Merges per-event-type IRF files into the `_FB`/`_PSF`/`_EDISP` files. |
| `make_irf_index.py` | Writes `irf_index.fits` from `irf_index.txt`. |
| `update_caldb_indx.py`, `correct_cbd.py`, `fix_calib_keywords.py` | Edit `caldb.indx` rows and CALDB header keywords. |
| `update_checksums.py`, `add_to_cvs.py` | Housekeeping (CVS-era). |

## Gotchas

- **Adding an IRF needs three edits:** the FITS files, rows in
  `caldb.indx` (with `CAL_DATE` after 2007-01-01 and `CAL_QUAL` 0), and a
  row in `irf_index.fits` `BITMASK_MAPPING`. Missing the last one makes
  the IRF disappear silently: it is listed from `caldb.indx`, but
  `EventTypeMapper::mapping()` throws and the loader swallows the error.
- **Epoch order comes from the index.** `HdCaldb::getFiles` returns rows
  in index order, and `EpochDep` assumes ascending start times (it uses
  `upper_bound`). Rows for a multi-epoch IRF must be in time order. No
  shipped IRF has more than one epoch.
- **Every component must be present** for every event type. `IrfHdus`
  queries every component name and `HdCaldb::getFiles` throws when a
  query finds nothing. A missing `PHI_DEP` or `EFFICIENCY_PARS` row makes
  the whole effective area fail to load, even though `Aeff` would cope
  with a missing phi table.
- `caldb/doc/release.notes` stops at 2017; the P8R3 `_V3` IRFs and the
  latest index fixes (git log: "Fix time boundary issue in caldb.indx")
  came later.
- The ROOT files, `v2update.out`, `crcif.par` and the DC1/Pass 4/Pass 5
  files are shipped and installed, but nothing can load them.
