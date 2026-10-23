# Loading IRFs and resolving names

Part of the [irfs architecture docs](README.md). Covers `irfLoader`,
`irfUtil` and the loader half of `latResponse::IrfLoader`.

## Three kinds of name

| Kind | Examples | Where it lives |
|------|----------|----------------|
| **Loader name** | `LATRESPONSE`, `DC1A`, `TEST` | `IrfRegistry` loader map. The argument of `irfLoader::Loader::go(name)`. |
| **Factory key** | `P8R3_SOURCE_V3::PSF2`, `P7REP_CLEAN_V15::FRONT`, `DC1A::Front`, `testIrfs::Back` | `IrfsFactory`. The argument of `IrfsFactory::create`. One key = one `Irfs` = one IRF and one event type. |
| **Response id** (`respIds` key) | `P8R3_SOURCE_V3`, `P8R3_SOURCE_V3 (PSF)`, `P8R3_SOURCE_V3 (EDISP)`, `P8R3_SOURCE_V3::PSF2`, `DC1A`, `DC1AF`, `TEST`, `TESTB` | `IrfRegistry::respIds()`. Maps a user-facing name to a list of factory keys. This is what `gtirfs` prints and what `Likelihood::ResponseFunctions::load` looks up. |

For a LAT IRF named `N`, `latResponse::IrfLoader::registerEventClasses`
creates:

| Response id | Factory keys |
|-------------|--------------|
| `N` | `N::FRONT`, `N::BACK` (partition `none` gets no suffix) |
| `N (PSF)` | `N::PSF0` … `N::PSF3` (Pass 8 only) |
| `N (EDISP)` | `N::EDISP0` … `N::EDISP3` (Pass 8 only) |
| `N::<type>` | itself, for each allowed event type |

The analytic loaders register `DC1A` → {`DC1A::Front`, `DC1A::Back`},
`DC1AF`, `DC1AB`, and `TEST` → {`testIrfs::Front`, `testIrfs::Back`},
`TESTF`, `TESTB`. `dc1aResponse` also adds `DC1A2::Front/Back` to the
factory without a response id.

`respIds` keys are ordered alphabetically (it is a `std::map`), so
`gtirfs` prints them in that order.

## Load sequence

### 1. Library load (static initialization)

`irfLoader/src/Loader.cxx` defines
`Registrar * Loader::s_registrar(new Registrar())`. Its constructor
(`Registrar.cxx`) runs when `libirfLoader` is loaded:

1. Reads `getenv("CALDB")`.
2. Calls `st_facilities::Environment::instance()`, which calls
   `facilities::commonUtilities::setupEnvironment()` once per process.
   That sets `CALDB`, `CALDBCONFIG` and `CALDBALIAS` only if they are
   unset, to paths under `$FERMI_INST_DIR/refdata/fermi/caldb`, which is
   not where the CALDB is installed (see
   `src/facilities/doc/architecture/`).
3. If `CALDB` was set before step 2, re-sets it with overwrite. Because
   `setupEnvironment` never overwrites, this changes nothing; it is a
   leftover from SCons builds.
4. Registers `new dc1aResponse::IrfLoader()` (unless compiled with
   `ST_USER`, which nothing defines), `new latResponse::IrfLoader()` and
   `new testResponse::IrfLoader()`.

`dc1Response`, `dc2Response` and `g25Response` loaders exist but are
**not registered**.

`latResponse::IrfLoader`'s constructor does real I/O:
- `new irfUtil::HdCaldb("GLAST", "LAT")` checks with
  `st_facilities::Util::file_ok` that `$CALDBALIAS` and `$CALDBCONFIG`
  exist (throws otherwise).
- `read_caldb_indx()` throws `runtime_error` if `CALDB` or `CALDBCONFIG`
  is unset or `caldb.config` has no `GLAST LAT` line. It then reads
  `caldb.indx` and lists IRF names (see [caldb.md](caldb.md#how-the-index-is-used)).
- `readCustomIrfNames()` throws if `CUSTOM_IRF_DIR` is set without
  `CUSTOM_IRF_NAMES`.

None of this is inside a `try`. An exception during static
initialization calls `std::terminate`, so a bad CALDB setup aborts any
program that loads `libirfLoader` before `main` runs (from reading the
code; not checked at runtime).

**Static-initialization order.** `latResponse` and `irfUtil` are static
libraries linked into `libirfLoader`, so their globals are initialized in
an unspecified order relative to `s_registrar`. Today the registrar path
uses only constant-initialized statics (`s_instance` pointers,
`s_interpolate_edisp`). `irfUtil::IrfHdus::s_aeff_cnames`/`s_psf_cnames`/
`s_edisp_cnames` are dynamically initialized `std::vector`s. Code that
used them from a loader constructor could see them still empty.

### 2. `Loader::go()`

```
Loader::go()                       // all registered loaders
  for name in IrfRegistry::irfNames():   go(name)
Loader::go(name)
  if name not a loader name → throw invalid_argument("Invalid IRF named ...")
  if name == "EGRET" and CALIB_DIR unset → skip
  IrfRegistry::loadIrfs(name)
    loader->loadIrfs()
    loader->registerEventClasses()
```

For `LATRESPONSE`:
- `loadIrfs()`: for every name from `caldb.indx`, gets the event-type map
  from `EventTypeMapper` and adds `new latResponse::Irfs(name, type)` to
  the factory under `name::type`, unless that key already exists. The
  `Irfs` constructor only records the names and sets `irfID` to the
  event-type bit. **No IRF file is opened.**
- `registerEventClasses()`: fills `respIds` as in the table above.
- Errors for individual IRFs (`runtime_error` from `EventTypeMapper`) are
  caught and ignored in both steps.

The first `EventTypeMapper::instance()` call reads `irf_index.fits`. It
happens here, unless `dataSubselector` got there first.

`go()` is idempotent and cheap to repeat. `Likelihood::ResponseFunctions::load`
calls `Loader_go()` every time.

`Loader_go()` and `Loader_respIds()` are free-function wrappers kept for
old Windows DLL linkage problems. `Likelihood` and `pyLikelihood` use
them.

### 3. First use of a component

`latResponse::Irfs::aeff()` (and `psf()`, `edisp()`,
`efficiencyFactor()`):
1. Builds `irfUtil::IrfHdus::aeff(name, type)`. For each component name
   (`EFF_AREA`, `PHI_DEP`, `EFFICIENCY_PARS` for aeff; `RPSF`,
   `PSF_SCALING` for psf; `EDISP`, `EDISP_SCALING` for edisp) it calls
   `HdCaldb::getFiles`, which runs `HDgtcalf` over `caldb.indx`. Each call
   rescans the index; there is no cache.
2. Passes the `IrfHdus` to `latResponse::IrfLoader::aeff/psf/edisp/efficiency_factor`,
   which builds a single-epoch object or an `*EpochDep` wrapper (see
   [lat-response.md](lat-response.md)).
3. Stores it in the `Irfs` (via `setAeff`, …).

`efficiencyFactor()` builds a second `IrfHdus::aeff`, so it repeats the
three aeff index scans.

Because factory prototypes stay empty, every `IrfsFactory::create()`
result loads its own copies. Two `create("P8R3_SOURCE_V3::PSF3")` calls
read the files twice.

## `irfUtil`

| Class | Role |
|-------|------|
| `EventTypeMapper` | Singleton built from `irf_index.fits` (see [caldb.md](caldb.md#irf_indexfits)). `mapping(irf)` → `{event type → (bit, partition)}` limited to the event types the IRF allows; `bitPos(irf, type)`; `getPartitions(irf)`; `full_bitmasks_by_partition(irf)`. Unknown IRF names throw `runtime_error("IRFs named X not found in CALDB.")`. `deleteInstance()` exists for tests. |
| `IrfHdus` | (file, HDU) pairs per component, one per epoch. `numEpochs()` is the size of the **alphabetically first** component's list (`EFFICIENCY_PARS` for aeff, `PSF_SCALING` for psf, `EDISP` for edisp), so all components must have the same number of epochs. `convType()` is 1 only for `BACK` (used by pre-Pass 8 PSFs). |
| `HdCaldb` | Wraps `HDgtcalf`. `getFiles(files, extnums, detName, respName, irfName)` returns all epochs; `operator()(det, resp, expr, filter, date, time)` returns one file for a date. Allocates 2 × 100 × 1024-byte buffers per instance. |
| `Util` | `joinPaths("a b c")` → `a/b/c`; `getCaldbFile(...)` (legacy single-date lookup that prints a hint and rethrows on failure). |

`irfUtil/irfUtil/RootTable.h` and `src/RootTable.cxx` (ROOT) and
`irfUtil/dgaus8.h` are installed or present but not compiled and not used.

## Resolving `irfs=CALDB` in consumers

irfs does not choose the IRF for a data file. The consumers do, and they
lean on `EventTypeMapper` and the `respIds` naming above.

**`dataSubselector::Cuts::CALDB_implied_irfs()`** (used by Likelihood):
1. Reads `BITMASK_MAPPING` from `irf_index.fits` itself
   (IRF name → `1 << BITPOSITION`).
2. With `irfs=CALDB`, chooses, among the IRFs whose event-class bit
   equals the file's `EVENT_CLASS` DSS mask and whose pass version matches
   the file (`P8R3`, …), the one with the highest `_V<n>`. With an
   explicit name it validates it and warns if a newer version exists.
3. `append_event_type_partition` takes `floor(log2(EVENT_TYPE mask))`
   and, if that bit's partition is not `none`, appends `" (PSF)"` or
   `" (EDISP)"`. For a mask with bits in several partitions, only the
   highest bit decides.

**`Likelihood::AppHelpers`** forbids `::` qualifiers in the `irfs`
parameter, adds an `EVENT_TYPE` cut from the `evtype` parameter if the
file has none, calls `CALDB_implied_irfs()`, and passes the result
(for example `P8R3_SOURCE_V3 (PSF)`) to `ResponseFunctions::load`, which
creates every factory key in `respIds[name]` and keeps those whose
`irfID()` is among the selected event-type bits.

**`gtobssim`** builds the response id itself as `irfs + " (" + evtype + ")"`
(its `evtype` parameter is a partition name such as `PSF`), looks it up in
`Loader::respIds()` directly, and uses
`EventTypeMapper::full_bitmasks_by_partition` to write the matching
`EVENT_TYPE` DSS cut.

## Gotchas

- `Loader::go(name)` expects a loader name. Passing an IRF name such as
  `P8R3_SOURCE_V3` throws `invalid_argument("Invalid IRF named ...")`,
  which reads as though the IRF did not exist.
- Silent drops: IRFs listed in `caldb.indx` but missing from
  `irf_index.fits` (all Pass 6 IRFs, `P7SOURCE_V6MC`) never reach the
  factory and produce no message. The only symptom is "Invalid response
  function choice" from Likelihood, or their absence from `gtirfs`.
- `CAL_DATE` ≤ `2007-01-01` rows are ignored entirely; the cutoff is
  hard-coded in `read_caldb_indx`.
- `irf_index.fits` is always read from `$CALDB/data/glast/lat/bcf/`, but
  `caldb.indx` comes from `CALDBCONFIG`. Pointing `CALDBCONFIG` at a
  different CALDB without changing `CALDB` mixes the two.
- `CUSTOM_IRF_DIR`/`CUSTOM_IRF_NAMES` are parsed but do nothing
  (`loadCustomIrfs` is commented out and never called). Old user docs
  that describe custom IRFs through these variables are stale.
- `latResponse::IrfLoader::loadIrfs(irfName)` builds its error message as
  `"IRF " + irfName + "not found in CALDB"` (missing space). It is always
  swallowed anyway.
- `IrfRegistry` and `IrfsFactory` are never cleared by `Loader`. Tests
  that call `IrfsFactory::delete_instance()` must call `Loader::go()`
  again.
- `EventTypeMapper` and `IrfsFactory` singletons are never deleted at
  exit (by design; avoids destruction-order problems).
