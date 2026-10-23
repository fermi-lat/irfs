# irfInterface: the abstract IRF interface

Part of the [irfs architecture docs](README.md).

`irfInterface` (shared library, headers in `irfInterface/irfInterface/`)
is what nearly every client includes. It defines the component interfaces,
the `Irfs` container, two global registries, and generic numerical
fallbacks (integration, sampling, containment) that the concrete classes
may override.

## Classes

### `Irfs`

`irfInterface/Irfs.h`, header-only.

- Holds raw pointers to `IAeff`, `IPsf`, `IEdisp` and an optional
  `IEfficiencyFactor`, plus `int irfID`.
- **`irfID` is the event-type bit position** for Pass 7/8 IRFs (FRONT=0,
  BACK=1, PSF0–3=2–5, EDISP0–3=6–9). It is set from
  `EventTypeMapper::bitPos` in `latResponse::Irfs`. `Likelihood` uses it as
  the key of its `respPtr` map and as the event type of each event. The
  analytic `DC1A`/`TEST` responses use 0 (front) and 1 (back).
- The copy constructor **deep-copies** every component with `clone()`.
  `clone()` is the copy constructor.
- The destructor deletes all four components.
- `aeff()`, `psf()`, `edisp()` and `efficiencyFactor()` are `virtual`, so
  `latResponse::Irfs` can load components lazily.
- `setAeff`/`setPsf`/`setEdisp` store the new pointer and take ownership,
  but **do not delete the old one**. `setEfficiencyFactor` stores a clone
  of its argument (the caller keeps ownership of the argument) and also
  leaks any previous factor.
- `efficiencyFactor()` is `const` but returns a `const IEfficiencyFactor*`.
  The LAT implementation uses `const_cast` to load it lazily.

### `IAeff`

`IAeff.h`, header-only.

| Method | Notes |
|--------|-------|
| `value(E, srcDir, scZ, scX, t)` | Effective area in cm² for a sky direction and spacecraft attitude. |
| `value(E, theta, phi, t)` | Same in instrument coordinates (degrees). |
| `value(vector E, vector theta, vector phi, t)` | Loop over the scalar version; throws if sizes differ. |
| `upperLimit()` | Maximum over the table, used by `observationSim` for rejection sampling. |
| `setPhiDependence(bool)` / `usePhiDependence()` | Toggles azimuthal modulation. `Likelihood` turns it off for livetime-cube exposure (`ExposureCube`, `gtexpcube2`, `Drm`) and sets it per source in `PointSource`. |

`IAeff` has no constructor, so `m_usePhiDependence` starts uninitialized
in any subclass that does not set it. `latResponse::Aeff` sets it; the
analytic `dc1a`/`test` `Aeff` classes and `latResponse::AeffEpochDep` do
not (from reading the code; harmless unless a caller reads the flag
without setting it).

### `IPsf`

`IPsf.h`, `src/IPsf.cxx`.

| Method | Default implementation |
|--------|------------------------|
| `value(appDir, E, srcDir, scZ, scX, t)`, `value(sep, E, theta, phi, t)` | Pure virtual. `sep` and angles in degrees; the result is a density per steradian. |
| `angularIntegral(E, theta, phi, radius, t)` | Gaussian quadrature (`st_facilities::GaussianQuadrature::integrate`, f2c `dgaus8`) from 0 to `radius` degrees, relative error 1e-5. |
| `angularIntegral(E, srcDir, scZ, scX, cones, t)` and `(E, srcDir, theta, phi, cones, t)` | Integral of the PSF over the ROI cone, for a source at angle `psi` from the cone centre. Splits into a full-circle part (`mu` from `cos(r-psi)` to 1, skipped when that is ≥ 0.99) and a partial-arc part. **Uses only `cones.front()`**; any further cones are ignored. |
| `angularContainment(E, theta, phi, frac, t, rtol)` | Root-finds the radius whose integral is `frac` × the integral to 180°. |
| `appDir(E, srcDir, scZ, scX, t)` | Draws a reconstructed direction: tabulates the radial integrand on a fixed grid of 101 offsets (0 and 100 log-spaced values from 1e-4° to 90°), builds a piecewise-linear CDF, inverts it, and adds a uniform azimuth. Offsets beyond 90° can never be drawn. `USE_OLD_IPSF_SAMPLER` selects the older trapezoid/linear inversion. |

Every scalar `scZ`/`scX`-based method in the base class computes `theta`
from `scZ` and **sets `phi = 0`**; `scX` is ignored. The concrete LAT PSF
and energy dispersion ignore `phi` anyway (see
[lat-response.md](lat-response.md)).

**Static state.** All the integrals above pass `energy`, `theta`, `phi`,
`time` and the `this` pointer to the integrand through **static** members
(`s_energy`, `s_self`, `s_cp`, …), because the f2c integrator takes a
plain function pointer. Consequences:
- Not reentrant and not thread safe.
- An integrand that itself calls another `IPsf` integral would overwrite
  the statics. `Psf3::normalize_pars` calls `IPsf::angularIntegral` for
  energies below 120 MeV during construction. That is safe only because
  the call does not nest.
- `s_psi_values` (the `appDir` grid) is filled by the first `IPsf`
  constructor.

### `IEdisp`

`IEdisp.h`, `src/IEdisp.cxx`.

| Method | Default implementation |
|--------|------------------------|
| `value(appE, E, srcDir, scZ, scX, t)`, `value(appE, E, theta, phi, t)` | Pure virtual. Density per MeV of measured energy. |
| `integral(emin, emax, E, …)` | `dgaus8` over `[emin, emax]` with relative error 1e-7. If that fails with error code 2, it retries with the tolerance relaxed by 1000×. If the result is still 0 or the integrator throws, it sums `dgaus8` over 4 log-spaced pieces of ratio `(emax/emin)^(1/3)` starting at `max(emin, 1e-5)`. Four pieces of that ratio end at `emax·(emax/emin)^(1/3)`, so this fallback **integrates past `emax`** (suspected off-by-one, from reading the code; not checked at runtime). The `theta/phi` overload returns 0 instead of throwing when `E` is more than 1000× outside `[emin, emax]`. |
| `appEnergy(E, srcDir, scZ, scX, t)` | Samples a measured energy from a piecewise-linear CDF on 200 log-spaced points over `[E/10, 10E]`. Measured energies outside that range are never drawn. |
| `meanAppEnergy(E, …)` | Ratio of `dgaus8` integrals over `[0, 10E]`; normalizes explicitly because "the energy dispersion is not guaranteed to be properly normalized". |
| `meanTrueEnergy(appE, …)` | Same over true energy `[30 MeV, min(176 GeV, 10·appE)]`. The hard-coded limits suit the LAT only. |

Unlike `IPsf`, `IEdisp` integrates through functor objects
(`EdispIntegrand`, …) and the templated `GaussianQuadrature::dgaus8`, so it
has no static state.

### `IEfficiencyFactor`

`IEfficiencyFactor.h`, `src/IEfficiencyFactor.cxx`.

The livetime-fraction efficiency correction. The effective area in the
CALDB is for a nominal livetime fraction; the actual efficiency depends on
energy and the instantaneous livetime fraction (dead time from the trigger
rate). The LAT implementation is linear:
`factor = p0(logE)·livetimefrac + p1(logE)`.

| Method | Notes |
|--------|-------|
| `value(E, ltfrac, front, met)` | Pure virtual. Pass 8 implementations ignore `front`. |
| `value(E, ltfrac, met)` | Average of the front and back values. Returns 1 if `OMIT_EFFICIENCY_FACTOR` is set. This is the call `Likelihood` makes (`PointSource`, `EventContainer`, exposure cubes through `setEfficiencyFactor`). |
| `getLivetimeFactors(E, f1, f2, met)` | Returns `p1` and `p0` so callers can fold the correction into livetime-weighted sums. No caller outside irfs. |
| `operator()(E, met)` | Looks up the livetime fraction at `met` from FT2 data loaded with `readFt2File`. |
| `readFt2File(ft2)`, `clearFt2Data()` | Fill or clear the `START`/`STOP`/`livetime/(STOP-START)` vectors. |

### `IrfsFactory`

`IrfsFactory.h`, `src/IrfsFactory.cxx`. Process-wide singleton
(`instance()`, `delete_instance()`).

- `addIrfs(name, Irfs*)` takes ownership of a **prototype**. If the name
  exists, it replaces it; the old prototype is deleted and a warning
  printed only when `verbose=true`. Otherwise the old one leaks.
- `create(name)` returns `prototype->clone()`, owned by the caller. An
  unknown name throws `std::invalid_argument` listing every valid name.
  With all the LAT IRFs loaded that list is several hundred lines.
- `irfNames()`/`getIrfsNames()` return names in insertion order.
- Factory keys are `"<IRF>::<event type>"`, for example
  `P8R3_SOURCE_V3::PSF2`, `DC1A::Front`, `testIrfs::Back`.

### `IrfRegistry` and `IrfLoader`

`IrfRegistry.h`, `IrfLoader.h`, `src/IrfRegistry.cxx`. Singleton.

- `registerLoader(IrfLoader*)` stores the loader under `loader->name()`.
  Ownership is never released (loaders live for the whole process).
- `loadIrfs(loaderName)` calls `loader->loadIrfs()` and then
  `loader->registerEventClasses()`. Unknown names throw `runtime_error`.
- `registerEventClasses(name, keys)` and `registerEventClass(name, key)`
  fill the **`respIds` map**: user-facing name → list of factory keys.
  This map is what `irfLoader::Loader::respIds()` returns and what
  `Likelihood::ResponseFunctions::load` looks names up in.
- `irfNames()` lists **loader** names, not IRF names.

`IrfLoader` is the plug-in interface: `name()`, `loadIrfs()`,
`registerEventClasses()`, all `const`.

### `AcceptanceCone`

`AcceptanceCone.h`, header-only. A centre `SkyDir` and a radius in
degrees.
- `operator==` uses tolerances of 1e-5° on the radius and 1e-3 **radians**
  (about 0.06°) on the centre. `Psf3` uses `!=` to decide whether its ROI
  integral cache is stale, so moving an ROI centre by less than about
  0.06° reuses the old cache.
- `inside()` and `clone()` are non-const.

## Gotchas

- `irfInterface` links `st_facilities` `PUBLIC` and `astro`, `tip`, CLHEP
  `PRIVATE`. Its headers forward-declare `astro::SkyDir`, except
  `AcceptanceCone.h`, which includes `astro/SkyDir.h`.
- `irfInterface/data/P6_V3_DIFFUSE_eff.txt` is installed to the root of
  `share/fermitools/data/`. Nothing reads it.
- The test suite (`src/test/`) defines toy `Aeff`/`Psf`/`Edisp` classes and
  a `MyIrfLoader`. These are good templates for a minimal implementation.
