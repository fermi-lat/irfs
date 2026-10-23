# latResponse: the LAT IRF implementation

Part of the [irfs architecture docs](README.md).

`latResponse` (static library) implements the irfInterface classes from
the CALDB FITS files. Public headers are in `latResponse/latResponse/`;
several of the classes it instantiates (`Irfs`, `Psf`, `Psf2`, `Edisp`,
`Edisp2`, `EfficiencyFactor`, the `*EpochDep` wrappers, `PsfIntegralCache`,
`CaldbDate`) have headers only in `latResponse/src/`.

## Component selection

`latResponse::IrfLoader` has static factory functions that turn an
`irfUtil::IrfHdus` into a component:

| Function | Single epoch | More than one epoch |
|----------|--------------|---------------------|
| `aeff(hdus)` | `Aeff(hdus, 0)` | `AeffEpochDep` holding one `Aeff` per epoch |
| `psf(hdus)` | by `PSFVER` in the `RPSF` header: 1 (or missing) → `Psf`, 2 → `Psf2`, 3 → `Psf3`; anything else throws | `PsfEpochDep` |
| `edisp(hdus)` | by `EDISPVER`: 1 (or missing) → `Edisp`, 2 → `Edisp2`, 3 → `Edisp3` | `EdispEpochDep` |
| `efficiency_factor(hdus)` | `EfficiencyFactor(hdus, 0)` | `EfficiencyFactorEpochDep` |

Each version lookup opens the FITS table once just to read the keyword.
Pass 8 and P7REP PSFs are `PSFVER=3`; Pass 8 energy dispersion is
`EDISPVER=3`, Pass 7 is 2. `Psf`, `Psf2`, `Edisp` and `Edisp2` exist for
older files. `Psf3` reuses `Psf2`'s static `psf_base_function` and
`psf_base_integral`.

The epoch start of each component comes from `EpochDep::epochStart`, which
reads `CVSD0001`/`CVST0001` from that HDU's header and converts them to
MET with `astro::JulianDate`.

## Tables and interpolation

### `FitsTable` and `Bilinear`

`FitsTable(file, extname, column, nrow)` reads one parameter column:
- `ENERG_LO/HI` become log10 bin centres (`log10(sqrt(lo*hi))`) and log10
  bin edges.
- `CTHETA_LO/HI` become bin centres and edges.
- `minCosTheta()` is the lowest `CTHETA_LO` (0.2 for Pass 8, about 78.5°).
- `maximum()` is the largest table value.

`value(logE, costh, interpolate)`:
- `interpolate=true`: **bilinear** in (log10 E, cos θ) between bin
  centres. `Bilinear` pads the grid with extra rows and columns at
  log10 E = 0 and 10 and cos θ = −1 and 1 that **repeat the edge values**
  ("replicate nasty THF2 and RootEval::Table behavior from
  handoff_response"). So values are flat outside the tabulated centres,
  and anything above 10^10 MeV throws `invalid_argument`. cos θ above the
  last bin centre is clamped to it.
- `interpolate=false`: the value of the bin containing the point,
  clamped to the edge bins.

`ParTables` holds one `FitsTable` per parameter column (every column after
the four bounds columns, which must come first and in order).

`st_facilities` and `skymaps` have their own copies of `FitsTable`,
`Bilinear` and `Aeff`/`EfficiencyFactor` code. They are not kept in sync.

### `getVectorData`

`FitsTable::getVectorData(table, column, values, nrow)` reads row `nrow`
of a vector column as `float`, then widens to `double`. All IRF
parameters therefore have single precision, whatever the file stores.

## Effective area: `Aeff`

`latResponse/Aeff.h`, `src/Aeff.cxx`.

```
value(E, theta, phi) =
    0                                    if cos θ < minCosTheta
    EFFAREA(log10 E, min(cos θ, 0.99999)) [bilinear] × 1e4      (m² → cm²)
  × phi_modulation(log10 E, cos θ, phi)  [nearest-bin parameters]
```

- **Phi modulation.** With parameters `p0 = PHIDEP0`, `p1 = PHIDEP1` taken
  from the bin (no interpolation):
  `norm·(1 + p0·x^p1)` with `x = |2·(2/π·|φ' | − 1/2)|`,
  `φ' = fmod(φ, 180°) − 90°` and `norm = 1/(1 + p0/(1 + p1))`. The
  modulation is periodic in φ with period 180°, symmetric about its
  midpoint, and averages to 1 over φ.
  Negative φ is shifted by 360°.
- `setPhiDependence(false)` turns the modulation off. The constructor sets
  it on if the `PHI_DEP` HDU could be read. (With `IrfHdus` the HDU must
  exist anyway, see [caldb.md](caldb.md#gotchas).)
- `value(E, srcDir, scZ, scX, t)` computes θ from `scZ` and φ from `scX`
  and `scZ × scX`. This is the only LAT component that uses φ.
- `time` is ignored; epochs are handled by `AeffEpochDep`.
- `upperLimit()` = `maximum() × 1e4 × max_phi_modulation()`.
- `pars()`, `logEnergies()`, `costhetas()`, `values()`, `setValues()` give
  the Python bindings access to the raw table.

## Livetime efficiency: `EfficiencyFactor`

`src/EfficiencyFactor.h/.cxx`.

- Reads the `EFFICIENCY_PARS` HDU: two rows, each a 6-element vector
  that becomes an `EfficiencyParameter`, a continuous three-segment
  piecewise-linear function of log10 E (slopes `a0`, `a1`, `a2`, intercept
  `b0`, break points `logEb1`, `logEb2`).
- If every coefficient is zero, `m_havePars=false` and every method
  returns 1 (or `factor1=1, factor2=0`).
- `value(E, ltfrac, front, met) = p0(log10 E)·ltfrac + p1(log10 E)`.
  `front` and `met` are ignored for Pass 8.
- `operator()(E, met)` needs `readFt2File()` first. Without FT2 data it
  returns 1. With FT2 data it throws if `met` is outside the FT2 range.
- Another constructor reads `EFFICIENCY_PAR_FILE` (text with two lines of
  coefficients, or a FITS file).

Likelihood uses `IEfficiencyFactor::value(E, ltfrac, met)`, the front/back
average defined in the base class. Because Pass 8 ignores `front`, it is
the same as the single value. `Likelihood::ResponseFunctions::efficiencyFactor()`
returns the factor of the **first** loaded event type and applies it to
all event types (a Likelihood design choice, noted here because it
depends on the per-type factors being similar).

## PSF: `Psf3`

`latResponse/Psf3.h`, `src/Psf3.cxx`, base `PsfBase`.

**Model.** For each (energy, cos θ) bin the six parameters `NCORE`,
`NTAIL`, `SCORE`, `STAIL`, `GCORE`, `GTAIL` define a sum of two King
functions in the scaled separation:

```
S(E)  = sqrt((c0·(E/100 MeV)^beta)² + c1²)          # PSFSCALE = (c0, c1, beta)
u_c   = (sep / (SCORE·S))² / 2,  u_t = (sep / (STAIL·S))² / 2
K(u,g)= (1 − 1/g)·(1 + u/g)^(−g)                    # g==1 is nudged to 1.001
P     = NCORE·K(u_c, GCORE) + NTAIL·NCORE·K(u_t, GTAIL)
```

Note that `NTAIL` is relative to `NCORE`.

**Normalization.** After reading, `normalize_pars()` divides `NCORE` in
every bin by the integral to 90°. Above 120 MeV it uses the analytic King
integral. Below 120 MeV it calls the numerical
`IPsf::angularIntegral` (static-variable integrator). `setParams()`
renormalizes all bins.

**Interpolation.** `value(sep, E, θ)` evaluates the PSF at the four
surrounding grid corners, **each at its own corner energy** (so the
energy scaling `S(E)` uses the corner energy), and bilinearly
interpolates the four values in (log10 E, cos θ). Parameters are never
interpolated. The grid is padded like `Bilinear` (log10 E 0–10,
cos θ −1–1, edge values repeated); energies above 10^10 MeV throw
`invalid_argument` from `findIndex` (after printing the range to
`std::cout`).

**φ and time are ignored**, and the `scZ`/`scX` overload uses only θ.

**Integrals.**
- `angularIntegral(E, θ, φ, radius)`: analytic King integrals at the four
  corners, bilinearly interpolated; below 120 MeV it falls back to the
  numerical base-class integral.
- `angularIntegral(E, srcDir, θ, φ, cones)` (the ROI integral used by
  Likelihood): uses a `PsfIntegralCache` for `cones.at(0)`. The cache is
  rebuilt whenever the cone differs (by `AcceptanceCone::operator!=`) from
  the cached one. It tabulates the ROI integral of a unit King function
  on a grid of 501 source offsets ψ (0 plus 500 log-spaced values) × σ ×
  γ, filling grid points lazily on demand. Values outside the σ/γ grid
  are integrated directly. Each corner integral is linearly interpolated
  in ψ, then the corners are combined bilinearly. A source farther from
  the ROI centre than the last ψ grid value throws `runtime_error`.
- `appDir` (sampling) is the base-class tabulated inversion.

## Energy dispersion: `Edisp3`

`latResponse/Edisp3.h`, `src/Edisp3.cxx`.

**Model.** With `x = (E' − E)/E` scaled by
`s(log10 E, |cos θ|) = a0·L² + a1·c² + a2·L + a3·c + a4·L·c + a5`
(`EDISPSCALE`, 6 coefficients):

```
D(E'|E) = [F·g(x/s; S1, K1, BIAS, PINDEX1) + (1−F)·g(x/s; S2, K2, BIAS2, PINDEX2)] / (E·s)
g(x; σ, k, b, p) = p/(σ·Γ(1/p)) · k/(1+k²) · exp(−(k·max(x−b,0)/σ)^p) · exp(−(max(b−x,0)/(k·σ))^p)
```

This is the two-sided generalized-Gaussian ("Thibaut") function
referenced in the source. `Γ` uses a single-precision-era `gammln`
approximation. A negative value throws `runtime_error`.

**Evaluation.**
- Default: an `EdispInterpolator` built lazily on the first `value()`
  call (it rereads the FITS table). Like `Psf3`, it evaluates the function
  at the four (energy, θ) grid corners and interpolates bilinearly. Before
  evaluating it maps the requested `E'` to each corner's own scaled
  variable so that the interpolation follows the dispersion's shape, and
  it rescales each corner value by the corner `E·s`.
- `DISABLE_EDISP_INTERP` set: parameters of the containing bin
  (`ParTables::getPars(..., interpolate=false)`), cached for the last
  (log10 E, cos θ) pair.
- cos θ is clamped to the last bin centre. **φ and time are ignored.**
- `Edisp3::renormalize` is a deliberate no-op: the Pass 8 parameters are
  taken to be normalized. `IEdisp::meanAppEnergy` still normalizes
  numerically.
- `IrfLoader::set_edisp_interpolation()`/`interpolate_edisp()` affect only
  `Edisp2`; the corresponding code in `Edisp3` is commented out.

## Epoch-dependent wrappers

`AeffEpochDep`, `PsfEpochDep`, `EdispEpochDep` and
`EfficiencyFactorEpochDep` derive from the interface and from `EpochDep`.
They hold one clone per epoch and dispatch on `EpochDep::index(met)`:
- `met` before the first epoch start throws "Requested MET not covered by
  selected IRFs".
- `met` at or after the last start uses the last epoch (no upper bound).
- `m_curr_index` caches the last epoch (mutable, not thread safe).

**No IRF in the shipped CALDB has more than one epoch**, so these classes
run only in `test_latResponse` (with `latResponse/data/*_epoch_*.fits`).
They have latent bugs (from reading the code; not checked at runtime):
- `EpochDep::index` uses `it - m_epochStart.end() - 1` where it means
  `begin()`. Going back to an earlier, non-last epoch after a later one
  yields a huge index. The unit test queries epoch 0 and then epoch 1 on
  separate objects, so it does not catch this.
- `EpochDep::epochStart` computes hours as `h + 60·(m + 60·s)` instead of
  `h + m/60 + s/3600`. Correct only for `00:00:00`.
- `EfficiencyFactorEpochDep::getLivetimeFactors` calls `value()` on the
  component, so its two output arguments are never set.
- `AeffEpochDep` does not forward `setPhiDependence` to its component
  `Aeff`s, so turning phi dependence off has no effect on a multi-epoch
  effective area. Its `m_usePhiDependence` is also never initialized.
- `EpochDep` requires epochs in ascending order of start time, which
  depends on the row order in `caldb.indx`.

## Caches and ownership

| Object | Cache | Invalidated by |
|--------|-------|----------------|
| `latResponse::Irfs` | lazily built components | never (one per `Irfs` object) |
| `Psf3` | `PsfIntegralCache` for one ROI cone | a different cone; `setParams`. Not copied by the copy constructor (the clone rebuilds it). |
| `Edisp3`/`Edisp2` | `EdispInterpolator`; last (log10 E, cos θ) parameters | never |
| `EpochDep` | last epoch index | — |
| `IPsf` (base) | static integrand state; static ψ grid | every integral call |

`Edisp3` and `Edisp2` have no user-defined copy constructor, so `clone()`
copies the raw `m_interpolator` pointer. If the original has been
evaluated (interpolator built), the clone and the original both delete it
(suspected double free, from reading the code; not checked at runtime).
Factory prototypes are never evaluated, so `IrfsFactory::create` is safe;
cloning an `Irfs` that has already been used is not.

## Other gotchas

- `Aeff::max_phi_modulation()` evaluates the modulation at `M_PI/4` and
  `M_PI/2`, but `phi_modulation` takes degrees. `upperLimit()`, which
  `observationSim` uses as the rejection-sampling envelope, may therefore
  underestimate the true maximum (suspected bug, from reading the code;
  not checked at runtime).
- `EfficiencyFactor::operator()(E, met)` leaves `ltfrac` uninitialized if
  `met` falls in a gap between FT2 intervals (suspected bug, from reading
  the code).
- `Psf3::findIndex` writes diagnostic numbers to `std::cout` before
  throwing.
- Fixed-size parameter buffers: `Edisp3`/`Edisp2` keep parameters in
  `double m_pars[10]`. A dispersion table with more than 10 parameter
  columns would overflow it.
- `latResponse/work/TH2_test.C` is a ROOT macro left from development.
- `src/combine_irfs/` builds no target; `pfiles/combine_irfs.par` is still
  installed.
- The `latResponse/data/` files (`aeff_epoch_*.fits`, `psf_epoch_*.fits`,
  `edisp_epoch_*.fits`, `psf_Pass5_v0_front.fits`) are test inputs only.

## Extension recipe: a new PSF parameterization

1. Write `PsfN : public PsfBase` (or `IPsf`), implementing both `value`
   overloads and `clone()`. Override `angularIntegral(E, θ, φ, radius)`
   and the cone version if an analytic integral exists; the base-class
   versions work but are slow and use static state.
2. Construct it from `(const irfUtil::IrfHdus &, size_t iepoch)` like
   `Psf3`, using `hdus("RPSF")` and `hdus("PSF_SCALING")`.
3. Add `case N:` to `IrfLoader::psf(hdus, iepoch)` and write `PSFVER = N`
   in the CALDB files.
4. If Python access is wanted, add the header to `pyIrfLoader.i` with a
   `cast` helper like the existing ones.
5. Add checks to `latResponse/src/test/test.cxx` (normalization,
   zero separation, ROI integral).
