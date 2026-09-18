# 2026_her.27.1-24a514a_assessment

2026 - Herring (Clupea harengus) in subareas 1, 2, 5 and divisions 4.a and 14.a, Norwegian spring-spawning herring (the Northeast Atlantic and Arctic Ocean) - WGWIDE (ADGWIDE)

## Project structure

The project are structured as a ices TAF project with some additions. See description underneath

📦 TAF Project

```         
── 📁 boot/ 
    ── 📄 DATA.bib
    ── 📄 SOFTWARE.bib
    ── 📁 initial/data/.      (input data for the assessment)
    ── 📁 rscript/.           (ansamble of utility script)
── 📄 boot.R
    -> prepeare the environment as specified in SOFTWARE.bib, and transfare data by DATA.bib
── 📄 data.R
    -> Read the data files
    -> Store relevant information in the folder data/
    -> Preparing the data for the SAM stockassessment model, and store this as .Rda in data/
── 📄 model.R
    -> Define the configuration
    -> Run SAM model
    -> Run dignostic
    -> Save model and dignostic output in the folder model/
── 📄 output.R
    -> Make figures and tables for dignostic purposes and store these in output/
── 📄 report.R
    -> Make figures and tables for WGWIDE report and store this in report/
── 📄 forecast.R (Not initially in the TAF official structure)
    -> Do the forecast
    -> Save outputs from the different scenarios as .Rda in the folder forecast/
    -> Save relevant figures in the folder forecast/figures
    -> Save relevant tables in the folder forecast/tables
── 📄 test.R (Not initially in the TAF official structure)
    -> Test if files generated localy matches those produced in TAF
── 📄 publish.R (Not initially in the TAF official structure)
    -> Make the presentations and store them in the folder publish/
```

📦 TAF Project artifacts

```
All files made by the infrastructure are available as artifakt, and is stored as: 
── 📁 data/            
    -> data files and data objects
── 📁 model/           
    -> model object 
── 📁 output/          
    -> Figures and tables
── 📁 report/          
    -> Official files for assessment report
── 📁 forecast/          
    -> Official files for forecast report
── 📁 pubslih/          
    -> Presentations relevant for the WG
```

## Utility scripts (`boot/rscript/`)

These are not TAF pipeline steps themselves, but are `source()`d by several of the
scripts above. Kept separate from the pipeline structure since they are shared
building blocks, not something you run directly:

- **`utils.r`** - general-purpose helpers sourced by *every* script in the
  pipeline (`data.R`, `model.R`, `output.R`, `report.R`, `forecast.R`):
  observation weighting (`computeExternalWeight()`), and the leave-one-out /
  retrospective / fit-comparison plotting functions used for diagnostics.
- **`forecast_utils.R`** - forecast-specific helpers, sourced by `data.R` and
  `forecast.R`: the density-dependent maturity prediction wrapper
  (`get_mopred()`), the scenario runner (`run_F_scenario()`), the catch-options
  table builder (`build_output_row()`), and the functions behind each forecast
  figure (`plot_fig_4_8_1_1()`, `plot_fig_4_8_1_2()`, etc.).
- **`DD_model.r`** - the cohort-based density-dependent maturity-at-age model
  (Vatnehol et al. 2024, ICES JMS 81(10): 2028-2036), used by `get_mopred()` to
  predict maturity for the most recent (back-calculation-lagged) years and the
  forecast years.
- **`SAM_forecast_from_GIT.R`** - a vendored copy of the SAM package's forecast
  engine (`forecast_XSAM()`), pulled from
  [fishfollower/SAM](https://github.com/fishfollower/SAM) on GitHub rather than
  installed as a package. Not written for this assessment - leave as-is unless
  specifically updating it against upstream.
- **`WGWIDE2025.Rda`, `fit_official_benchmark.Rda`, `forecast_2024.Rda`,
  `forecast_2025.Rda`** - saved fit/forecast objects from previous assessment
  years, loaded for year-over-year comparison figures in `output.R`/
  `forecast.R` (e.g. this year's fit vs. last year's).

## Updating data for a new assessment year

When starting a new assessment year, update things in this order:

1. **Drop the new/updated data files into `boot/initial/data/`.**
   Keep the same filenames the pipeline already expects — `caa.txt`, `cw.txt`,
   `matprop.txt`, `nm.txt`, `pf.txt`, `pm.txt`, `west.txt`, the survey and RFID
   files, `prop_mat_new.txt` — including `survey_T2025-02-05.txt`, the
   multi-fleet survey file, despite the date in its name.

2. **Register every new/changed file in `boot/DATA.bib`.**
   Either by hand, or by re-running `prepare_bibfiles.R` (update the
   `year`/`period` arguments in there first). This step is easy to forget —
   double-check the `period` field actually matches the file's real year
   range before committing, since a stale `DATA.bib` silently misrepresents
   what data the assessment used.

3. **Add any newly-used R package to `boot/SOFTWARE.bib`.**
   Via `prepare_bibfiles.R`'s `TAF::draft.software()` calls, whenever a
   script starts using a package that wasn't there before. `SOFTWARE.bib` is
   what `taf.boot()` uses to install dependencies on a fresh checkout — a
   package only present in your own local R library won't be installed for
   anyone else unless it's listed here.

4. **Run `boot.R`** (`TAF::taf.boot()`) to copy the new files from
   `boot/initial/data/` into `boot/data/` and refresh the bootstrap.

5. **Run the pipeline in order:**
   `data.R` → `model.R` → `output.R` → `report.R` → `forecast.R` → `test.R` →
   `publish.R`.

### Things worth double-checking whenever the underlying data changes

- **Survey fleet order** — `data.R`'s `surveys` list, mirrored in `model.R`'s
  `conf$key...` matrices and in `output.R`/`report.R`'s fleet-name lookups —
  is assumed fixed, not re-derived. If a fleet is ever added, removed or
  reordered, all three places need updating together.
- **`test.R`** compares a fresh run against whatever `data/`, `model/`,
  `output/`, `report/` and `forecast/` looked like at the last commit. After a
  real, intentional change to data or code, the reproducibility check will
  fail until you commit the newly produced output folders as the new
  reference.

---

## License

Released under the [MIT License](LICENSE).

## AI assistance

Parts of this repository's code were structured, refactored, documented and
reviewed with the assistance of AI tooling (Claude Code). All scientific
methods, model configuration and final results remain the responsibility of
the assessment authors.
