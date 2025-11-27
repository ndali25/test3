Below is a **clean, machine-friendly, C#-parser-friendly README** designed so a program can reliably extract:

* Software version
* Required folders
* Execution order
* Docker instructions
* Script list
* Data files

I avoided long sentences, avoided ambiguous formatting, and used **strict key–value sections** with **stable tag markers** so your C# program can parse it deterministically.

Feel free to adjust tag names as needed.

---

# **PROJECT TITLE**

What drives the recent surge in inflation? The historical decomposition roller coaster

---

# **METADATA**

* LANGUAGE: MATLAB
* VERSION_REQUIRED: R2022b
* CPU_RECOMMENDED: quad-core
* RAM_RECOMMENDED_GB: 16
* DISK_MIN_GB: 2.5
* RUNTIME_FAST: 2h
* RUNTIME_SLOW: 48h

---

# **FOLDER_STRUCTURE**

```
/Data
/Functions
/Figures
/ (root scripts)
```

---

# **DATA_FILES**

* US bivariate.xlsx
* US large.xlsx
* EuroArea FRED.xlsx
* IMF data different countries.xlsx

---

# **SCRIPTS_EXECUTION_ORDER**

1. macro VAR diffuse.m
2. macro VAR diffuse 83 19.m
3. macro VAR *.m (all remaining macro VAR scripts)
4. simulation VAR*.m (all simulation scripts)

---

# **SCRIPT_TAGS**

* TYPE_MAIN: macro VAR diffuse.m, macro VAR diffuse 83 19.m
* TYPE_VARIANTS: macro VAR *.m
* TYPE_SIMULATION: simulation VAR*.m

---

# **REPRODUCIBILITY**

* RANDOM_SEED: fixed at start of every script
* INTERACTION_REQUIRED: false
* OUTPUT_AUTOSAVE: true
* OUTPUT_FOLDER: Figures

---

# **EXECUTION_INSTRUCTIONS**

1. Open MATLAB R2022b.
2. Set working directory to repository root.
3. Run scripts in **SCRIPTS_EXECUTION_ORDER**.
4. All figures saved automatically to `/Figures`.

---

# **DOCKER_INSTRUCTIONS**

The following instructions are designed for automatic generation by a C# program.

## **DOCKER_REQUIREMENTS**

* MATLAB_VERSION: R2022b
* MATLAB_RUNTIME_ALLOWED: true
* FOLDERS_TO_COPY: Data, Functions, Figures(optional), root scripts
* DISK_MIN_GB: 2.5

## **DOCKER_COMMANDS**

To run scripts inside Docker (executed sequentially):

```
matlab -batch "macro VAR diffuse"
matlab -batch "macro VAR diffuse 83 19"
matlab -batch "run_all_macro_VAR_variants"
matlab -batch "run_all_simulation_VAR"
```
---

# **PROGRAMMATIC_HINTS**

This section is intentionally structured to support easy parsing.

* TAG_START:DATA_FILES
* TAG_END:DATA_FILES
* TAG_START:SCRIPTS_EXECUTION_ORDER
* TAG_END:SCRIPTS_EXECUTION_ORDER
* TAG_START:DOCKER_REQUIREMENTS
* TAG_END:DOCKER_REQUIREMENTS

