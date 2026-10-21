# Provenance: vendored/KaryoScope-ISCN/src (KaryoScope-ISCN analysis modules)

Source: barthel-lab/KaryoScope-ISCN at commit `5507369` (2026-08-12), on branch
`feat/refine-pivot-and-k-sweep`, merged into `main` on 2026-09-11 as `fd2aab0`.
The A6 Robertsonian parameter sweep and the A2 step 3 k-sweep ran against this commit.

Files are copied unmodified (sha256 below). They are the modules imported by
`revision/a6-parameter-sweeps/rob/scripts/sweep_rob.py` and
`revision/a2-k-sensitivity/step3/scripts/{sweep_refine,binsweep_analyse}.py`
(`core.bedio`, `core.cumulative`, `core.detect`, `core.refine`, `reference`), the modules
those import in turn (found by importing them and listing every `karyoscope_iscn` module
loaded), and the `chm13v2.0` reference tables that `reference.load_reference("chm13v2.0")`
reads. They use only the Python standard library.

The pin matters: `core/cumulative.py` here defines `DEFAULT_WINDOW_PAD`, which
`sweep_rob.py` imports. It was removed from KaryoScope-ISCN on 2026-09-11 (`75a9ebf`), so
the sweep fails against later versions.

Set `KS_ISCN_SRC` to this `src/` directory before running those scripts.

Checked 2026-10-20 against this copy, with a filled ROB configuration: `sweep_rob.py`
reproduces `rob_param_sweep.tsv` byte for byte, and `sweep_refine.py` reproduces the
k = 31 rows of `sst1_breakpoints_by_k.tsv` for all three samples.

One pin serves both analyses: the A2 k = 31 chain run with `karyoscope-iscn analyze` at this
commit (version 0.2.1) for the baseline and these modules for `sweep_refine.py` gives the same
three k = 31 rows. The scaffold map that `analyze` writes is the output of `karyoscope scaffold
--mode bed --feature-set region` (public KaryoScope); run directly with the same input label it is
byte-identical on all three samples.

Copyright: the KaryoScope-ISCN authors (from `git log` of these files at `5507369`:
Rhyker Ranallo-Benavidez (62 commits)). Distributed here under GPL-3.0 (see `../LICENSE`).

| file | sha256 |
|---|---|
| `karyoscope_iscn/__init__.py` | `b4eea860bef3f027e90dad26c832d8bd3b557c932409fb8a6da48ac420d61010` |
| `karyoscope_iscn/__pycache__/__init__.cpython-312.pyc` | `5081faca34c54fc62f1efa344defce0adf8e578ead453091714bb668192cdca9` |
| `karyoscope_iscn/__pycache__/_version.cpython-312.pyc` | `aef54261c61a4685b53fb373973630c5916d22614e26b46f06111983bf654bb1` |
| `karyoscope_iscn/__pycache__/exceptions.cpython-312.pyc` | `be3a3abdd3959aea91cb979d257a1ae86d1bc1aa2645bd881a7655718f8cd0ed` |
| `karyoscope_iscn/__pycache__/reference.cpython-312.pyc` | `26da5646ba55587f2ef66227e2fd62f3e03c43696ec6050c274a6016ece3ca1e` |
| `karyoscope_iscn/_version.py` | `dfce7115d09e4b35aa4ca8ea044ed974d9b11eaa91a454ab3333b39347712399` |
| `karyoscope_iscn/core/__init__.py` | `51db4ea80e9cdfae264ba9d343425d23941b92d5aaaeff8081edb42f9867b4bf` |
| `karyoscope_iscn/core/__pycache__/__init__.cpython-312.pyc` | `52245f22951d6ca5d6591275fccd1602a1d8d81d10718aa56b39f6ed2886256a` |
| `karyoscope_iscn/core/__pycache__/bedio.cpython-312.pyc` | `401b0d8ff52525b452acd7c6cd1bf346732ef86892db3a7a662813038219e417` |
| `karyoscope_iscn/core/__pycache__/cumulative.cpython-312.pyc` | `97a2cf35254b08b9b4cf5fddbf0408ac76bf97bffc2fc715a8d7a80e63ea73b9` |
| `karyoscope_iscn/core/__pycache__/cytoband.cpython-312.pyc` | `bb863531f30d8244607a61d8c111f4228cb0343509ada60c0969e9080bed00da` |
| `karyoscope_iscn/core/__pycache__/detect.cpython-312.pyc` | `92be1bba12dc78c1beb2ccc975775c04089e94fbf1fb5a949b2bcb9a1d236165` |
| `karyoscope_iscn/core/__pycache__/refine.cpython-312.pyc` | `6b2a3d80e180e39c5fd4e6f47dd68fccae86138c76665ef71b90f75ae7fe4b89` |
| `karyoscope_iscn/core/__pycache__/segments.cpython-312.pyc` | `e58315a75aeac2c5cd5e0e38695ea1284530f5eb8d680bd95fc653dcd9a18e2c` |
| `karyoscope_iscn/core/bedio.py` | `6a506e66d9c765402876702472c71d2488e4693b7571d8d8b4c9229bebc8eccc` |
| `karyoscope_iscn/core/cumulative.py` | `e532db4ed02ead7fd63563c408ac539fa52f96d6fd09d6d69864256378604a37` |
| `karyoscope_iscn/core/cytoband.py` | `0c94f136dc11c358bc5ed0a42238624c4a943c089ec1e890f1ce58c18d327897` |
| `karyoscope_iscn/core/detect.py` | `cd1d7ac45f2320a365373ccb4997151083f9c887a5ca7d6a70fe26f8ebc030b0` |
| `karyoscope_iscn/core/refine.py` | `1c52569a2c98a5fa1fc14b7bba0399da92b60c1e8ce627cdc45843cc75038922` |
| `karyoscope_iscn/core/segments.py` | `f42fbb882eb2649eda247abab991b294ab210f01dc6b8a0182a0c4380863e47d` |
| `karyoscope_iscn/data/references/chm13v2.0/chromosomes.tsv` | `c714012c59ece7455c108555f32838f66c0d489b6d0d2e27cde9d390c9f16e23` |
| `karyoscope_iscn/data/references/chm13v2.0/cytobands.bed` | `bc4272f4a866717e4d155de3856e969721be41772ebaf324555061a2c522aafd` |
| `karyoscope_iscn/data/references/chm13v2.0/mappability.tsv` | `4576fd07bb21152ca84616b97ce23b0f374ee28afc7a3daf94a442d15284e382` |
| `karyoscope_iscn/exceptions.py` | `b2c8e750ec679c2e8ad56a8156790467ea4557a4bbce5b3e34f84ad8d8528d10` |
| `karyoscope_iscn/reference.py` | `25a59ed2592d6dfab817ba0e9f05f38b617df44afb39b98b81fafa10a43f8a37` |
