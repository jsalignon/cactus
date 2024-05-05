
# jsalignon/cactus: Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

<!-- Guiding Principles

    Changelogs are for humans, not machines.
    There should be an entry for every single version.
    The same types of changes should be grouped.
    Versions and sections should be linkable.
    The latest version comes first.
    The release date of each version is displayed.
    Mention whether you follow Semantic Versioning.

Types of changes

    Added for new features.
    Changed for changes in existing functionality.
    Deprecated for soon-to-be removed features.
    Removed for now removed features.
    Fixed for any bug fixes.
    Security in case of vulnerabilities. -->

## [v1.0.0](https://github.com/jsalignon/cactus/releases/tag/v1.0.0) - Aztekium, 05-05-2024

### Added
- [#21, commit 021de6f](https://github.com/jsalignon/cactus/commit/021de6f05c25ed834aee493d3128de774c79dff7) - Added PA filters to keep peaks at less or more than 10 or 100 kilobases (lt10kb, lt100kb, mt10kb, mt100kb).
- [#21, commit 021de6f](https://github.com/jsalignon/cactus/commit/021de6f05c25ed834aee493d3128de774c79dff7) - Added PA filters to keep peaks at less than arbitrary defined cutoffs (ltXkb, mtYkb), with cutoffs defined with the new parameters custom_distance__less_than_X_b and custom_distance__more_than_Y_b.

### Changed

### Fixed
- [#21](https://github.com/jsalignon/cactus/pull/21); commits [5ac9710](https://github.com/jsalignon/cactus/commit/5ac9710ebac7341dd378db5cbfa0ec187e413b13#diff-6401496ba455b9488ffa902a6e4d7732b2c60ff2d77c5c3ef96b28a7ac7d3b28L2772) and [827df11](https://github.com/jsalignon/cactus/commit/827df11f8af185213762426fc8c9cbd08e1e15a2#diff-6401496ba455b9488ffa902a6e4d7732b2c60ff2d77c5c3ef96b28a7ac7d3b28L2775) - Fixed a small typo that was making the previous release (0.9.0) crash during the DA_ATAC__saving_detailed_results_tables process.

### Deprecated

### Removed
- [#21](https://github.com/jsalignon/cactus/pull/21) - Removed the PA filters: 3kb, 8kb, and 30kb.

### Dependencies


## [v0.9.0](https://github.com/jsalignon/cactus/releases/tag/v0.9.0) - Gymnocalycium, 22-04-2024

### Added

 - [#20](https://github.com/jsalignon/cactus/pull/20) - Added a tutorial section in the documentation and an associated script.
 - [#20](https://github.com/jsalignon/cactus/pull/20) - Added figures for the manuscript: genome track plots, examples of figures and tables from the worm test dataset.
 - [#20, commit 8ab3377](https://github.com/jsalignon/cactus/commit/8ab33772fd9493a542478ebc1746ae3eb9af0cc3) - Added Peak Annotation filters to keep only peaks at less than 3, 8 or 30 kilobases to their closest gene.
 - [#20, commit 7f797f8](https://github.com/jsalignon/cactus/commit/7f797f892db6d83f9409eea33d5cf658ebe6ad4b) - Added ChIPseeker options to ignore upstream or downstream of peak for closest gene annotation.
 - [#19](https://github.com/jsalignon/cactus/pull/19) - Added quickstart scripts.

### Changed

### Fixed

 - [#20, commit 250afb6](https://github.com/jsalignon/cactus/commit/250afb6174c3e976ec9eda72da480732daf8a938#diff-6401496ba455b9488ffa902a6e4d7732b2c60ff2d77c5c3ef96b28a7ac7d3b28) - Fixing a typo which prevented bigWig files from being saved.

### Deprecated

### Removed

### Dependencies


## [v0.8.6](https://github.com/jsalignon/cactus/releases/tag/v0.8.6) - Pygmaeocereus, 30-11-2023

### Added

 - [#18](https://github.com/jsalignon/cactus/pull/18) - Added the code for the case study.

### Changed

### Fixed

### Deprecated

### Removed

### Dependencies


## [v0.8.5](https://github.com/jsalignon/cactus/releases/tag/v0.8.5) - Harrisia, 28-11-2023

### Added

### Changed

- Changed the code to make Cactus work with the latest versions of Singularity, Docker, conda and Mamba. [commit](https://github.com/jsalignon/cactus/commit/c8757995185317b0a828f61375bdb7605ce111e4)
- Changed the documentation of the Quickstart and dependencies. [567e048](https://github.com/jsalignon/cactus/commit/567e048a75e83ee062c9ec32996693ce24b165e8) [f8b2027](https://github.com/jsalignon/cactus/commit/f8b2027b105f823ef7a460a17677e190388f056c)
- Changed the documentation of the new containers. [e89cb32](https://github.com/jsalignon/cactus/commit/e89cb32f42d833a0352cb1a570778b7729aee461)
- Changed the Nextflow version manifest to give an error message if a wrong version is used. [51981f3](https://github.com/jsalignon/cactus/commit/51981f39e3db3050c8a2958767cd89269dc578b5)

### Fixed

- Fixed the conda and mamba bugs from the previous release by creating new mulled containers. [c875799](https://github.com/jsalignon/cactus/commit/c8757995185317b0a828f61375bdb7605ce111e4) [6cfcf16](https://github.com/jsalignon/cactus/commit/6cfcf167a820bfc5f2e21b69fbdf81156b95c4db)
- Fixed bugs to run the main.nf and the download.nf scripts with Docker. [716adc2](https://github.com/jsalignon/cactus/commit/716adc222acad0a9587f194461410f8943a678c7) [c947936](https://github.com/jsalignon/cactus/commit/c94793649e82172b1d3fcfcc4a067b7cc748d3c9)

### Deprecated

### Removed

### Dependencies
 - upgraded dependencies: Docker 24.0.5, build 24.0.5-0ubuntu1~20.04., conda: 23.7.4, Mamba: 1.4.2. [567e048](https://github.com/jsalignon/cactus/commit/567e048a75e83ee062c9ec32996693ce24b165e8)



## [v0.8.4](https://github.com/jsalignon/cactus/releases/tag/v0.8.4) - Matucana, 09-06-2023

### Added

### Changed

- Changed the code to make Cactus work with the latest versions of Nextflow, Singularity and Docker. [commit](https://github.com/jsalignon/cactus/commit/ce2b6dd3d8ca8beee5479e85a4d24e4c4b022641)
- Changed the documentation to clarify the dependencies. [commit](https://github.com/jsalignon/cactus/commit/0f45173c61a78dd8ce742ff6413c90d1c699575e)
- Changed the CITATIONS.md and the README.md to add a link to Cactus preprint. [commit](https://github.com/jsalignon/cactus/commit/0e7b6ba4cb4269708897f8ad27adc57048dc9229)

### Fixed

- Fixed a bug happening in conda/Mamba when users have a .Renviron or a .Rprofile file. [commit](https://github.com/jsalignon/cactus/commit/a188dbb31dfe34547aaab427175dde28237ac36a)
- Fixed a bug with conda/Mamba for merging pdf. [commit](https://github.com/jsalignon/cactus/commit/a188dbb31dfe34547aaab427175dde28237ac36a)

### Deprecated

- Deprecated the use of conda and Mamba temporarily until a bug in conda is resolved. [commit](https://github.com/jsalignon/cactus/commit/0f45173c61a78dd8ce742ff6413c90d1c699575e)

### Removed

### Dependencies



## [v0.8.3](https://github.com/jsalignon/cactus/releases/tag/v0.8.3) - Espostoopsis, 15-05-2023

### Added

### Changed

- Changed the README.md file to add a link to the preprint and to the LICENCE: [1](https://github.com/jsalignon/cactus/commit/c6cae1efa6a1d6d94a8fe31d4f1a3c73a8046b26), [2](https://github.com/jsalignon/cactus/commit/b4ae5c359f6aea80ab236c5064372742e76ff34e)

### Fixed

- [#11](https://github.com/jsalignon/cactus/pull/11) - Fixed a bug in the heatmap process that made the fly test dataset fail.

### Deprecated

### Removed

### Dependencies



## [v0.8.2](https://github.com/jsalignon/cactus/releases/tag/v0.8.2) - Mammilloydia, 12-05-2023

### Added

- [#6](https://github.com/jsalignon/cactus/pull/6) - Added options to keep unique DA and NDA genes or not in the splitting process

### Changed

- [#5](https://github.com/jsalignon/cactus/pull/5) - Changed the docs

### Fixed

- [#8](https://github.com/jsalignon/cactus/pull/8) - Fixed a bug in the script to download test datasets and references

### Deprecated

### Removed

### Dependencies


## [v0.8.1](https://github.com/jsalignon/cactus/releases/tag/v0.8.1) - Cephalocereus, 19-04-2023

### Added

 - [adding a functionality to manually set the size of the numbers indicating the overlap in the cells for the heatmap figures](https://github.com/jsalignon/cactus/commit/bf0f8aa2925d03760d54c35475415bae538625a3)

### Changed

 - Updated the docs to match with latest version of manuscript in preparation: [1](https://github.com/jsalignon/cactus/commit/bc11d3a39c05686c9dc6ceb7555539334043ce70), [2](https://github.com/jsalignon/cactus/commit/ac0eafa643f85de6efa2cbd8d1aae62059683325), [3](https://github.com/jsalignon/cactus/commit/f2739517621504c756f1304139dab01b99f0ecae), [4](https://github.com/jsalignon/cactus/commit/cf934e22270e7a18122ea3ad44a0f570e8610e53), [5](285d579023f7cbac7771405854d17508030420a4), [6](https://github.com/jsalignon/cactus/commit/285d579023f7cbac7771405854d17508030420a4), [7](https://github.com/jsalignon/cactus/commit/f04cd11eb6b5992f771e72364c019e7f1c251b68), [8](https://github.com/jsalignon/cactus/commit/f68018bf82bd5165056f9cfaa26b3e2f29e1cb08), [9](https://github.com/jsalignon/cactus/commit/37ad16c21c79ca4674c11de84254e65019605f8b)

 - [making clearer labels for the heatmap](https://github.com/jsalignon/cactus/commit/1fd6d2054839fa4978d8c8e2f765c90611f36b02)

### Fixed

 - [fixing bugs in the heatmap process](https://github.com/jsalignon/cactus/commit/eb3de68366fd21c90220cb2a9155af06da04da55)

 - [fixing a bug with the plot_FDR_by_PA_filters function](https://github.com/jsalignon/cactus/commit/4e5cafc61b54958f2144f6f3be4637611e5b8de7)


### Deprecated

### Removed

### Dependencies



## [v0.8.0](https://github.com/jsalignon/cactus/releases/tag/v0.8.0) - Lophophora, 30-12-2022

First release.

<!-- ## [2021.11.16](https://github.com/nf-core/sarek/releases/tag/2.7.1) - Pårtejekna

Pårtejekna is one of glaciers of the Pårte Massif. 

check out this Changelog for a formatting example.
https://raw.githubusercontent.com/veeso/ATtila/main/CHANGELOG.md
-->


