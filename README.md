CMGG Structure experiments
================

## No boilerplate

- Viash expects you to write a script
  (e.g. [`src/wgd_score/script.R`](src/wgd_score/script.R)) and add some
  metadata
  (e.g. [`src/wgd_score/config.vsh.yaml`](src/wgd_score/config.vsh.yaml)).

- In return, Viash will generate a lot of boiler plate code for you!
  This includes:

  - Argument parsing
  - A custom Dockerfile for a given component
  - A standalone executable (run from CLI)
  - A Nextflow pipeline (run from CLI)
  - A Nextflow module (include in other Nextflow pipelines)

Looking at the [R script](src/wgd_score/script.R), we can see it becomes
a lot easier to read what the actual functionality of the component is.

We can let Viash build all components it can detect in the repo as
follows:

``` bash
viash ns build --setup cachedbuild
```

    Exporting wgd_score  =docker=> /home/rcannood/workspace/cmgg/cmgg-structural/target/docker/wgd_score
    [notice] Building container 'ghcr.io/rcannood/wgd_score:dev' with Dockerfile
    Exporting wgd_score  =nextflow=> /home/rcannood/workspace/cmgg/cmgg-structural/target/nextflow/wgd_score
    [32mAll 2 configs built successfully[0m

## Native support for Bash / R / Python scripts

Viash supports multiple scripting languages (Bash / Python / R) and
behaves the same way in every language.

## Viash generates standalone executables

By running the build command above, Viash generated a [Bash
executable](target/docker/wgd_score/wgd_score) which contains the
`wgd_score` script. Behind the screens, this script is actually being
run inside a Docker container. This means that you can send it to
someone, the script will automatically fetch the required Docker
container and everything should just work from there.

``` bash
target/docker/wgd_score/wgd_score --help
```

    wgd_score dev

    Arguments:
        --matrix
            type: file, required parameter, file must exist

        --scoring_mask
            type: file, required parameter, file must exist

        --scores
            type: file, required parameter, output, file must exist
            example: scores.txt.gz

        --plot
            type: file, output, file must exist
            example: plot.pdf

        --gzip
            type: boolean_true
            Gzip output files

        --outliers
            type: boolean_true
            label outliers on dosage score visualization

Example (not run due to lack of data):

``` bash
WDIR="/Users/rlc/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_analysis/WGD_version2_Dec2017/"

target/docker/wgd_score/wgd_score \
  --matrix "$WDIR/WGD_training.6F_adjCov.WGD_scoring_masked.100bp.matrix.bed.gz" \
  --scoring_mask "$WDIR/WGD_scoring_mask.6F_adjusted.100bp.h37.bed" \
  --scores "$HOME/scratch/WGDscore_testing/scores.txt.gz" \
  --plot "$HOME/scratch/WGDscore_testing/plot.pdf" \
  --gzip
```

## Each component has its own container

This allows to specify component-specific dependencies.

This component doesn’t need any additional dependencies, but I added
`funkyheatmap` as a dependency for demonstration purposes:

``` bash
target/docker/wgd_score/wgd_score ---dockerfile
```

    FROM rocker/tidyverse:4.2

    RUN Rscript -e 'if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")' && \
      Rscript -e 'remotes::install_cran(c("funkyheatmap"), repos = "https://cran.rstudio.com")'

    LABEL org.opencontainers.image.description="Companion container for running component wgd_score"
    LABEL org.opencontainers.image.created="2022-11-30T08:36:53+01:00"
    LABEL org.opencontainers.image.source="https://github.com/rcannood/cmgg-structural"
    LABEL org.opencontainers.image.version="dev"

## Viash generates Nextflow pipelines

Viash also generated a [Nextflow pipeline](target/nextflow/wgd_score/).

``` bash
nextflow run target/nextflow/wgd_score/main.nf --help
```

    N E X T F L O W  ~  version 21.10.6
    Launching `target/nextflow/wgd_score/main.nf` [cheesy_bernard] - revision: cf00380714
    wgd_score dev

    Arguments:
        --id
            type: string
            default: run
            A unique id for every entry.

        --matrix
            type: file, required parameter, file must exist

        --scoring_mask
            type: file, required parameter, file must exist

        --scores
            type: file, required parameter, output, file must exist
            default: $id.$key.scores.gz
            example: scores.txt.gz

        --plot
            type: file, output, file must exist
            default: $id.$key.plot.pdf
            example: plot.pdf

        --gzip
            type: boolean_true
            default: false
            Gzip output files

        --outliers
            type: boolean_true
            default: false
            label outliers on dosage score visualization

    Nextflow input-output arguments:
        Input/output parameters for Nextflow itself. Please note that both
        publishDir and publish_dir are supported but at least one has to be
        configured.

        --publish_dir
            type: string, required parameter
            example: output/
            Path to an output directory.

        --param_list
            type: string
            example: my_params.yaml
            Allows inputting multiple parameter sets to initialise a Nextflow
            channel. A `param_list` can either be a list of maps, a csv file, a json
            file, a yaml file, or simply a yaml blob.
            
            * A list of maps (as-is) where the keys of each map corresponds to the
            arguments of the pipeline. Example: in a `nextflow.config` file:
            `param_list: [ ['id': 'foo', 'input': 'foo.txt'], ['id': 'bar', 'input':
            'bar.txt'] ]`.
            * A csv file should have column names which correspond to the different
            arguments of this pipeline. Example: `--param_list data.csv` with
            columns `id,input`.
            * A json or a yaml file should be a list of maps, each of which has keys
            corresponding to the arguments of the pipeline. Example: `--param_list
            data.json` with contents `[ {'id': 'foo', 'input': 'foo.txt'}, {'id':
            'bar', 'input': 'bar.txt'} ]`.
            * A yaml blob can also be passed directly as a string. Example:
            `--param_list "[ {'id': 'foo', 'input': 'foo.txt'}, {'id': 'bar',
            'input': 'bar.txt'} ]"`.
            
            When passing a csv, json or yaml file, relative path names are
            relativized to the location of the parameter file. No relativation is
            performed when `param_list` is a list of maps (as-is) or a yaml blob.

You can run the pipeline in much the same way as before:

``` bash
WDIR="/Users/rlc/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_analysis/WGD_version2_Dec2017/"

nextflow run target/nextflow/wgd_score/main.nf \
  -profile docker \
  --matrix "$WDIR/WGD_training.6F_adjCov.WGD_scoring_masked.100bp.matrix.bed.gz" \
  --scoring_mask "$WDIR/WGD_scoring_mask.6F_adjusted.100bp.h37.bed" \
  --scores "scores.txt.gz" \
  --plot "plot.pdf" \
  --publish_dir "$HOME/scratch/WGDscore_testing/"
  --gzip
```

There are multiple profiles (docker, singularity, podman) which should
allow you to run the pipeline on HPC systems which don’t have Docker
installed.

## … Which are also Nextflow modules

``` groovy
include { wgd_score } from './target/nextflow/wgd_score/main.nf'

def wdir = "/Users/rlc/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_analysis/WGD_version2_Dec2017/"

workflow {
  data = [
    matrix: file("$wdir/WGD_training.6F_adjCov.WGD_scoring_masked.100bp.matrix.bed.gz"), 
    scoring_mask: file("$wdir/WGD_scoring_mask.6F_adjusted.100bp.h37.bed")
  ]

  Channel.value(["input", data])
    | wgd_score
}
```

## VDSL3 modules are dynamic

We refer to Nextflow modules generated by Viash as VDSL3 modules. This
is because they enable a completely different way of writing pipelines.

One of the main annoyances I have with standard Nextflow modules is that
they are very difficult to reuse across projects because they are quite
static. For instance, depending on the pipeline a component is part of,
it might be necessary to change some of the directives dynamically.

``` groovy
include { wgd_score } from './target/nextflow/wgd_score/main.nf'

def wdir = "/Users/rlc/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_analysis/WGD_version2_Dec2017/"

workflow {
  data = [
    matrix: file("$wdir/WGD_training.6F_adjCov.WGD_scoring_masked.100bp.matrix.bed.gz"), 
    scoring_mask: file("$wdir/WGD_scoring_mask.6F_adjusted.100bp.h37.bed")
  ]

  Channel.value(["input", data])
    | wgd_score.run(
      directives: [
        cache: false,
        label: ["bigmem", "bigcpu", "gpu"],
        maxRetries: 3
      ],
      auto: [
        publish: true
      ]
    )
}
```

You can also find out more information about the Viash config from
within the Nextflow pipeline:

``` groovy
include { wgd_score } from './target/nextflow/wgd_score/main.nf'

println("Fetch a field from the config: {wgd_score.config.functionality.name}")
```
