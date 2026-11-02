//
// Subworkflow with functionality specific to the biobakery/hvp pipeline.
//
// Laid out as the nf-core pipeline template's utils_nfcore_<name>_pipeline is,
// and built on the standard nf-core subworkflows (subworkflows/nf-core/, pinned
// in modules.json, updated with `nf-core subworkflows update`):
//
//   UTILS_NEXTFLOW_PIPELINE   --version, params dump to pipeline_info/, conda checks
//   UTILS_NFSCHEMA_PLUGIN     --help / --help_full, parameter summary, validation
//                             against nextflow_schema.json (nf-schema plugin)
//   UTILS_NFCORE_PIPELINE     config / profile sanity checks, completionSummary
//
// plus what is HVP's own: the input read channel and the component registry.

include { UTILS_NFSCHEMA_PLUGIN   } from '../../nf-core/utils_nfschema_plugin'
include { UTILS_NFCORE_PIPELINE   } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE } from '../../nf-core/utils_nextflow_pipeline'
include { completionSummary       } from '../../nf-core/utils_nfcore_pipeline'
include { read_input              } from '../read_input/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    readsdir          //  string: Folder of fastq files
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    //
    // Print version and exit if required, and dump pipeline parameters to JSON
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and print the parameter summary
    //
    def before_text = "\nbiobakery/hvp ${workflow.manifest.version}\n"
    def after_text  = "\nComponents: conf/components.json  ·  Architecture: docs/architecture.md\n"
    def command     = "nextflow run biobakery/HVP-nextflow -profile <site> --readsdir <fastqs> --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Create the read channel from the input folder
    //
    if (!readsdir) {
        error "ERROR: --readsdir is required. Example: --readsdir /path/to/fastqs"
    }
    reads = read_input(readsdir, 'hvp')

    emit:
    reads // channel: [ [id: sample, paired_end: bool], reads ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    monochrome_logs // boolean: Disable ANSI colour codes in log output

    main:
    workflow.onComplete {
        completionSummary(monochrome_logs)
    }

    workflow.onError {
        log.error "Pipeline failed. See .nextflow.log, and the work directory of the failed task."
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPONENT REGISTRY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// conf/components.json lists every component, the stage it belongs to,
// whether it is implemented, and any Nextflow plugins it needs. It is read
// here at launch, and by bin/make_diagrams.py to draw the architecture
// diagram, so the two cannot disagree about what exists.

def loadRegistry() {
    def path = "${workflow.projectDir}/conf/components.json"
    return new groovy.json.JsonSlurper().parse(new File(path))
}

// Turn --components / --only into the list of components to run.
//
//   --components a,b,c   run exactly these (the default is every ready one)
//   --only x             run component x on its own; the preprocessing
//                        components named in --components still run first
//
// Either way, a component whose input nothing selected produces gets the
// default component of the stage that produces it added, with a log line: e.g.
// --only genomad needs contigs, so the 'assembly' component is added.
def resolveComponents(String components_param, String only) {
    def registry = loadRegistry()
    def known    = registry.components.collectEntries { [(it.name): it] }
    def stages   = registry.stages.collectEntries { [(it.id): it] }

    def asked = (components_param ?: '').tokenize(',')*.trim().findAll { it }
    if (only) {
        asked = asked.findAll { known[it]?.stage == 'preprocess' } + [only.trim()]
    }

    def unknown = asked.findAll { !known.containsKey(it) }
    if (unknown) {
        error "ERROR: unknown component(s): ${unknown.join(', ')}\n" +
              "Known components: ${known.keySet().join(', ')}"
    }

    def selected = asked.unique() as List

    // Is this data type produced, given the current selection? A chain stage
    // with nothing selected passes its input through, so it produces whenever
    // its own input is produced.
    def produced
    produced = { String data ->
        if (data == 'raw_reads') return true
        def stage = registry.stages.find { it.produces == data }
        if (selected.any { known[it].stage == stage.id }) return true
        return stage.mode == 'chain' && produced(stage.consumes)
    }

    def changed = true
    while (changed) {
        changed = false
        for (name in new ArrayList(selected)) {
            def needs = stages[known[name].stage].consumes
            if (produced(needs)) continue
            def producer = registry.stages.find { it.produces == needs }
            if (!producer.default) {
                error "ERROR: component '${name}' needs ${registry.data[needs]}, and no component produces it."
            }
            log.info "[hvp] '${name}' needs ${registry.data[needs]}: adding '${producer.default}' (${producer.title})."
            selected << producer.default
            changed = true
        }
    }

    // Report in pipeline order, and say plainly which ones are placeholders.
    def ordered = registry.components*.name.findAll { it in selected }
    log.info "[hvp] Components: ${ordered.join(' -> ')}" + (only ? "  (--only ${only})" : '')
    def skeletons = ordered.findAll { known[it].status == 'skeleton' }
    if (skeletons) {
        log.warn "[hvp] Placeholder component(s) selected: ${skeletons.join(', ')}. They are wired in " +
                 "but have no tool behind them yet, so they produce no output."
    }
    return ordered
}
