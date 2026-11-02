#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Pipeline plumbing around the orchestrator, in the shape nf-core pipelines use
// (PIPELINE_INITIALISATION / PIPELINE_COMPLETION), plus the component registry.
//
// Parameter validation, --help and the parameter summary come from the nf-schema
// plugin and nextflow_schema.json rather than hand-written checks.

include { validateParameters } from 'plugin/nf-schema'
include { paramsSummaryLog   } from 'plugin/nf-schema'
include { read_input         } from '../read_input/main.nf'

// ── Initialisation: validate params, build the input channel ─────────────────
workflow PIPELINE_INITIALISATION {

    take:
    validate_params // boolean: run nf-schema validation
    readsdir        // string:  folder of fastq files

    main:
    if (validate_params) {
        validateParameters()
    }
    log.info paramsSummaryLog(workflow)

    if (!readsdir) {
        error "ERROR: --readsdir is required. Example: --readsdir /path/to/fastqs"
    }

    reads = read_input(readsdir, 'hvp')

    emit:
    reads // Channel: [ [id: sample, paired_end: bool], reads ]
}

// ── Completion: one line saying where the results went ───────────────────────
workflow PIPELINE_COMPLETION {

    main:
    workflow.onComplete {
        if (workflow.success) {
            log.info "[hvp] Completed. Results: ${params.outdir}"
        } else {
            log.error "[hvp] Failed. See .nextflow.log and the work directory of the failed task."
        }
    }
}

// ── Component registry ──────────────────────────────────────────────────────
// conf/components.json lists every component, the stage it belongs to and
// whether it is implemented. It is read here at launch, and by
// bin/make_diagrams.py to draw the architecture diagram, so the two cannot
// disagree about what exists.

def loadRegistry() {
    def path = "${workflow.projectDir}/conf/components.json"
    return new groovy.json.JsonSlurper().parse(new File(path))
}

// Turn --components / --plugin into the list of components to run.
//
//   --components a,b,c   run exactly these (the default is every ready one)
//   --plugin x           run component x on its own; the preprocessing
//                        components named in --components still run first
//
// Either way, a component whose input nothing selected produces gets the
// default component of the stage that produces it added, with a log line: e.g.
// --plugin genomad needs contigs, so the 'assembly' component is added.
def resolveComponents(String components_param, String plugin) {
    def registry = loadRegistry()
    def known    = registry.components.collectEntries { [(it.name): it] }
    def stages   = registry.stages.collectEntries { [(it.id): it] }

    def asked = (components_param ?: '').tokenize(',')*.trim().findAll { it }
    if (plugin) {
        asked = asked.findAll { known[it]?.stage == 'preprocess' } + [plugin.trim()]
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
    log.info "[hvp] Components: ${ordered.join(' -> ')}" + (plugin ? "  (plugin mode: ${plugin})" : '')
    def skeletons = ordered.findAll { known[it].status == 'skeleton' }
    if (skeletons) {
        log.warn "[hvp] Skeleton component(s) selected: ${skeletons.join(', ')}. They are wired in " +
                 "but have no tool behind them yet, so they produce no output."
    }
    return ordered
}

// Called by every skeleton component, so a run says what it skipped.
def skeletonNotice(String name) {
    log.warn "[hvp] ${name}: skeleton component, nothing to run yet (see docs/architecture.md, 'Adding a component')."
}
