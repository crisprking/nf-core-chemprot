process CHEMPROT_ABPP {
    tag "$meta.id"
    label "process_medium"
    publishDir "${params.outdir}/chemprot/abpp", mode: "copy"

    input:
    tuple val(meta), path(msstats_file)
    val condition_col
    val fdr_thresh

    output:
    tuple val(meta), path("*.abpp_results.tsv"), emit: results
    path "versions.yml"                        , emit: versions

    script:
    def out_prefix = msstats_file.baseName
    """
    Rscript \$projectDir/bin/chemprot_abpp.R \        $msstats_file \
        "$condition_col" \
        $fdr_thresh \
        $out_prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(Rscript -e "cat(as.character(R.version.major, R.version.minor))")
    END_VERSIONS
    """
}
