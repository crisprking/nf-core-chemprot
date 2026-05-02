process CHEMPROT_TPP {
    tag "$meta.id"
    label "process_medium"
    publishDir "${params.outdir}/chemprot/tpp", mode: "copy"

    input:
    tuple val(meta), path(msstats_file)
    val temperatures
    val min_r2
    val min_peptides

    output:
    tuple val(meta), path("*.tpp_results.tsv"), emit: results
    path "versions.yml"                       , emit: versions

    script:
    def out_prefix = msstats_file.baseName
    """
    Rscript \$projectDir/bin/chemprot_tpp.R \        $msstats_file \
        "$temperatures" \
        $min_r2 \
        $min_peptides \
        $out_prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(Rscript -e "cat(as.character(R.version.major, R.version.minor))")
    END_VERSIONS
    """
}
