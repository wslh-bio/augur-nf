#!/usr/bin/env nextflow

//Description: Nextflow implementation of Nextstrain's Augur pipeline
//Available:
//Authors of this Nextflow: Abigail Shockey
//Email: abigail.shockey@slh.wisc.edu

// Input channels

Channel
    .fromPath( "${params.reference}")
    .into { reference_alignment; reference_translate }

Channel
    .fromPath( "${params.metadata}")
    .ifEmpty { exit 1, "Cannot find metadata file in: ${params.metadata}" }
    .into { filter_metadata; refine_tree_metadata; traits_metadata; metadata_export }

Channel
    .fromPath( "${params.colors}")
    .set { colors_export }

Channel
    .fromPath( "${params.lat_long}")
    .set { lat_long_export }

Channel
    .fromPath( "${params.auspice_config}")
    .set { config_export }

// Filter sequences if sequences to drop are included
if (params.dropped) {
    Channel
        .fromPath( "${params.sequences}")
        .ifEmpty { exit 1, "Cannot find sequences in: ${params.sequences}" }
        .set { sequences }
        
    Channel
        .fromPath( "${params.dropped}")
        .ifEmpty { exit 1, "Cannot find reference sequence in: ${params.reference}" }
        .set { dropped_strains }

    process filter{
      publishDir "${params.outdir}/filter", mode:'copy'
    
      input:
      file(fasta) from sequences
      file(metadata) from filter_metadata
      file(exclude) from dropped_strains
    
      output:
      file "filtered.fasta" into input_sequences
    
      shell:
      """
      augur filter \
        --sequences ${fasta} \
        --metadata ${metadata} \
        --output filtered.fasta \
        --group-by country year month \
        --sequences-per-group ${params.seq_per_group} \
        --min-date ${params.min_date}
      """
    }
}

else {
    Channel
        .fromPath( "${params.sequences}")
        .ifEmpty { exit 1, "Cannot find sequences in: ${params.sequences}" }
        .set { input_sequences }
}

process align{
  publishDir "${params.outdir}/align", mode:'copy'

  input:
  file(filtered_fasta) from input_sequences
  file(ref) from reference_alignment

  output:
  file "aligned.fasta" into raw_tree_alignment, refine_tree_alignment, ancestral_alignment

  shell:
  """
  augur align \
    --sequences ${filtered_fasta} \
    --reference-sequence ${ref} \
    --output aligned.fasta \
    --fill-gaps
  """
}

process tree{
  publishDir "${params.outdir}/tree", mode:'copy'

  input:
  file(msa) from raw_tree_alignment

  output:
  file "raw_tree.newick" into raw_tree

  shell:
  """
  augur tree \
    --alignment ${msa} \
    --output raw_tree.newick
  """
}

process refine{
  publishDir "${params.outdir}/refine", mode:'copy'

  input:
  file(tree) from raw_tree
  file(msa) from refine_tree_alignment
  file(metadata) from refine_tree_metadata

  output:
  file "refined_tree.newick" into refined_tree_traits, refined_tree_ancestral, refined_tree_translate, refined_tree_export
  file "branch_lengths.json" into branch_lengths_export

  shell:
  """
  augur refine \
    --tree ${tree} \
    --alignment ${msa} \
    --metadata ${metadata} \
    --output-tree refined_tree.newick \
    --output-node-data branch_lengths.json \
    --timetree \
    --coalescent opt \
    --date-confidence \
    --date-inference marginal \
    --clock-filter-iqd 4
  """
}

process ancestral{
  publishDir "${params.outdir}/ancestral_nt", mode:'copy'

  input:
  file(tree) from refined_tree_ancestral
  file(msa) from ancestral_alignment
  output:
  file "nt_muts.json" into ancestral_nt_translate, ancestral_nt_export

  shell:
  """
  augur ancestral \
    --tree ${tree} \
    --alignment ${msa} \
    --output nt_muts.json \
    --inference joint
  """
}

process translate{
  publishDir "${params.outdir}/ancestral_aa", mode:'copy'

  input:
  file(tree) from refined_tree_translate
  file(nt_muts) from ancestral_nt_translate
  file(ref) from reference_translate

  output:
  file "aa_muts.json" into ancestral_aa_muts

  shell:
  """
  augur translate \
    --tree ${tree} \
    --ancestral-sequences ${nt_muts} \
    --reference-sequence ${ref} \
    --output aa_muts.json
  """
}

process traits{
  publishDir "${params.outdir}/traits", mode:'copy'

  input:
  file(tree) from refined_tree_traits
  file(metadata) from traits_metadata

  output:
  file "traits.json" into traits_export

  when:
  params.traits == true

  shell:
  """
  augur traits \
    --tree ${tree} \
    --metadata ${metadata} \
    --output traits.json \
    --columns region country \
    --confidence
  """
}

if (params.traits) {
    process export{
      publishDir "${params.outdir}", mode:'copy'
    
      input:
      file(tree) from refined_tree_export
      file(metadata) from metadata_export
      file(branch_lengths) from branch_lengths_export
      file(nt_muts) from ancestral_nt_export
      file(aa_muts) from ancestral_aa_muts
      file(traits) from traits_export
      file(colors) from colors_export
      file(lat_long) from lat_long_export
      file(config) from config_export
    
      output:
      file "zika.json"
    
      shell:
      """
      augur export v2 \
        --tree ${tree} \
        --metadata ${metadata} \
        --node-data ${branch_lengths} \
                    ${nt_muts} \
                    ${nt_muts} \
                    ${traits} \
        --colors ${colors} \
        --lat-longs ${lat_long} \
        --auspice-config ${config} \
        --output zika.json
      """
    }
}

else {
    process export{
      publishDir "${params.outdir}", mode:'copy'
    
      input:
      file(tree) from refined_tree_export
      file(metadata) from metadata_export
      file(branch_lengths) from branch_lengths_export
      file(nt_muts) from ancestral_nt_export
      file(aa_muts) from ancestral_aa_muts
      file(colors) from colors_export
      file(lat_long) from lat_long_export
      file(config) from config_export
    
      output:
      file "auspice.json"
    
      shell:
      """
      augur export v2 \
        --tree ${tree} \
        --metadata ${metadata} \
        --node-data ${branch_lengths} \
                    ${nt_muts} \
                    ${nt_muts} \
        --colors ${colors} \
        --lat-longs ${lat_long} \
        --auspice-config ${config} \
        --output auspice.json
      """
    }
}
