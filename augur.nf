#!/usr/bin/env nextflow

//Description: Nextflow implementation of Nextstrain's Augur pipeline
//Available: https://github.com/nextstrain/augur
//Authors of this Nextflow: Abigail Shockey
//Email: abigail.shockey@slh.wisc.edu

// Input channels

Channel
    .fromPath( "${params.reference}")
    .ifEmpty { exit 1, "Cannot find reference sequence in: ${params.reference}" }
    .into { reference_alignment; reference_translate }

Channel
    .fromPath( "${params.metadata}")
    .ifEmpty { exit 1, "Cannot find metadata file in: ${params.metadata}" }
    .into { filter_metadata; refine_tree_metadata; traits_metadata; metadata_export; metadata_export_traits }

Channel
    .fromPath( "${params.colors}")
    .ifEmpty { exit 1, "Cannot find colors file in: ${params.colors}" }
    .into { colors_export; colors_export_traits }

Channel
    .fromPath( "${params.lat_long}")
    .ifEmpty { exit 1, "Cannot find latitude and longitude file in: ${params.lat_long}" }
    .into { lat_long_export; lat_long_export_traits }

Channel
    .fromPath( "${params.config}")
    .ifEmpty { exit 1, "Cannot find Auspice config file in: ${params.auspice_config}" }
    .into { config_export; config_export_traits }

// Filter sequences if sequences to drop are included
if (params.filter) {
    Channel
        .fromPath( "${params.sequences}")
        .ifEmpty { exit 1, "Cannot find input sequences in: ${params.sequences}" }
        .set { sequences }
    Channel
        .fromPath( "${params.filter}")
        .ifEmpty { exit 1, "Cannot find filter file in: ${params.filter}" }
        .set { sequences }

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
        --min-date ${params.min_date}
      """
    }
}

else {
    Channel
        .fromPath( "${params.sequences}")
        .ifEmpty { exit 1, "Cannot find input sequences in: ${params.sequences}" }
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
  file "refined_tree.newick" into refined_tree_traits, refined_tree_ancestral, refined_tree_translate, refined_tree_export, refined_tree_export_traits
  file "branch_lengths.json" into branch_lengths_export, branch_lengths_export_traits

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
  file "nt_muts.json" into ancestral_nt_translate, ancestral_nt_export, ancestral_nt_export_traits

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
  file "aa_muts.json" into ancestral_aa_export, ancestral_aa_export_traits

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
  params.traits == "TRUE"

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


process export{
  publishDir "${params.outdir}", mode:'copy'

  input:
  file(tree) from refined_tree_export
  file(metadata) from metadata_export
  file(branch_lengths) from branch_lengths_export
  file(nt_muts) from ancestral_nt_export
  file(aa_muts) from ancestral_aa_export
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
                ${aa_muts} \
    --colors ${colors} \
    --lat-longs ${lat_long} \
    --auspice-config ${config} \
    --output auspice.json
  """
}

if (params.traits) {
    process export_with_traits{
        publishDir "${params.outdir}", mode:'copy'

        input:
        file(tree) from refined_tree_export_traits
        file(metadata) from metadata_export_traits
        file(branch_lengths) from branch_lengths_export_traits
        file(nt_muts) from ancestral_nt_export_traits
        file(aa_muts) from ancestral_aa_export_traits
        file(traits) from traits_export
        file(colors) from colors_export_traits
        file(lat_long) from lat_long_export_traits

        output:
        file "auspice.json"
        
        shell:
        """
        augur export v2 \
        --tree ${tree} \
        --metadata ${metadata} \
        --node-data ${branch_lengths} \
                    ${nt_muts} \
                    ${aa_muts} \
                    ${traits} \
        --colors ${colors} \
        --lat-longs ${lat_long} \
        --auspice-config ${config} \
        --output auspice.json
        """
    }
}
