rule clean_assembly:
    shell:
        "rm -rf results/assembly/"

def determine_contigs(assembly, config):
    genome = lookup_in_config(config, ["assembly", assembly, "genome"])
    all_contigs = lookup_in_config(config, ["reference", genome, "genbanks"])
    return all_contigs

def assemblies(config):
    """
    Get the unique sample names from the config file
    """
    return list(config['assembly'].keys())
 

def determine_all_tx_gene_associations(config):
    '''
    Figure out which contigs are in each assembly and get the associated
    tx_gene.tsv files

    Input:
        config - the config object from Snakemake
    Output:
        a list of tx_gene.tsv association files to be created for all assemblies
    '''

    outf = []
    for assembly in assemblies(config):
        for contig in determine_contigs(assembly, config):
            outf.append("results/assembly/%s/tx_gene_association/%s_tx_gene.tsv"%(assembly, contig))
    return outf
        
        
rule run_assembly:
    input:
        expand("results/assembly/{assembly}/assembled_transcriptome.fa", assembly = assemblies(config)),
        determine_all_tx_gene_associations(config)

## GENERAL USE FOR ASSEMBLY

def contig_to_genbank(config, assembly, contig):
    genome = lookup_in_config(config, ["assembly", assembly, "genome"])
    genbank = lookup_in_config(config, ["reference", genome, "genbanks", contig])
    return genbank


def convert_contig_to_rh_fullname(contig):
    contig = convert_contig_to_rh_name(contig)
    return "'gi|xxxxxxxx|gb|%s.x|%s|'"%(contig, contig)

def convert_contig_to_rh_name(contig):
    return contig.replace(".", "_")

rule process_genbank_rh:
    message: "Processing genbank for {wildcards.assembly} for {wildcards.contig}"
    input:
        lambda wildcards: contig_to_genbank(config, wildcards.assembly, wildcards.contig)
    output:
        "results/assembly/{assembly}/process_genbank/{contig}/{contig}.{outfmt}"
    log:
        stderr="results/assembly/logs/{assembly}/process_genbank/rockhopper/{contig}_{outfmt}.err"
    threads: 1
    conda:
        "../envs/assembly.yaml"
    params:
        rh_name = lambda wildcards: convert_contig_to_rh_fullname(wildcards.contig),
        features =  lambda wildcards: lookup_in_config(config, ["assembly", wildcards.assembly, "process_genbank", "features", "value"], "CDS tRNA rRNA ncRNA"),
        qual_name = lambda wildcards: lookup_in_config(config, ["assembly", wildcards.assembly, "process_genbank", "qual_name", "value"], "gene")
    shell:
         "python3 workflow/scripts/parse_genbank.py {input} --outfmt {wildcards.outfmt} "
         "--features {params.features} "
         "--qual_name {params.qual_name} "
         "--chrm {params.rh_name}  "
         " > {output} 2> {log.stderr}"


def merge_bed_input(assembly, contig):
    assembler = lookup_in_config(config, ['assembly', assembly, 'assembler'])
    if assembler == "rockhopper":
        groupings = lookup_in_config(config, ['assembly', assembly, 'samples'])
        outpath = "results/assembly/%s/rockhopper/rockhopper_%s/%s_transcriptome.bed"
        outfiles = [outpath%(assembly,grouping,contig) for grouping in groupings]
    else:
        raise ValueError("Don't know how to merge beds for assembler %s"%(assembler))
    return outfiles

rule assembly_merge_beds:
    message: "Merging bed files for assembly {wildcards.assembly} genome {wildcards.contig}"
    input:
        lambda wildcards: merge_bed_input(wildcards.assembly, wildcards.contig)
    output:
        "results/assembly/{assembly}/merged_beds/{contig}_transcriptome.bed"
    threads: 1
    conda:
        "../envs/assembly.yaml"
    params:
        min_diff = lambda wildcards: lookup_in_config(config, ['assembly', wildcards.assembly, "merge_dist"], 100)
    log:
        stdout = "results/assembly/logs/{assembly}/merge_beds/{contig}.log",
        stderr = "results/assembly/logs/{assembly}/merge_beds/{contig}.err"
    shell:
        "python workflow/scripts/merge_beds.py {output} {input} --min_diff {params.min_diff} > {log.stdout} 2> {log.stderr}"

def find_beds(assembly, contig):
    assembler = lookup_in_config(config, ['assembly', assembly,'assembler'])
    outstr = "results/assembly/%s/%s/%s_transcriptome.bed"
    if assembler == "rockhopper":
        path= "merged_beds"
    elif assembler == "parse_genbank":
        path = "parse_genbank"
    elif assembler == "stringtie" or assembler == "trinity":
        path = "gtf_to_beds"
    elif assembler == "from_beds":
        path = lookup_in_config(config, ["assembly", assembly, "path"])
        return path + "%s_transcriptome.bed"%(contig)
    else:
        raise ValueError("Don't know where to find beds to convert for assembler"%assembler)
    return outstr%(assembly, path, contig)

def find_nonrh_fasta(assembly, contig):
    genome = lookup_in_config(config, ["assembly", assembly, "genome"])
    out = "results/alignment/process_genbank/%s/%s.fna"%(contig, contig)
    return out
    

rule assembly_bed_to_fasta:
    message: "Converting bed file into fasta file"
    input:
        to_convert=lambda wildcards: find_beds(wildcards.assembly, wildcards.contig),
        contig= lambda wildcards: find_nonrh_fasta(wildcards.assembly, wildcards.contig)
    output:
        "results/assembly/{assembly}/converted_beds/{contig}_transcriptome.fa"
    threads: 1
    conda:
        "../envs/assembly.yaml"
    log:
        stdout = "results/assembly/logs/{assembly}/bed_to_fasta/{contig}.log",
        stderr = "results/assembly/logs/{assembly}/bed_to_fasta/{contig}.err"
    shell:
        "python workflow/scripts/bed_to_fasta.py {input.to_convert} {input.contig} {output} > {log.stdout} 2> {log.stderr}"

def merge_fasta_input(assembly):
    all_contigs = determine_contigs(assembly, config)
    outpath = "results/assembly/%s/converted_beds/%s_transcriptome.fa"
    outfiles = [outpath%(assembly,contig) for contig in all_contigs]
    return outfiles

rule assembly_merge_fasta:
    message: "Combining fasta files together"
    input:
        lambda wildcards: merge_fasta_input(wildcards.assembly)
    output:
        "results/assembly/{assembly}/assembled_transcriptome.fa"
    threads: 1
    conda:
        "../envs/assembly.yaml"
    log:
        stdout = "results/assembly/logs/{assembly}/merge_fasta.log",
        stderr = "results/assembly/logs/{assembly}/merge_fasta.err"
    shell:
        "python workflow/scripts/merge_fasta.py {output} {input} > {log.stdout} 2> {log.stderr}"

## ROCKHOPPER SPECIFIC RULES

## Rockhopper input helper functions
def pull_processed_fastqs(samplename, config):
    """
    Figure out processed fastqname(s) from samplename
    """
    read1 = "results/preprocessing/trimmomatic/%s_trim_paired_R1.fastq.gz"%(samplename)
    read2 = "results/preprocessing/trimmomatic/%s_trim_paired_R2.fastq.gz"%(samplename)
    return([read1, read2])

def determine_grouping_samples(assembly, grouping, config):
    these_samples = config["assembly"][assembly]["samples"][grouping]
    # deal with a single replicate group
    if type(these_samples) is str:
        raise ValueError("Error in Assembly %s. Grouping %s. Grouping must have more than one sample."%(assembly,grouping))
    return these_samples

def determine_fastqs_bygroup(assembly, grouping, config):
    these_samples = determine_grouping_samples(assembly, grouping, config)
    outfiles = []
    for sample in these_samples:
        outfiles.extend(pull_processed_fastqs(sample, config))
    return outfiles

def configure_contigs_for_rockhopper(assembly, grouping, config):
    all_contigs = determine_contigs(assembly, config)
    outfiles = []
    outfiles.extend(["results/assembly/%s/process_genbank/%s/%s.rnt"%(assembly, contig,contig) for contig in all_contigs])
    outfiles.extend(["results/assembly/%s/process_genbank/%s/%s.ptt"%(assembly, contig,contig) for contig in all_contigs])
    outfiles.extend(["results/assembly/%s/process_genbank/%s/%s.fna"%(assembly, contig,contig) for contig in all_contigs])
    return outfiles

# rockhopper output helper functions

def determine_output_from_rockhopper(assembly, grouping, config):
    all_contigs = determine_contigs(assembly, config)
    outfiles = []
    for contig in all_contigs:
        outfiles.extend(["results/assembly/%s/rockhopper_%s/%s_transcripts.txt"%(assembly,grouping,contig) for contig in all_contigs])
        outfiles.extend(["results/assembly/%s/rockhopper_%s/%s_operons.txt"%(assembly,grouping,contig) for contig in all_contigs])
    return outfiles

# rockhopper command line parameter helper functions

def configure_contig_input_for_rockhopper(assembly,grouping,config):
    contigs = determine_contigs(assembly, config)
    return ",".join(["results/assembly/%s/process_genbank/%s/"%(assembly, contig) for contig in contigs])

def configure_librarytype_for_rockhopper(assembly,grouping,config, pep):
    these_samples = determine_grouping_samples(assembly, grouping,config)
    # check that samples all have the same library type
    library_types = set()
    for sample in these_samples:
        library_types.add(lookup_sample_metadata(sample, "library_type", pep))
    if len(library_types) > 1:
        raise ValueError("Samples must have same library type for rockhopper. Consider a smaller grouping based on library")
    library = library_types.pop()
    if library == "rf-stranded":
        out = "-rf"
    elif library == "fr-stranded":
        out = "-fr"
    elif library == "unstranded":
        out = "-s false" 
    else:
        raise ValueError("Only unstranded, rf-stranded, and fr-stranded are supported. Not %s"%(library))
    return out

def configure_sample_input_for_rockhopper(assembly, grouping, config):
    """
    Convert all samples to a list for rockhopper. Groupings are biological
    replicates that are treated as independent samples here.
    """
    these_samples = determine_grouping_samples(assembly, grouping, config)
    outstr = ""
    for sample in these_samples:
        this_set = "%".join(pull_processed_fastqs(sample, config))
        outstr = outstr + " " + this_set
    return outstr

def configure_labels_for_rockhopper(assembly, grouping, config):
    """
    Convert all samples to a list for rockhopper. Treating each replicate as an
    individual condition
    """
    these_samples = determine_grouping_samples(assembly, grouping, config)
    outstr = ",".join(these_samples)
    return outstr

rule rockhopper_assembly:
    message: "Running rockhopper for assembly {wildcards.assembly}, grouping {wildcards.grouping}"
    input:
        lambda wildcards: determine_fastqs_bygroup(wildcards.assembly,
                                                    wildcards.grouping,
                                                    config),
        lambda wildcards: configure_contigs_for_rockhopper(wildcards.assembly,
                                                        wildcards.grouping,
                                                        config)
    output:
        touch("results/assembly/{assembly}/rockhopper/rockhopper_{grouping}.done")
    params:
        path_to_rockhopper = "~/src/",
        contigs = lambda wildcards: configure_contig_input_for_rockhopper(wildcards.assembly, wildcards.grouping, config),
        files = lambda wildcards: configure_sample_input_for_rockhopper(wildcards.assembly, wildcards.grouping, config),
        labels = lambda wildcards: configure_labels_for_rockhopper(wildcards.assembly, wildcards.grouping, config),
        library = lambda wildcards: configure_librarytype_for_rockhopper(wildcards.assembly, wildcards.grouping, config, pep),
        outdir = "results/assembly/{assembly}/rockhopper/rockhopper_{grouping}"
    threads: 10
    resources:
        mem_mb=5000
    conda:
        "../envs/assembly.yaml"
    log:
        stdout="results/assembly/logs/{assembly}/rockhopper/rockhopper_{grouping}.log",
        stderr="results/assembly/logs/{assembly}/rockhopper/rockhopper_{grouping}.err"
    shell:
        "java -Xmx{resources.mem_mb}m -cp {params.path_to_rockhopper}Rockhopper.jar "
        "Rockhopper {params.files} -p {threads} {params.library} -L {params.labels} "
        "-g {params.contigs} "
        "-o {params.outdir} > {log.stdout} 2> {log.stderr}"

rule parse_rockhopper:
    message: "Parsing rockhopper output for assembly {wildcards.assembly}, grouping {wildcards.grouping}, chrm {wildcards.chrm}"
    input:
        "results/assembly/{assembly}/rockhopper/rockhopper_{grouping}.done"
    output:
        "results/assembly/{assembly}/rockhopper/rockhopper_{grouping}/{chrm}_transcriptome.bed"
    params:
        operons = lambda wildcards: "results/assembly/%s/rockhopper/rockhopper_%s/%s_operons.txt"%(wildcards.assembly,wildcards.grouping, convert_contig_to_rh_name(wildcards.chrm)),
        transcripts = lambda wildcards: "results/assembly/%s/rockhopper/rockhopper_%s/%s_transcripts.txt"%(wildcards.assembly,wildcards.grouping, convert_contig_to_rh_name(wildcards.chrm))
    conda:
        "../envs/assembly.yaml"
    threads: 1
    log:
        stderr="results/assembly/logs/{assembly}/rockhopper/rockhopper_{grouping}/{chrm}_transcriptome.err"
    shell:
        "python3 workflow/scripts/rh_tools.py {params.operons} {params.transcripts} "
        "{wildcards.chrm} --min_size 30 --skip_unstranded > {output} 2> {log.stderr}"


def tx_gene_association_input(assembly, contig):
    bed = find_beds(assembly, contig)
    ptt = "results/assembly/{assembly}/process_genbank/{contig}/{contig}.ptt"
    rnt = "results/assembly/{assembly}/process_genbank/{contig}/{contig}.rnt"
    outfiles = [bed, ptt, rnt]
    return outfiles

rule tx_gene_association:
    message: "Associating transcripts with genes for assembly {wildcards.assembly}, chrm {wildcards.contig}"
    input:
        lambda wildcards: tx_gene_association_input(wildcards.assembly, wildcards.contig)
    output:
        "results/assembly/{assembly}/tx_gene_association/{contig}_tx_gene.tsv"
    log:
        stdout="results/assembly/logs/{assembly}/tx_gene_association/{contig}_tx_gene.log",
        stderr="results/assembly/logs/{assembly}/tx_gene_association/{contig}_tx_gene.err"
    conda:
        "../envs/assembly.yaml"
    threads: 1
    shell:
        "Rscript workflow/scripts/associate_gene_transcripts.R {input} "
        "{output} > {log.stdout} 2> {log.stderr}"



################### WIP BELOW HERE

## STRINGTIE specific rules


def get_assembly_fastas(config, assembly):
    out_fastas = []
    for contig in lookup_in_config(config, ["assembly", assembly, "contigs"]):
        this_fasta = lookup_in_config(config, \
        ["contigs", contig, "fasta"], \
        "results/assembly/%s/process_genbank/%s/%s.fna"%(assembly,contig, contig))
        out_fastas.append(this_fasta) 
    return out_fastas

def masked_regions_file_for_combine_fastas(config, assembly):
    outfile = lookup_in_config(config, ["assembly", assembly, "masked_regions"], "none")
    if outfile is not "none":
        out = "--masked_regions %s"%(outfile)
    else:
        out = ""
    return out

rule assembly_combine_fastas:
    message: "Generating genome fasta for assembly {wildcards.assembly}"
    input:
        lambda wildcards: get_assembly_fastas(config, wildcards.assembly)
    output:
        "results/assembly/{assembly}/combine_fasta/{assembly}_genome.fa",
        "results/assembly/{assembly}/combine_fasta/{assembly}_genome_contig_sizes.tsv",
        "results/assembly/{assembly}/combine_fasta/{assembly}_genome_mappable_size.txt"
    log:
        stdout="results/assembly/{assembly}/logs/combine_fasta/{assembly}_genome.log",
        stderr="results/assembly/{assembly}/logs/combine_fasta/{assembly}_genome.err"
    threads: 1
    params:
        masked_regions = lambda wildcards: masked_regions_file_for_combine_fastas(config, wildcards.assembly)
    conda:
        "../envs/assembly.yaml"
    shell:
         "python3 workflow/scripts/combine_fasta.py "
         "results/assembly/{wildcards.assembly}/combine_fasta/{wildcards.assembly}_genome "
         "{params.masked_regions} "
         "{input} > {log.stdout} 2> {log.stderr}"

rule assembly_bowtie2_index:
    input:
        "results/assembly/{assembly}/combine_fasta/{assembly}_genome.fa"
    output:
        "results/assembly/{assembly}/bowtie2_index/{assembly}_genome.1.bt2"
    threads:
        5
    log:
        stdout="results/assembly/{assembly}/logs/bowtie2_index/{assembly}_genome.log",
        stderr="results/assembly/{assembly}/logs/bowtie2_index/{assembly}_genome.err" 
    conda:
        "../envs/stringtie.yaml"
    shell:
        "bowtie2-build --threads {threads} "
        "{input} "
        "results/assembly/{wildcards.assembly}/bowtie2_index/{wildcards.assembly}_genome "
        "> {log.stdout} 2> {log.stderr}"

rule assembly_bowtie2_map:
    input:
        in1="results/preprocessing/trimmomatic/{sample}_trim_paired_R1.fastq.gz",
        in2="results/preprocessing/trimmomatic/{sample}_trim_paired_R2.fastq.gz",
        bt2_index_file= "results/assembly/{assembly}/bowtie2_index/{assembly}_genome.1.bt2"
    output:
        temp("results/assembly/{assembly}/bowtie2/{sample}.bam")
    log:
        stderr="results/assembly/{assembly}/logs/bowtie2/{sample}_bt2.log" 
    params:
        bt2_index= "results/assembly/{assembly}/bowtie2_index/{assembly}_genome",
        bowtie2_param_string= lambda wildcards: lookup_in_config_persample(config,\
        pep, ["assembly", "bowtie2_map", "bowtie2_param_string"], wildcards.sample,\
        "--end-to-end --very-sensitive --phred33"),
        samtools_view_param_string = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["assembly", "bowtie2_map", "samtools_view_param_string"],\
        wildcards.sample, "-b")
    threads: 
        5
    conda:
        "../envs/stringtie.yaml"
    shell:
        "bowtie2 -x {params.bt2_index} -p {threads} "
        "-1 {input.in1} -2 {input.in2} --phred33  "
        "{params.bowtie2_param_string} 2> {log.stderr} "
        "| samtools view {params.samtools_view_param_string} > {output}"

rule assembly_bam_sort:
    input:
        "results/assembly/{assembly}/bowtie2/{sample}.bam"
    output:
        "results/assembly/{assembly}/bowtie2_srt/{sample}_sorted.bam"
    log:
        stderr="results/assembly/{assembly}/logs/bowtie2/{sample}_bt2_sort.log"
    conda:
        "../envs/stringtie.yaml"
    shell:
        "samtools sort {input} > {output} 2> {log.stderr}"

rule assembly_bam_index:
    input:
        "results/assembly/{assembly}/bowtie2_srt/{sample}_sorted.bam"
    output:
        "results/assembly/{assembly}/bowtie2_srt/{sample}_sorted.bam.bai"
    log:
        stdout="results/assembly/{assembly}/logs/bowtie2/{sample}_bt2_index.log",
        stderr="results/assembly/{assembly}/logs/bowtie2/{sample}_bt2_index.log"
    conda:
        "../envs/stringtie.yaml"
    shell:
        "samtools index {input} {output} > {log.stdout} 2> {log.stderr}"

def configure_librarytype_for_stringtie(sample, pep):
    library_type = lookup_sample_metadata(sample, "library_type", pep)
    if library_type == "rf-stranded":
        out = "--rf"
    elif library_type == "fr-stranded":
        out = "--fr"
    elif library_type == "unstranded":
        out = ""
    else: 
        raise ValueError("Only unstranded, rf-stranded, and fr-stranded are supported. Not %s"%(library_type))
    return out


def get_assembly_gffs(assembly, config):
    out_gffs = []
    for contig in lookup_in_config(config, ["assembly", assembly, "contigs"]):
        this_genbank = lookup_in_config(config, \
        ["contigs", contig, "genbank"], "none")
        if this_genbank is not "none":
            out_gffs.append("results/assembly/process_genbank/%s/%s"%(contig,contig) + ".gff")
    return out_gffs

rule combine_ref_gffs:
    input:
        gffs = lambda wildcards: get_assembly_gffs(wildcards.assembly, config)
    output:
        "results/assembly/{assembly}/combined_ref_gff/all_genes.gff"
    shell:
        "for FILE in {input.gffs}; do cat ${{FILE}} "
        " >> results/assembly/{wildcards.assembly}/combined_ref_gff/all_genes.gff; done"

rule stringtie:
    input:
        bam = "results/assembly/{assembly}/bowtie2_srt/{sample}_sorted.bam",
        reference = lambda wildcards: lookup_in_config(config, ['assembly', wildcards.assembly, 'reference_gff'], "results/assembly/%s/combined_ref_gff/all_genes.gff"%(wildcards.assembly))
    output:
        gtf = "results/assembly/{assembly}/stringtie/{sample}/{sample}.gtf"
    log:
        stdout = "results/assembly/{assembly}/logs/stringtie/{sample}.log",
        stderr = "results/assembly/{assembly}/logs/stringtie/{sample}.err"
    conda:
        "../envs/stringtie.yaml"
    threads:
        5
    params:
        library = lambda wildcards: configure_librarytype_for_stringtie(wildcards.sample, pep),
        stringtie_param_string = lambda wildcards: 
            lookup_in_config_persample(config, pep, 
            ["assembly", wildcards.assembly, "stringtie_param_string"], 
            wildcards.sample, default = "-m 30")
    shell:
        "stringtie {input.bam} -o {output.gtf} -p {threads} {params.library} " 
        "{params.stringtie_param_string} "
        #"-G {input.reference} "
        "> {log.stdout} 2> {log.stderr}" 

def get_assembly_gtfs(assembly, config):
    out_gtfs = []
    assembler = lookup_in_config(config, ["assembly", assembly, "assembler"])
    if assembler == "stringtie":
        outpath = "results/assembly/%s/stringtie/%s/%s.gtf"
    elif assembler == "trinity":
        outpath = "results/assembly/%s/gmap/%s/%s.gtf"
    else:
        raise ValueError("Don't know how to find gtfs for assembler %s"%assembler)
    for sample in lookup_in_config(config, ["assembly", assembly, "samples"]):
        out_gtfs.append(outpath%(assembly, sample, sample))
    return out_gtfs
        

rule stringtie_merge:
    input:
        gtfs = lambda wildcards: get_assembly_gtfs(wildcards.assembly, config)
    output:
        gtf = "results/assembly/{assembly}/stringtie_merge/all_transcripts.gtf"
    log:
        stdout = "results/assembly/{assembly}/logs/stringtie_merge/stringtie_merge.log",
        stderr = "results/assembly/{assembly}/logs/stringtie_merge/stringtie_merge.err"
    conda:
        "../envs/stringtie.yaml"
    shell:
        "stringtie --merge -o {output.gtf} {input.gtfs} > {log.stdout} "
        "2> {log.stderr}"

def feature_type_to_filter(config, assembly):
    assembler = lookup_in_config(config, ["assembly", assembly, "assembler"])
    if assembler == "stringtie" or assembler == "trinity":
        out = "transcript"
    else:
        raise ValueError("Don't know feature type for assembler %s"%(assembler))
    return out

rule gtf_to_beds:
    input:
        "results/assembly/{assembly}/stringtie_merge/all_transcripts.gtf"
    output:
        "results/assembly/{assembly}/gtf_to_beds/{contig}_transcriptome.bed"
    log:
        stdout = "results/assembly/{assembly}/logs/gtf_to_beds/{contig}.log",
        stderr = "results/assembly/{assembly}/logs/gtf_to_beds/{contig}.err"
    conda:
        "../envs/assembly.yaml"
    params:
        ftype = lambda wildcards: feature_type_to_filter(config, wildcards.assembly)
    shell:
        "Rscript workflow/scripts/convert_gtf_to_beds.R {input} {output} "
        "{params.ftype} {wildcards.contig}  > {log.stdout} 2> {log.stderr} "

## RULES FOR Trinity De Novo assembler

def configure_librarytype_for_trinity(sample, pep):
    library_type = lookup_sample_metadata(sample, "library_type", pep)
    if library_type == "rf-stranded":
        out = "--SS_lib_type RF"
    elif library_type == "fr-stranded":
        out = "--SS_lib_type FR"
    elif library_type == "unstranded":
        out = ""
    else: 
        raise ValueError("Only unstranded, rf-stranded, and fr-stranded are supported. Not %s"%(library_type))
    return out

rule trinity:
    input:
        bam = "results/assembly/{assembly}/bowtie2_srt/{sample}_sorted.bam"
    output:
        "results/assembly/{assembly}/trinity/trinity_{sample}/Trinity-GG.fasta"
    log:
        stdout = "results/assembly/{assembly}/logs/trinity/{sample}.log",
        stderr = "results/assembly/{assembly}/logs/trinity/{sample}.err"
    conda:
        "../envs/trinity.yaml"
    threads: 10
    resources:
        mem_mb=10000
    params: 
        library = lambda wildcards: configure_librarytype_for_trinity(wildcards.sample, pep)
    shell:
        "Trinity --genome_guided_bam {input.bam} "
        "--genome_guided_max_intron 10000 --max_memory 10G --CPU {threads} "
        "{params.library} "
        "--output 'results/assembly/{wildcards.assembly}/trinity/trinity_{wildcards.sample}' "
        "--jaccard_clip "
        "--full_cleanup > {log.stdout} 2> {log.stderr}"

rule gmap_build:     
    input:
        lambda wildcards: get_assembly_fastas(config, wildcards.assembly)
    output:
        "results/assembly/{assembly}/gmap_build/{assembly}_genome/{assembly}_genome.chromosome"
    log:
        stdout = "results/assembly/{assembly}/logs/gmap_build/{assembly}_genome.log",
        stderr = "results/assembly/{assembly}/logs/gmap_build/{assembly}_genome.err"
    conda:
        "../envs/trinity.yaml"
    shell:
         "gmap_build -d {wildcards.assembly}_genome {input} " 
         "-D 'results/assembly/{wildcards.assembly}/gmap_build' "
         "> {log.stdout} 2> {log.stderr} "

rule gmap:
    input:
        in_idx="results/assembly/{assembly}/gmap_build/{assembly}_genome/{assembly}_genome.chromosome",
        in_trinity ="results/assembly/{assembly}/trinity/trinity_{sample}/Trinity-GG.fasta"
    output:
        "results/assembly/{assembly}/gmap/{sample}/{sample}.gtf"
    log:
        stderr = "results/assembly/{assembly}/logs/gmap/{sample}.err"
    params:
        idx_dir = "results/assembly/{assembly}/gmap_build/"
    conda:
        "../envs/trinity.yaml"
    shell:
        "gmap -D {params.idx_dir} -d {wildcards.assembly}_genome {input.in_trinity} -f "
        "gff3_gene | sed '3,${{/^#/d}}' > {output} 2> {log.stderr}"
