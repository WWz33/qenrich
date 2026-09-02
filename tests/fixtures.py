"""Minimal fixture files for every supported format."""

import textwrap

EGGNOG = textwrap.dedent("""\
    ## eggNOG-mapper v2
    #query\tseed_ortholog\tevalue\tscore\tGOs\tKEGG_ko\tKEGG_Pathway\tPFAMs\tCOG_category
    Gene01\tK00001\t1e-50\t500\tGO:0000001,GO:0000002\tko:K00001\tko00010,ko01100\tPF00001.20\tK
    Gene02\tK00002\t1e-40\t400\tGO:0000001,GO:0000003\tko:K00002\tko00010\tPF00002.5\tE
    Gene03\tK00003\t1e-30\t300\tGO:0000002\tko:K00003\tko00020\tPF00001.30\tG
    Gene04\tK00004\t1e-20\t200\t\t\t\tPF00003\t
    Gene05\tK00005\t1e-10\t100\tGO:0000001\tko:K00005\t\t\tT
    Gene06\tK00006\t1e-10\t100\tGO:0000004\t\t\t\t
""")

IPRSCAN = textwrap.dedent("""\
    Gene01\t0123456789abcdef0123456789abcdef\t150\tPfam\tPF00001\tPkinase\t1\t200\t1.5e-30\tT\t20240101\tIPR000001\tDomain\tGO:0000001,GO:0000002\t-
    Gene02\t0223456789abcdef0123456789abcdef\t180\tPANTHER\tPTHR10000\tKinase\t5\t150\t2.2e-10\tT\t20240101\tIPR000002\tRepeat\tGO:0000001\t-
    Gene03\t0323456789abcdef0123456789abcdef\t200\tPfam\tPF00001\tPkinase\t10\t300\t1e-05\tT\t20240101\tIPR000001\tDomain\tGO:0000003\t-
""")

B2G_ANNOT = textwrap.dedent("""\
    Gene01\tGO:0000001;GO:0000002\tInterPro\tIPR000001
    Gene02\tGO:0000001\tInterPro\tIPR000002
    Gene03\tGO:0000003\tInterPro\tIPR000001
""")

B2G_TABULAR = textwrap.dedent("""\
    Sequence Name\tGO IDs\tEnzyme Codes
    Gene01\tGO:0000001, GO:0000002\t-
    Gene02\tGO:0000001\t-
    Gene03\tGO:0000003\t-
""")

PANNZER = textwrap.dedent("""\
    qpid\trank\tgo_id\tannotation\treliability
    Gene01\t1\tGO:0000001\tkinase\t0.9
    Gene01\t2\tGO:0000002\tbinding\t0.8
    Gene02\t1\tGO:0000001\tkinase\t0.7
    Gene03\t1\tGO:0000003\thydrolase\t0.6
""")

TRINOTATE = textwrap.dedent("""\
    #gene_id\ttranscript_id\tsprot_Top_BLASTX_hit\tRNAMMER\tProtId\tProtein_Length\tTop_BLASTP_hit\tPFAM\tSignalP\tTmHMM\tEGGNOG\tgene_ontology_blast\tgene_ontology_pfam\tgene_ontology_all
    Gene01\tGene01.t1\t-\t-\tp1\t300\t-\tPF00001\t-\t-\t-\tGO:0000001^molecular_function^kinase\t-\tGO:0000001^molecular_function^kinase,GO:0000002^biological_process^growth
    Gene02\tGene02.t1\t-\t-\tp2\t200\t-\tPF00002\t-\t-\t-\tGO:0000001^molecular_function^kinase\t-\tGO:0000003^cellular_component^membrane
""")

GAF = textwrap.dedent("""\
    !gaf-version: 2.2
    UniProtKB\tGene01\tG00001\t\tGO:0000001\tPMID:1\tIMP\t\tP\tkinase\t\tgene\t\tprotein\ttaxon:1\t20240101\tTAIR
    UniProtKB\tGene02\tG00002\t\tGO:0000001\tPMID:2\tIMP\t\tP\tkinase\t\tgene\t\tprotein\ttaxon:1\t20240101\tTAIR
    UniProtKB\tGene03\tG00003\t\tGO:0000002|GO:0000003\tPMID:3\tIEA\t\tP\tgrowth\t\tgene\t\tprotein\ttaxon:1\t20240101\tTAIR
    UniProtKB\tGene04\tG00004\tNOT\tGO:0000001\tPMID:4\tIEA\t\tP\tbad\t\tgene\t\tprotein\ttaxon:1\t20240101\tTAIR
""")

GFF3 = textwrap.dedent("""\
    ##gff-version 3
    chr1\tannot\tgene\t100\t500\t.\t+\t.\tID=gene:Gene01;Ontology_term=GO:0000001,GO:0000002;Dbxref=InterPro:IPR000001
    chr1\tannot\tmRNA\t100\t500\t.\t+\t.\tID=Gene01;Parent=gene:Gene01;Dbxref=InterPro:IPR000002
    chr1\tannot\tgene\t600\t900\t.\t-\t.\tID=gene:Gene02;Ontology_term=GO:0000001
    chr1\tannot\tgene\t1000\t1300\t.\t-\t.\tID=gene:Gene03;Ontology_term=GO:0000003
""")

KOFAM = textwrap.dedent("""\
    # Gene\tKO\tThreshold\tScore\tE-value
    Gene01\tK00001\t100.5\t200.3\t1e-30
    Gene02\tK00002\t90.0\t95.1\t1e-10
    Gene03\tK00001\t100.5\t150.0\t1e-05
    Gene04\t-\t-\t-\t-
""")

GENERIC = textwrap.dedent("""\
    gene\tgo_terms\tko\tipr
    Gene01\tGO:0000001,GO:0000002\tK00001\tIPR000001
    Gene02\tGO:0000001\tK00002\tIPR000002
    Gene03\tGO:0000003\t-\tIPR000001
""")

NET = textwrap.dedent("""\
    source\ttarget
    GO:0000001\tGene01
    GO:0000001\tGene02
    GO:0000001\tGene03
    GO:0000002\tGene01
    GO:0000002\tGene03
    GO:0000002\tGene04
    GO:0000003\tGene01
    GO:0000003\tGene05
    GO:0000003\tGene06
""")

GENELIST = textwrap.dedent("""\
    up\tdown
    Gene01\tGene04
    Gene02\tGene05
    Gene03\tGene06
""")

GENELIST_NOHEADER = textwrap.dedent("""\
    Gene01\tGene04
    Gene02\tGene05
    Gene03\tGene06
""")

OBO = textwrap.dedent("""\
    format-version: 1.2

    [Term]
    id: GO:0000001
    name: parent process
    namespace: biological_process
    alt_id: GO:0000099

    [Term]
    id: GO:0000002
    name: child process
    namespace: biological_process
    is_a: GO:0000001

    [Term]
    id: GO:0000003
    name: grandchild process
    namespace: cellular_component
    is_a: GO:0000002
    relationship: part_of GO:0000001

    [Term]
    id: GO:0000004
    name: obsolete thing
    namespace: biological_process
    is_a: GO:0000001
    is_obsolete: true
""")
