# reading in data from files
counts <- read.table(
    'psoriasis.all.exp.tsv',
    sep = '\t',
    quote = '',
    header = TRUE,
    row.names = 'X'
)

meta <- read.table(
    'psoriasis.all.meta.tsv',
    sep = '\t',
    quote = '',
    header = TRUE,
    row.names = 'X'
)

# applying recombatseq correction using disease as wanted covariate
corrected <- reComBat.seq(
    t(counts), 
    batch=meta$sra_study_acc, 
    wanted.variation=meta['Disease']
)
