# reading in data from files
counts <- read.table(
    'muscular_dystrophy.exp.tsv',
    sep = '\t',
    quote = '',
    header = TRUE,
    row.names = 'X'
)

meta <- read.table(
    'muscular_dystrophy.meta.tsv',
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