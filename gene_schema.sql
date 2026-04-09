-- =============================================================================
-- Baylor Genetics Technical Assessment - Part 2
-- SQL schema + seed data + queries based on gene_output.csv
-- =============================================================================


-- -----------------------------------------------------------------------------
-- SCHEMA
-- -----------------------------------------------------------------------------

-- Core gene record (one row per HGNC entry)
CREATE TABLE gene (
    id            INTEGER      PRIMARY KEY,
    hgnc_id       TEXT         NOT NULL UNIQUE,   -- e.g. HGNC:13394
    gene_name     TEXT         NOT NULL,           -- HGNC approved name
    hg38_coords   TEXT,                            -- chr:start-end (GRCh38)
    hg19_coords   TEXT                             -- chr:start-end (GRCh37/hg19)
);

-- One row per alias (alias_symbol + prev_symbol from HGNC)
CREATE TABLE gene_alias (
    id       INTEGER PRIMARY KEY,
    gene_id  INTEGER NOT NULL REFERENCES gene(id) ON DELETE CASCADE,
    alias    TEXT    NOT NULL
);

-- Normalised disease vocabulary (one row per unique disease name)
CREATE TABLE disease (
    id           INTEGER PRIMARY KEY,
    disease_name TEXT    NOT NULL UNIQUE
);

-- Many-to-many bridge: a gene can link to multiple diseases and vice-versa
CREATE TABLE gene_disease (
    gene_id    INTEGER NOT NULL REFERENCES gene(id)    ON DELETE CASCADE,
    disease_id INTEGER NOT NULL REFERENCES disease(id) ON DELETE CASCADE,
    PRIMARY KEY (gene_id, disease_id)
);


-- -----------------------------------------------------------------------------
-- SEED DATA  (sourced from gene_output.csv / gene_output.json)
-- -----------------------------------------------------------------------------

INSERT INTO gene (id, hgnc_id, gene_name, hg38_coords, hg19_coords) VALUES
    (1, 'HGNC:7060',  'matrix Gla protein',
        'chr12:14880864-14887639',  'chr12:15034115-15038860'),
    (2, 'HGNC:13394', 'NPHS2 stomatin family member, podocin',
        'chr1:179550494-179575952', 'chr1:179519674-179545087'),
    (3, 'HGNC:2204',  'collagen type IV alpha 3 chain',
        'chr2:227164624-227314792', 'chr2:228029281-228179508'),
    (4, 'HGNC:19903', 'Ras related GTP binding D',
        'chr6:89364616-89412735',   'chr6:90074355-90121989'),
    (5, 'HGNC:11621', 'HNF1 homeobox A',
        'chr12:120978543-121002512','chr12:121416346-121440315'),
    (6, 'HGNC:618',   'apolipoprotein L1',
        'chr22:36253071-36267530',  'chr22:36649056-36663576'),
    (7, 'HGNC:644',   'androgen receptor',
        'chrX:67543353-67730619',   'chrX:66764465-66950461'),
    (8, 'HGNC:4868',  'HECT and RLD domain containing E3 ubiquitin protein ligase 2',
        'chr15:28111040-28322179',  'chr15:28356186-28567298');

-- Aliases (pipe-separated values from HGNC alias_symbol + prev_symbol)
INSERT INTO gene_alias (gene_id, alias) VALUES
    -- HGNC:7060 MGP  — no aliases in HGNC
    -- HGNC:13394 NPHS2
    (2, 'SRN1'),
    (2, 'PDCN'),
    -- HGNC:2204 COL4A3 — no aliases in HGNC
    -- HGNC:19903 RRAGD
    (4, 'DKFZP761H171'),
    (4, 'bA11D8.2.1'),
    -- HGNC:11621 HNF1A
    (5, 'HNF1'),
    (5, 'LFB1'),
    (5, 'HNF1α'),
    (5, 'MODY3'),
    (5, 'TCF1'),
    -- HGNC:618 APOL1
    (6, 'APOL'),
    -- HGNC:644 AR
    (7, 'AIS'),
    (7, 'NR3C4'),
    (7, 'SMAX1'),
    (7, 'HUMARA'),
    (7, 'DHTR'),
    (7, 'SBMA'),
    -- HGNC:4868 HERC2
    (8, 'jdf2'),
    (8, 'p528'),
    (8, 'D15F37S1');

-- Disease vocabulary
INSERT INTO disease (id, disease_name) VALUES
    ( 1, 'total joint arthroplasty'),
    ( 2, 'genetic disorder'),
    ( 3, 'Short palm'),
    ( 4, 'nephrotic syndrome'),
    ( 5, 'nephrotic syndrome, type 2'),
    ( 6, 'familial idiopathic steroid-resistant nephrotic syndrome'),
    ( 7, 'focal segmental glomerulosclerosis'),
    ( 8, 'Nephrotic range proteinuria'),
    ( 9, 'congenital nephrotic syndrome, Finnish type'),
    (10, 'Proteinuria'),
    (11, 'autosomal dominant Alport syndrome'),
    (12, 'Alport syndrome 3b, autosomal recessive'),
    (13, 'autosomal recessive Alport syndrome'),
    (14, 'Alport syndrome'),
    (15, 'hypomagnesemia 7, renal, with or without dilated cardiomyopathy'),
    (16, 'cardiomyopathy'),
    (17, 'nephrocalcinosis'),
    (18, 'heart disease'),
    (19, 'heart failure'),
    (20, 'Abnormal nasolacrimal system morphology'),
    (21, 'MODY'),
    (22, 'type 2 diabetes mellitus'),
    (23, 'diabetes mellitus'),
    (24, 'maturity-onset diabetes of the young type 3'),
    (25, 'type 1 diabetes mellitus'),
    (26, 'monogenic diabetes'),
    (27, 'maturity-onset diabetes of the young'),
    (28, 'chronic kidney disease'),
    (29, 'kidney disease'),
    (30, 'kidney failure'),
    (31, 'partial androgen insensitivity syndrome'),
    (32, 'developmental delay with autism spectrum disorder and gait instability'),
    (33, 'eye color');

-- Gene–disease associations
INSERT INTO gene_disease (gene_id, disease_id) VALUES
    -- MGP (1)
    (1,  1), (1,  2), (1,  3),
    -- NPHS2 (2)
    (2,  4), (2,  5), (2,  6), (2,  7), (2,  2), (2,  8), (2,  9), (2, 10),
    -- COL4A3 (3)
    (3, 11), (3, 12), (3, 13), (3, 14),
    -- RRAGD (4)
    (4, 15), (4, 16), (4, 17), (4, 18), (4, 19), (4, 20), (4,  2),
    -- HNF1A (5)
    (5, 21), (5, 22), (5, 23), (5, 24), (5, 25), (5, 26), (5, 27),
    -- APOL1 (6)
    (6,  7), (6, 28), (6, 29), (6, 10), (6, 30), (6,  8),
    -- AR (7)
    (7, 31),
    -- HERC2 (8)
    (8, 32), (8, 33);


-- =============================================================================
-- QUERIES
-- =============================================================================

-- -----------------------------------------------------------------------------
-- Query 1: HGNC ID and disease connection
-- Each row = one gene–disease pair, ordered by HGNC ID then disease name.
-- -----------------------------------------------------------------------------

SELECT
    g.hgnc_id,
    d.disease_name
FROM gene          g
JOIN gene_disease  gd ON gd.gene_id    = g.id
JOIN disease       d  ON d.id          = gd.disease_id
ORDER BY
    g.hgnc_id,
    d.disease_name;


-- -----------------------------------------------------------------------------
-- Query 2: HGNC Gene Name and all aliases (collapsed to one row per gene)
-- Uses GROUP_CONCAT (SQLite/MySQL) — replace with STRING_AGG for PostgreSQL.
-- -----------------------------------------------------------------------------

SELECT
    g.hgnc_id,
    g.gene_name                              AS hgnc_gene_name,
    GROUP_CONCAT(ga.alias, ' | ')            AS gene_aliases
FROM gene       g
LEFT JOIN gene_alias ga ON ga.gene_id = g.id
GROUP BY
    g.id,
    g.hgnc_id,
    g.gene_name
ORDER BY
    g.hgnc_id;

-- PostgreSQL variant (comment out the query above and use this instead):
-- SELECT
--     g.hgnc_id,
--     g.gene_name                                          AS hgnc_gene_name,
--     STRING_AGG(ga.alias, ' | ' ORDER BY ga.alias)       AS gene_aliases
-- FROM gene       g
-- LEFT JOIN gene_alias ga ON ga.gene_id = g.id
-- GROUP BY g.id, g.hgnc_id, g.gene_name
-- ORDER BY g.hgnc_id;
