# Data sources and assay linkage

The Inflammatory Bowel Disease Multi’omics Database (HMP2/IBDMDB) provides the stool metagenomic profiles, stool metabolomic measurements and questionnaires used in this study. The source manifest identifies each download by filename, URL and SHA-256 digest. [IBDMDB data portal](https://ibdmdb.org/results) · [Source manifest](../data/sources.csv)

| Public data product | Measurements used |
| --- | --- |
| `hmp2_metadata_2018-08-20.csv` | Participant and assay identifiers, diagnosis, stool receipt dates, dietary recency, medication and symptom reports |
| `taxonomic_profiles_3.tsv.gz` | Relative abundances of the 187 taxonomic entries ending at genus rank |
| `HMP2_metabolomics_w_metadata.biom.gz` | Processed metabolomic peak intensities and annotations, including the HILIC-negative peak annotated as butyrate (`HILn_QI36`) |

## Linking the assays

The reconstruction groups stool assay records by participant, collection identifier (`site_sub_coll`) and collection interval (`interval_days`). It includes collections from participants with Crohn’s disease or without IBD that contain both assays, consistent responses for four foods (whole grains, vegetables, fruit and red meat), and a consistent exported antibiotic flag. An identical external assay ID takes precedence when matching the two assays. Otherwise the collection must identify a unique metagenome. A collection-level match documents a shared collection record. Physical aliquot identity is unknown.

The resulting table contains 341 collections from 76 people. In 279 collections the assay IDs match exactly. In the other 62, the collection key identifies the pair. Of the linked collections, 218 come from 50 people with Crohn’s disease and 123 from 26 people without IBD. The genus matrix retains two all-zero profiles with explicit flags. The metabolomics matrix retains all 81,867 peaks and the implicit zeros of the sparse BIOM file. The reconstruction preserves source units. In the linked data, `antibiotics` stores the exported child flag. Of 306 collections with a negative flag, 19 had a blank parent medication question on both assay records. These collections remain eligible for descriptive linkage. [Aggregate reconstruction counts](../results/tables/reference_counts.csv)

## Analysis cohorts

Each analysis starts from its relevant assay and questionnaires. The within-person genus comparison uses consecutive metagenomic sampling intervals from participants with Crohn’s disease. The wider-diet analysis requires complete answers on all 21 dietary items at both endpoints. The butyrate analysis uses metabolomic visits from participants with usable measurements in both whole-grain reporting categories. Each can therefore include visits outside the dataset linked across assays. For clinical analyses, the medication screen is also required to interpret a negative antibiotic flag. [Eligibility criteria and transformations](HMP2-analysis-plan.md)

The cohort builder reads the source assays directly. It matches the 21 dietary variables by their exact metadata headers. A missing answer on another food can exclude an interval from the wider dietary analysis while retaining it for the primary comparison. The antibiotic sensitivity analysis requires a valid medication screen and a negative antibiotic flag at both endpoints. [Dietary item definitions](../data/diet_items.csv) · [Cohort builder](../R/build_cohorts.R)

The public flow table gives the observation unit, number of observations and distinct participants at each step. Rows prefixed `mgx_` refer to metagenomics, and `mbx_` to metabolomics. One person or kit can contribute to several rows. The local file `data/processed/study_cohorts.rds` contains the metagenomic interval membership and fold assignment (`mgx_intervals`), the 22 full taxonomic identifiers (`genus_panel`), and the eligible metabolomics visits (`butyrate_visits`). Raw inputs and participant-level files stay outside Git. [Aggregate cohort flow](../results/tables/cohort_flow.csv) · [Reproduction instructions](../README.md)
