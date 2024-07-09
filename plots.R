palette <- list("Control_Fecis"="#2d9092", "MS_Fecis"="#be5e1f")

metadata <- read_tsv("~/pipeline_metagenomics/metadata.tsv")

counts <- read_tsv("~/pipeline_metagenomics/kraken2/counts/bacteria_species/normalized_counts_unique.tsv") %>% 
  gather(sample, counts, -Gene) %>% 
  mutate(log2_counts=log2(counts),
         Gene=gsub("\\[|\\]", "", Gene)) %>% 
  left_join(metadata) %>% 
  filter(!grepl("brush", condition))

bam_associated <- xlsx::read.xlsx(file="~/documents/list of BAM producers for heatmap.xlsx", sheetName="Foglio1") %>% 
  mutate(Gene=gsub("\\[|\\]", "", Gene)) %>% 
  drop_na(Gene) %>% 
  pull(Gene)

dcalca_associated <- xlsx::read.xlsx(file="~/documents/heatmap DCALCA producers.xlsx", sheetName="Foglio1", startRow = 2, header = F) %>% 
  filter(X7<0.1) %>% 
  mutate(X1=gsub("\\[|\\]", "", X1)) %>% 
  pull(X1) %>% 
  print()

## read number ----
metadata <- read_tsv("~/pipeline_metagenomics/metadata.tsv")

data <- read_tsv("~/pipeline_metagenomics/qc/multiqc_report_data/multiqc_fastqc.txt") %>% 
  mutate(sample=gsub("\\..*", "", Sample),
         read=gsub(".*\\.", "", Sample),
         "Total sequences (M)"=`Total Sequences`/1e6) %>% 
  left_join(metadata) %>% 
  filter(tissue=="Fecis")

ggplot(data) +
  aes(x = read, y = `Total sequences (M)`, fill = read) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  facet_wrap(vars(sample)) +
  theme_l +
  scale_y_continuous(labels = scales::label_comma()) +
  ggeasy::easy_remove_legend()

## heatmaps ----
counts %>% 
  filter(Gene %in% bam_associated) %>% 
  arrange(condition) %>% 
  tidyheatmaps::tidy_heatmap(rows=Gene, columns=sample, values=log2_counts, fontsize=20,
                             show_rownames=T, show_colnames=F, 
                             cluster_rows=T, cluster_cols=F,
                             annotation_col=condition, 
                             gaps_col=condition,
                             scale="row")

counts %>% 
  filter(Gene %in% dcalca_associated) %>% 
  arrange(condition) %>% 
  tidyheatmaps::tidy_heatmap(rows=Gene, columns=sample, values=log2_counts, fontsize=20,
                             show_rownames=T, show_colnames=F, 
                             cluster_rows=T, cluster_cols=F,
                             annotation_col=condition, 
                             gaps_col=condition,
                             scale="row")

## violin plots ----
counts %>% 
  filter(Gene %in% bam_associated) %>% 
  ggplot() +
  aes(x = condition, y = log2_counts, fill = condition) +
  geom_violin() +
  theme_minimal() +
  facet_wrap(vars(Gene), scales = "free_y") +
  scale_fill_manual(values = palette) +
  ggeasy::easy_remove_x_axis() +
  theme_l

counts %>% 
  filter(Gene %in% dcalca_associated) %>% 
  ggplot() +
  aes(x = condition, y = log2_counts, fill = condition) +
  geom_violin() +
  theme_minimal() +
  facet_wrap(vars(Gene), scales = "free_y") +
  scale_fill_manual(values = palette) +
  ggeasy::easy_remove_x_axis() +
  theme_l

## diversity ----
diversity <- read_tsv("~/pipeline_metagenomics/kraken2/diversity/bacteria_species/shannon.tsv") %>% 
  filter(tissue=="Fecis")

ggplot(diversity) +
  aes(x = condition, y = Shannon_index, fill = condition) +
  geom_violin() +
  geom_jitter() +
  theme_minimal() +
  scale_fill_manual(values = palette) +
  ggeasy::easy_remove_x_axis() +
  labs(y="Shannon diversity index") +
  theme_l

## multivariate ----
openxlsx::read.xlsx(xlsxFile="~/documents/multivariata_valori epatici e clinici.xlsx") %>% 
  rename(patient=1) %>% 
  select(-"patient") %>% 
  set_names(gsub("\\.", " ", colnames(.))) %>% 
  ggstatsplot::ggcorrmat(type="nonparametric", output="plot", 
                         matrix.type="upper", p.adjust.method="fdr",
                         colors=c("#2166ac", "white", "#b2182b"), 
                         ggcorrplot.args=list(tl.srt=90, pch.cex=12),
                         pch=0)
