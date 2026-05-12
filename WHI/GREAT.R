library(dplyr)
library(tidyverse)

###### Prepare input for Genomic Regions Enrichment of Annotations Tool (GREAT) #######
# Get 841 CpGs data
coef <- read.csv("~/Desktop/GREAT_input/meta_coef_sig_annotated_SLEs.csv")

coef$chr <- paste0("chr", coef$CHR)

df <- data.frame(chr = coef$chr, start = coef$MAPINFO, end = coef$MAPINFO)
rownames(df) <- coef$IlmnID

glimpse(df)

bed <- data.frame(
  chr   = df$chr,
  start = df$start - 1,  # BED is 0-based
  end   = df$end         # end is exclusive
)

head(bed)

write.table(
  bed,
  file = "~/Desktop/GREAT_input/WHI_MRS_CpGs.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# Data input into: http://great.stanford.edu/public/html/

###### Plot output from GREAT #######
install.packages("roperators")
library(roperators)

file = "~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Zannas_lab/MRS TOPMed/GREAT_output/greatExportAll.tsv"


total_lines <- length(readLines(file))
rows_to_read <- total_lines - 8

tsv <- read.table(
  file,
  sep = "\t",
  skip = 3,
  nrows = rows_to_read,
  header = TRUE,
  quote = "",
  fill = TRUE,
  comment.char = ""
)


tsv[1:5,2]
tsv[496:500,2]

# create figure
library(ggplot2)
file = "~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Zannas_lab/MRS TOPMed/GREAT_output/shown-GOBiologicalProcess.tsv"

# Read header only
header <- readLines(file, n = 2)[2]
cn <- strsplit(header, "\t")[[1]]
cn <- trimws(cn)


tsv1 <- read.table(
  file,
  sep = "\t",
  skip = 2,
  header = FALSE,
  quote = "",
  fill = TRUE,
  comment.char = ""
)

cn
colnames(tsv)
cn1 <- c("Desc", "BinomRank", "BinomP", "BinomFdrQ", "RegionFoldEnrich", 
         "ObsRegions", "SetCov", "HyperRank", "HyperFdrQ", 
         "GeneFoldEnrich","ObsGenes", "TotalGenes", "GeneSetCov")
tsv1[,14] <- NULL
colnames(tsv1) <- cn1

glimpse(tsv1)
tsv1 <- tsv1 %>% mutate(logBinomP = -log10(BinomFdrQ))
tsv1 <- tsv1 %>% mutate(logRFC = log2(RegionFoldEnrich))

tsv_arranged <- tsv1 %>% arrange(desc(logBinomP))
identical(tsv1$Desc, tsv_arranged$Desc)
# the output is already in descending order
o <- tsv1$Desc
rev0 <- rev(o)

tsv1$Desc <- factor(tsv1$Desc, levels = rev0)

# create bar plot
library(ggplot2)

p = ggplot(data = tsv1, aes(x = Desc, y = logRFC, fill = logBinomP)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(logRFC, 2)), hjust = -0.1, color = "black", size = 3.5) +
  coord_flip() + scale_fill_gradient(low = "blue", high = "red", name = "-log10(FDR p-value)") + 
  xlab("") + ylab("log2(Region Fold Enrichment)") + 
  ggtitle(label = "GO Biological Processes") + 
  theme(axis.text.x = element_text(color = "black", size = 10), 
        axis.text.y = element_text(color = "black", size = 10),
        plot.title = element_text(color = "black", face = "bold")) +
  scale_y_continuous(expand = expansion(add = c(0, 1))) + 
  scale_x_discrete(expand=c(0,0.01)) +
  theme_bw()
  
ggsave("~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Zannas_lab/MRS TOPMed/GREAT_output/GREAT_barplot.pdf", 
       p, width = 10, height = 5)

# create dot plot

p2 <- ggplot(tsv1, aes(
  x = logRFC,
  y = Desc,
  size = ObsGenes,
  color = logBinomP
)) +
  geom_point() +
  geom_segment(aes(x = 0, xend = logRFC, y = Desc, yend = Desc, color = logBinomP), 
               linewidth = 1, show.legend = FALSE) +
  scale_color_gradient(
    low = "blue",
    high = "red",
    name = "-log10(FDR p-value)"
  ) +
  scale_size(
    range = c(3, 8),
    name = "Observed Genes"
  ) +
  labs(
    x = "log2(Region Fold Enrichment)",
    y = NULL
  ) +
  theme_bw() +
  ggtitle(label = "GO Biological Processes") +
  scale_x_continuous(breaks = 1:7, expand = expansion(mult = c(0,0.05))) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.title = element_text(color = "black", face = "bold", size = 13)
  )

ggsave("~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Zannas_lab/MRS TOPMed/GREAT_output/GO_dotplot.pdf", 
       p2, width = 10, height = 5)
