##### Loading required libraries
library(dplyr)
library(pathlinkR)

##### Documentation link
# #####https://bioconductor.org/packages/release/bioc/vignettes/pathlinkR/inst/doc/pathlinkR.html#supplemental-materials

# Building network from RNA-seq results
network <- ppiBuildNetwork(
  rnaseqResult=res[[1]],
  filterInput=F,
  order="zero"
)
write.csv(network, file = "network_output.csv")  # Writing network to CSV file

# Setting working directory and creating a PDF plot for the network
setwd("#####/network_analysis")
pdf("#####/network_visualization.pdf", height = 10, width = 10)
ppiPlotNetwork(
  network=network,
  title= "Network Visualization",
  fillColumn=LogFoldChange,
  foldChangeColours = c("#f7cac9", "#92a8d1"),
  fillType="foldChange",
  label=TRUE,
  hubColour = "#e86af0",
  labelColour = "#00c2c7",
  labelSize = 3,
  labelColumn=hgncSymbol,
  legend=TRUE
)
dev.off()

### Performing pathway enrichment analysis
enrichedResults <- pathwayEnrichment(
  inputList=res,
  analysis="sigora",
  pCutoff = 0.05,
  fcCutoff = 0.5,
  filterInput=F
)

head(enrichedResults)
write.csv(enrichedResults, file = "pathway_enrichment_results.csv")  # Saving enriched results to CSV

# Loading pathway database and calculating pathway distances
data("sigoraDatabase")
pathwayDistances <- getPathwayDistances(pathwayData = sigoraDatabase)

# Defining starting pathways based on calculated distances
foundationPathways <- pathnetFoundation(
  mat=pathwayDistances,
  maxDistance=0.8
)

# Filtering enriched pathways for network creation
filteredPathwayInput <- enrichedResults %>% 
  filter(comparison == "Group1 vs Group2")

pathwayNetwork <- pathnetCreate(
  pathwayEnrichmentResult=filteredPathwayInput,
  foundation=foundationPathways,
  trimOrder = 2
)

# Creating a PDF for pathway network visualization
pdf("#####/pathway_network_graph.pdf", height = 10, width = 10)
pathnetGGraph(
  pathwayNetwork,
  nodeSizeRange = c(2, 10),
  labelProp = 0.1,
  nodeLabelSize = 3,
  nodeLabelOverlaps = 8,
  segColour = "red",
  themeBaseSize = 12
)
dev.off()

# Additional pathway enrichment analysis with a different method
additionalResults <- pathwayEnrichment(
  inputList=res,
  analysis="hallmark",
  filterInput=F,
  split=TRUE
)
head(additionalResults)
write.csv(additionalResults, file = "additional_enrichment_results.csv")  # Saving additional results

# Generating PDF with pathway plots
pdf("#####/pathway_plots.pdf", height = 10, width = 10)
pathwayPlots(
  pathwayEnrichmentResults=additionalResults,
  columns=2,
  maxPVal=30,
  colourValues=c("#f7cac9", "#92a8d1")
)
dev.off()
