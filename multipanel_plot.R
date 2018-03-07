#-----------------------------------
#     Multipanel plots 
#-----------------------------------

install.packages('devtools')
devtools::install_bitbucket("graumannlabtools/multipanelfigure")

MAplot
pca

library(multipanelfigure)

# Specify a layout
figure1 <- multi_panel_figure(
  width = 150,
  height = 75,
  columns = 2,
  rows = 2
)
figure1

# Layout for 8.5 x 11 sheet of paper:
# 7" wide by 10" high : width = 
inches_to_mm <- function(inches){
  inches * 25.4
}
wid = inches_to_mm(4)
height = inches_to_mm(6)



figure1 <- multi_panel_figure(
  width = wid,
  height = height,
  columns = 2,
  rows = 3
)
figure1

# Set path to png output files
PCA_dir <- paste(outputdir,'/','PBS_STM_PCAplot.png', sep = '')
MA_dir <- paste(outputdir,'/','Styphimurium_over_PBS_MAplot.png', sep = '')
heatmap_dir <- paste(outputdir, '/','Styphimurium_over_PBSTop25VarianceGenes_heatmap.png', sep = '')
figure1 %<>% fill_panel(PCA_dir, row = 1, column = 1, scaling = 'fit')
figure1 %<>% fill_panel(MA_dir, column = 2, row = 1, scaling = 'fit')
figure1 %<>% fill_panel(heatmap_dir, column = 1, row = c(2,3), scaling = 'fit')
figure1

# Add plots to figure
library(magrittr)
MAplot
heatmap
testplot <- capture_base_plot()
figure1 %<>% fill_panel(testplot, column = 1)
figure1 %<>% fill_panel(pca, column = 2)

figure1

testplot <- capture_base_plot(
  plotMA(res, 
       ylim = c(-2,4), 
       main = paste(sample1, '/', sample2,'\n', 'MA plot'))
)
testplot



MAplotpdf <- paste(outputdir,'/',sample1,'_over_',sample2,'_MAplot.pdf', sep = '')
figure1 %<>% fill_panel(MAplotpdf, pdf, column = 2)
