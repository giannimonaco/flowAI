# Retrieves the number of events for each FCS file and creates
# a bar plot.
#
flow_set_qc <- function(set){
  if (!is(set, "flowSet"))
    stop("'set' needs to be of class 'flowSet'")

  cellNumbers <- as.numeric(fsApply(set, nrow))
  cellNumbers[cellNumbers == 1] <- NA
  cnframe <- data.frame(sampleName = sampleNames(set), cellNumber = cellNumbers)
  return(cnframe)
}

# produce a bar plot
flow_set_plot <- function(N_cell_set, area){

  ggplot(N_cell_set, aes(x=sampleName, y = cellNumber)) +
  geom_bar(stat = "identity",fill= area) + theme_classic() +
  theme( legend.position="none", axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
}
