library(ggplot2)
library(tidyr)
library(tools)
library(grid)
library(gridExtra)
library(lattice)
pdf(NULL)

# get all csv files in current directory that end in "trace.csv"
csv_files <- list.files("traces/", pattern = "trace.csv$", full.names = TRUE)

# loop through each csv file and save an individual plot for each
for (file in csv_files) {
    # read in the csv file with the first column as rownames
    print(paste0("Drawing fracment ion traces for: ", file))
    df <- read.csv(file, row.names = 1)

    # pivot the data to long format
    if (!any(grepl("^X", colnames(df)))) {
        next
    }

    df_long <- pivot_longer(df, cols = starts_with("X"), names_to = "column", values_to = "value")

    # define labels for each column (i.e. sample) for plotting to remove the X
    df_long$column <- gsub("X", "", df_long$column)

    # divide all values in as.numeric.trace_rts. column by 60 to convert from seconds to minutes
    df_long$as.numeric.trace_rts. <- df_long$as.numeric.trace_rts. / 60

    # create the plot
    g <- ggplot(df_long, aes(x = as.numeric.trace_rts., y = value, color = column)) +
        geom_line() +
        theme_classic(base_size = 18) +
        labs(x = "Retention Time (min)", y = "Signal", color = "Fragment m/z", title = strsplit(file, split = "_")[[1]][1])

    # save the plot to a file with the same name as the csv file
    ggsave(g, filename = file_path_sans_ext(file) %>% paste0(".png"), width = 6, height = 4, units = "in", dpi = 300)
}

# this is the name of the plots in the order they are made first time through
plotnames <- c(
    "12-Deoxy*",
    "Digoxin-d3",
    "Withanoside IV",
    "Withanoside V",
    "Withaferin A",
    "Withanolide A",
    "Withanolide B",
    "Withanone"
)

# combine all plots into one 4 by 4 grid
i <- 0
glist <- list()
for (file in csv_files) {
    i <- i + 1
    # read in the csv file with the first column as rownames
    print(paste0("Drawing fracment ion traces for: ", file))
    df <- read.csv(file, row.names = 1)

    # pivot the data to long format
    if (!any(grepl("^X", colnames(df)))) {
        next
    }

    df_long <- pivot_longer(df, cols = starts_with("X"), names_to = "column", values_to = "value")

    # define labels for each column (i.e. sample) for plotting to remove the X
    df_long$column <- gsub("X", "", df_long$column)

    # divide all values in as.numeric.trace_rts. column by 60 to convert from seconds to minutes
    df_long$as.numeric.trace_rts. <- df_long$as.numeric.trace_rts. / 60

    # create the plot
    
    glist[[i]] <- ggplot(df_long, aes(x = as.numeric.trace_rts., y = value, color = column)) +
        geom_line() +
        theme_classic(base_size = 18) +
        # labs(x = "Retention Time (min)", y = "Signal", color = "Fragment m/z", title = strsplit(file, split = "_")[[1]][1]) # for first run through to get all plots
        labs(x = "Retention Time (min)", y = "Signal", color = "Fragment m/z", title = plotnames[i]) # for second run through when plots are in order

    # save the plot to a file with the same name as the csv file
    # ggsave(g, filename = file_path_sans_ext(file) %>% paste0(".png"), width = 6, height = 4, units = "in", dpi = 300)
}

# For the first run through to get all plots
# gridplots <- grid.arrange(grobs = glist, ncol = 4, nrow = 2)
# ggsave("prm_traces.png", gridplots, width = 25, height = 8, dpi = 300)


# reorder the plots
plotorder <- c(3,4,5,1,6,8,7,2)
gridplots <- grid.arrange(grobs = glist[plotorder], ncol = 4, nrow = 2)
ggsave(paste0("traces/", "prm_traces.png"), gridplots, width = 25, height = 8, dpi = 300)
