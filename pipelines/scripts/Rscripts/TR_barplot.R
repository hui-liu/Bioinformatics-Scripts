argv <- commandArgs(TRUE)
inFile <- argv[1]
Title <- argv[2]
outFile <- argv[3]

# load data
df <- read.table(inFile)
df <- df[,c(1,6,7)]
names(df) <- c("Chr", "Length", "Frequency")
df <- df[df$Length > 100,]
df <- df[grepl("Chr", df$Chr), ]
df$Chr <- factor(df$Chr, levels = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "Chr8",
                                    "Chr9", "Chr10", "Chr11", "Chr12"))

# plot
library(ggplot2)
pdf(outFile, width=8,height=6)
ggplot(data=df, aes(x=Length, y=Frequency)) +
  geom_bar(stat="identity", fill="blue") +
  geom_vline(xintercept = c(146, 146*2, 146*3), linetype="dotted", 
             color = "red", linewidth=0.5) +
  labs(title = Title) +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(family="serif", size=9, face = "bold", colour = "black"),
    panel.background = element_blank(),
    plot.title = element_text(family="serif", size=9, face = "bold.italic", colour = "black", hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(family="serif", size=9, face = "bold", colour = "black"),
    axis.text.x = element_text(family="serif", size=9, face = "bold", colour = "black"),
    axis.text.y= element_text(family="serif", size=9, face = "bold", colour = "black"),
    strip.text.x = element_text(family="serif", size=9, face = "bold", colour = "black")
  ) +
  facet_wrap(.~Chr, ncol = 4, scales = "free_y")
dev.off()

