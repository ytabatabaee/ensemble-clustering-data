library(ggplot2)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(plyr)
library(lubridate)
library(scales)

degree_df <- read.csv('./output/average_degree', stringsAsFactors=FALSE, na.strings=c("None"))

# p_degree <- ggplot(degree_df, aes(x=size, y=average_degree)) + geom_point() + theme_bw() + theme(axis.title.y=element_blank()) + theme(axis.title.x=element_blank()) + theme(text=element_text(size=10)) + theme(legend.position="bottom") + theme(legend.title=element_blank()) + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
# p_degree <- ggplot(degree_df, aes(x=size, y=average_degree)) + geom_point() + theme_bw() + ylab("Average degree") + xlab("Number of nodes") + theme(text=element_text(size=10)) + theme(legend.position="bottom") + theme(legend.title=element_blank()) + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) + scale_y_continuous(trans="log10") + scale_x_continuous(trans="log10", labels = label_comma()) + geom_hline(yintercept=20, linetype="dashed", color = "red") + geom_vline(xintercept=10000, linetype="dashed", color = "blue")
# p_degree <- ggplot(degree_df, aes(x=size, y=average_degree)) + geom_point() + theme_bw() + ylab("Average degree") + xlab("Number of nodes") + theme(text=element_text(size=10)) + theme(legend.position="bottom") + theme(legend.title=element_blank()) + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) + scale_y_continuous(labels=label_comma()) + scale_x_continuous(trans="log10", labels = label_comma()) + geom_hline(yintercept=20, linetype="dashed", color = "red") + geom_vline(xintercept=10000, linetype="dashed", color = "blue")
p_degree <- ggplot(degree_df, aes(x=size, y=average_degree)) + geom_point() + theme_bw() + ylab("Average degree") + xlab("Number of nodes") + theme(text=element_text(size=10)) + theme(legend.position="bottom") + theme(legend.title=element_blank()) + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) + scale_y_continuous(labels=label_comma()) + scale_x_continuous(labels = label_comma()) + geom_hline(yintercept=20, linetype="dashed", color = "red") + geom_vline(xintercept=10000, linetype="dashed", color = "blue")
# p_degree <- ggplot(degree_df, aes(x=size, y=average_degree)) + geom_point() + theme_bw() + ylab("Average degree") + xlab("Number of nodes") + theme(text=element_text(size=10)) + theme(legend.position="bottom") + theme(legend.title=element_blank()) + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) + scale_y_continuous(trans="log10") + scale_x_continuous(labels = label_comma()) + geom_hline(yintercept=20, linetype="dashed", color = "red") + geom_vline(xintercept=10000, linetype="dashed", color = "blue")


fig_name = paste('./output/degree.pdf', sep="")
pdf(fig_name)
fig_name = paste('./output/degree.eps', sep="")
ggsave(height=180, width=122, units="mm", file=fig_name)
print(p_degree)
dev.off()
