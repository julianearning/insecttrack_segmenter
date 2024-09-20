library(scales)
df_loess <- read.table("build/loess.csv", sep=",", header=F)
df_points <- read.table("build/downsampled.csv", sep=",", header=F)


plot(c(0,1280), c(0,720),col=c('grey', 0.4), type="n")
for(cluster_idx in unique(df_points$V4)) {
#cluster_idx=1
	if(cluster_idx!=-1) {
	print(cluster_idx)

	df_loess1 <- df_loess[df_loess$V4==cluster_idx,]
	n <- nrow(df_loess1)

	df1 <- df_points[df_points$V4==cluster_idx,]
	points(df1$V1, df1$V2, col=c('grey',0.4))

	segments(df_loess1$V1[1:(n-1)], df_loess1$V2[1:(n-1)], df_loess1$V1[2:n], df_loess1$V2[2:n], col=(df_loess1$V4+4), lwd=3)
	
	}
}
