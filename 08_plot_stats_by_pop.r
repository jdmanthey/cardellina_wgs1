options(scipen=999)
library(Hmisc)

dxy <- read.table("window_dxy.txt", header=T)
fst <- read.table("window_fst.txt", header=T)
pi <- read.table("window_pi.txt", header=T)
pi$calculated_stat[is.na(pi$calculated_stat)] <- 0
tajima <- read.table("window_tajima.txt", header=T)
# combine all
x <- rbind(dxy, fst, pi, tajima)

# determine order of windows
scaffolds <- c("CM027507.1","CM027536.1","CM027508.1","CM027509.1","CM027510.1","CM027537.1","CM027511.1","CM027512.1","CM027513.1","CM027514.1","CM027515.1","CM027516.1","CM027517.1","CM027518.1","CM027519.1","CM027520.1","CM027521.1","CM027522.1","CM027523.1","CM027524.1","CM027525.1","CM027526.1","CM027527.1","CM027528.1","CM027529.1","CM027531.1","CM027532.1","CM027533.1","CM027535.1")
scaffold_names <- c("Chr1", "Chr1A", "Chr2", "Chr3", "Chr4", "Chr4A", "Chr5", "Chr6", "Chr7", "Chr8","Chr9", "Chr10","Chr11", "Chr12","Chr13", "Chr14","Chr15", "Chr17","Chr18", "Chr19","Chr20", "Chr21","Chr22", "Chr23","Chr24", "Chr26", "Chr27","Chr28", "ChrZ")
window_order <- c()
for(a in 1:length(scaffolds)) {
	a_rep <- x[x$chr == scaffolds[a], ]
	a_rep <- a_rep[match(unique(a_rep$start), a_rep$start),]
	a_rep <- a_rep[order(a_rep$start),]
	window_order <- rbind(window_order, a_rep)
}
window_order <- window_order[,4:6]

stat_order <- rbind(c("SMOcc", "TMVB", "Dxy"),
					c("SMOcc", "SMS", "Dxy"),
					c("SMS", "TMVB", "Dxy"),
					c("OH", "TMVB", "Dxy"),
					c("OH", "SMS", "Dxy"),
					c("OH", "XISTH", "Dxy"),
					c("SMOcc", "TMVB", "Fst"),
					c("SMOcc", "SMS", "Fst"),
					c("SMS", "TMVB", "Fst"),
					c("OH", "TMVB", "Fst"),
					c("OH", "SMS", "Fst"),
					c("OH", "XISTH", "Fst"),
					c("SMOcc", "none", "pi"),
					c("TMVB", "none", "pi"),
					c("SMS", "none", "pi"),
					c("OH", "none", "pi"),
					c("XISTH", "none", "pi"),
					c("SMOcc", "none", "Tajima_D"),
					c("TMVB", "none", "Tajima_D"),
					c("SMS", "none", "Tajima_D"),
					c("OH", "none", "Tajima_D"),
					c("XISTH", "none", "Tajima_D"))
					
					
# calculate mean values for each stat of interest
for(a in 1:17) {
	if(a >= 1 & a <= 6) {
		a_rep <- x[x$pop1 == stat_order[a,1] & x$pop2 == stat_order[a,2] & x$stat == stat_order[a,3],]
		a_rep <- sum(a_rep$calculated_stat * a_rep$number_sites) / sum(a_rep$number_sites)
		writeLines(paste(stat_order[a,1], stat_order[a,2], stat_order[a,3], a_rep, sep="\t\t"))
	} else if(a >= 7 & a <= 12) {
		a_rep <- x[x$pop1 == stat_order[a,1] & x$pop2 == stat_order[a,2] & x$stat == stat_order[a,3],]
		a_rep <- sum(a_rep$calculated_stat * a_rep$number_variable_sites) / sum(a_rep$number_variable_sites)
		writeLines(paste(stat_order[a,1], stat_order[a,2], stat_order[a,3], a_rep, sep="\t\t"))
	} else if(a >= 13 & a <= 17) {
		a_rep <- x[x$pop1 == stat_order[a,1] & x$pop2 == stat_order[a,2] & x$stat == stat_order[a,3],]
		a_rep <- sum(a_rep$calculated_stat * a_rep$number_sites) / sum(a_rep$number_sites)
		writeLines(paste(stat_order[a,1], stat_order[a,2], stat_order[a,3], a_rep, sep="\t\t"))
	}
}


# loop for each window to summarize
stat_adds <- list()
for(a in 1:nrow(window_order)) {
	if(a %% 1000 == 0) { print(a) }
	a_rep <- x[x$chr == window_order$chr[a] & x$start == window_order$start[a] & x$end == window_order$end[a],]
	stat_adds_rep <- c()
	# for each stat
	for(b in 1:nrow(stat_order)) {
		b_rep <- a_rep[a_rep[,1] == stat_order[b,1] & a_rep[,2] == stat_order[b,2] & a_rep[,3] == stat_order[b,3],]
		if(nrow(b_rep) == 1) {
			stat_adds_rep <- c(stat_adds_rep, b_rep[,9])
		} else {
			stat_adds_rep <- c(stat_adds_rep, NA)
		}
	}
	stat_adds[[a]] <- stat_adds_rep
}

#unlist stat_adds
stat_adds2 <- do.call(rbind.data.frame, stat_adds)
colnames(stat_adds2) <- c("dxy_SMOcc_TMVB", "dxy_SMOcc_SMS", "dxy_SMS_TMVB", "dxy_TMVB_OH", "dxy_SMS_OH", "dxy_OH_XISTH", "fst_SMOcc_TMVB", "fst_SMOcc_SMS", "fst_SMS_TMVB", "fst_TMVB_OH", "fst_SMS_OH", "fst_OH_XISTH", "pi_SMOcc", "pi_TMVB", "pi_SMS", "pi_OH", "pi_XISTH", "TajD_SMOcc", "TajD_TMVB", "TajD_SMS", "TajD_OH", "TajD_XISTH")

# combine data frames
x <- cbind(window_order, stat_adds2)
x <- na.omit(x)
rownames(x) <- seq(from=1, to=nrow(x), by=1)

# correlations of DXY and FST
 rcorr(as.matrix(x[,10:15]))
 rcorr(as.matrix(x[,4:9]))

# windows
windows <- as.numeric(rownames(x))
# sliding window size
# = 20 50kbp windows = 1 Mbp
window_size <- 20
# set up row numbers for line plots
line_rows <- seq(from=1,to=nrow(x), by=1)[seq(from=1,to=nrow(x), by=1) %% window_size == 0]

# what are the unique chromosomes?
chr <- unique(x$chr)
scaffold <- chr


# define population order and naming conventions
population_order <- c("SMOcc", "TMVB", "SMS", "OH", "XISTH")
population_names <- c("SMOcc", "TMVB", "SMS", "OH", "XISTH")

# define population comparisons order and naming conventions
pop_split_order <- c("SMOcc_TMVB", "SMOcc_SMS", "SMS_TMVB", "TMVB_OH", "SMS_OH", "OH_XISTH")
pop_split_names <- c("SMOcc x TMVB", "SMOcc x SMS", "TMVB x SMS", "TMVB x OH", "SMS x OH", "OH x XISTH")

# ylim boundaries for each stat
pi_ylim <- c(0, 0.008)
fst_ylim <- c(0, 0.8)
dxy_ylim <- c(0,0.012)
tajd_ylim <- c(-2, 2)

#########################################################################
#########################################################################
#########################################################################
# Plot Dxy
#########################################################################
#########################################################################
#########################################################################

# make the plotting polygons
chr_polygons <- list()
for(a in 1:length(chr)) {
	a1 <- windows[x$chr == chr[a]]
	a2 <- a1[length(a1)]
	a1 <- a1[1]
	chr_polygons[[a]] <- rbind(c(a1, dxy_ylim[1]), c(a2, dxy_ylim[1]), c(a2, dxy_ylim[2]), c(a1, dxy_ylim[2]), c(a1, dxy_ylim[1]))
}

# loop for each plot
pdf("_window_DXY.pdf", width=6.5, height=6.5)
par(mfrow=c(length(pop_split_order),1))
par(mar=c(2, 2, 0.1, 0.1))
for(b in 1:length(pop_split_order)) {
	# plot the chromosomes
	plot(c(-1,-1), ylim=c(dxy_ylim[1], dxy_ylim[2]), xlim=c(1, length(windows)), xaxt="n", col="white", bty="n", cex.axis=1.1, cex.lab=1.3, ylab="", xlab="")
	odd <- 0
	for(a in 1:length(chr_polygons)) {
		if(odd == 1) {
			polygon(chr_polygons[[a]], col="snow2", border="white")
			odd <- 0	
		} else {
			odd <- 1
		}
	}

	# plot the stat
	# subset this individual or comparison
	b_rep <- data.frame(chromosome=as.character(x$chr), start=as.numeric(x$start), stat=as.numeric(x[,colnames(x) == paste("dxy_", pop_split_order[b], sep="")]), window=as.numeric(rownames(x)))
	# plot
	points(rownames(b_rep), b_rep$stat, pch=19, cex=0.05, col="gray65")
	title(paste(pop_split_names[b], sep=""), adj=0.01, line=-1)
	
	# plot sliding mean line plots
	line_x_axis <- line_rows
	line_y_axis <- c()
	line_scaffold <- c()
	for(e in line_rows) {
		line_y_axis <- c(line_y_axis, mean(b_rep$stat[b_rep$window %in% (b_rep$window[e] - (floor(window_size / 2))):(b_rep$window[e] + (floor(window_size / 2))) & b_rep$chromosome == b_rep$chromosome[e]][1:window_size]))
		line_scaffold <- c(line_scaffold, b_rep$chromosome[e])
	}
	line_plotting <- data.frame(line_x_axis=as.numeric(line_x_axis), line_y_axis=as.numeric(unlist(line_y_axis)), line_scaffold=as.character(unlist(line_scaffold)))
	# plot each scaffold at a time (so lines don't connect between scaffolds)
	for(e in 1:length(unique(scaffold))) {
		a_rep <- line_plotting[line_plotting[,3] == unique(scaffold)[e],]
		lines(a_rep[,1:2], lwd=0.8, col="gray29")
	}


}
dev.off()





#########################################################################
#########################################################################
#########################################################################
# Plot Fst
#########################################################################
#########################################################################
#########################################################################

# make the plotting polygons
chr_polygons <- list()
for(a in 1:length(chr)) {
	a1 <- windows[x$chr == chr[a]]
	a2 <- a1[length(a1)]
	a1 <- a1[1]
	chr_polygons[[a]] <- rbind(c(a1, fst_ylim[1]), c(a2, fst_ylim[1]), c(a2, fst_ylim[2]), c(a1, fst_ylim[2]), c(a1, fst_ylim[1]))
}

# loop for each plot
pdf("_window_FST.pdf", width=6.5, height=6.5)
par(mfrow=c(length(pop_split_order),1))
par(mar=c(2, 2, 0.1, 0.1))
for(b in 1:length(pop_split_order)) {
	# plot the chromosomes
	plot(c(-1,-1), ylim=c(fst_ylim[1], fst_ylim[2]), xlim=c(1, length(windows)), xaxt="n", col="white", bty="n", cex.axis=1.1, cex.lab=1.3, ylab="", xlab="")
	odd <- 0
	for(a in 1:length(chr_polygons)) {
		if(odd == 1) {
			polygon(chr_polygons[[a]], col="snow2", border="white")
			odd <- 0	
		} else {
			odd <- 1
		}
	}

	# plot the stat
	# subset this individual or comparison
	b_rep <- data.frame(chromosome=as.character(x$chr), start=as.numeric(x$start), stat=as.numeric(x[,colnames(x) == paste("fst_", pop_split_order[b], sep="")]), window=as.numeric(rownames(x)))

	# plot
	points(rownames(b_rep), b_rep$stat, pch=19, cex=0.05, col="gray65")
	title(paste(pop_split_names[b], sep=""), adj=0.01, line=-1)
		
	# plot sliding mean line plots
	line_x_axis <- line_rows
	line_y_axis <- c()
	line_scaffold <- c()
	for(e in line_rows) {
		line_y_axis <- c(line_y_axis, mean(b_rep$stat[b_rep$window %in% (b_rep$window[e] - (floor(window_size / 2))):(b_rep$window[e] + (floor(window_size / 2))) & b_rep$chromosome == b_rep$chromosome[e]][1:window_size]))
		line_scaffold <- c(line_scaffold, b_rep$chromosome[e])
	}
	line_plotting <- data.frame(line_x_axis=as.numeric(line_x_axis), line_y_axis=as.numeric(unlist(line_y_axis)), line_scaffold=as.character(unlist(line_scaffold)))
	# plot each scaffold at a time (so lines don't connect between scaffolds)
	for(e in 1:length(unique(scaffold))) {
		a_rep <- line_plotting[line_plotting[,3] == unique(scaffold)[e],]
		lines(a_rep[,1:2], lwd=0.8, col="gray29")
	}
}
dev.off()


#########################################################################
#########################################################################
#########################################################################
# Plot pi
#########################################################################
#########################################################################
#########################################################################

# make the plotting polygons
chr_polygons <- list()
for(a in 1:length(chr)) {
	a1 <- windows[x$chr == chr[a]]
	a2 <- a1[length(a1)]
	a1 <- a1[1]
	chr_polygons[[a]] <- rbind(c(a1, pi_ylim[1]), c(a2, pi_ylim[1]), c(a2, pi_ylim[2]), c(a1, pi_ylim[2]), c(a1, pi_ylim[1]))
}

# loop for each plot
pdf("_window_pi.pdf", width=6.5, height=6.5)
par(mfrow=c(length(population_order),1))
par(mar=c(2, 2, 0.1, 0.1))
for(b in 1:length(population_order)) {
	# plot the chromosomes
	plot(c(-1,-1), ylim=c(pi_ylim[1],pi_ylim[2]), xlim=c(1, length(windows)), xaxt="n", col="white", bty="n", cex.axis=1.1, cex.lab=1.3, ylab="", xlab="")
	odd <- 0
	for(a in 1:length(chr_polygons)) {
		if(odd == 1) {
			polygon(chr_polygons[[a]], col="snow2", border="white")
			odd <- 0	
		} else {
			odd <- 1
		}
	}

	# plot the stat
	# subset this individual or comparison
	b_rep <- data.frame(chromosome=as.character(x$chr), start=as.numeric(x$start), stat=as.numeric(x[,colnames(x) == paste("pi_", population_order[b], sep="")]), window=as.numeric(rownames(x)))

	# plot
	points(rownames(b_rep), b_rep$stat, pch=19, cex=0.05, col="gray65")
	title(paste(population_names[b], sep=""), adj=0.01, line=-1)
		
	# plot sliding mean line plots
	line_x_axis <- line_rows
	line_y_axis <- c()
	line_scaffold <- c()
	for(e in line_rows) {
		line_y_axis <- c(line_y_axis, mean(b_rep$stat[b_rep$window %in% (b_rep$window[e] - (floor(window_size / 2))):(b_rep$window[e] + (floor(window_size / 2))) & b_rep$chromosome == b_rep$chromosome[e]][1:window_size]))
		line_scaffold <- c(line_scaffold, b_rep$chromosome[e])
	}
	line_plotting <- data.frame(line_x_axis=as.numeric(line_x_axis), line_y_axis=as.numeric(unlist(line_y_axis)), line_scaffold=as.character(unlist(line_scaffold)))
	# plot each scaffold at a time (so lines don't connect between scaffolds)
	for(e in 1:length(unique(scaffold))) {
		a_rep <- line_plotting[line_plotting[,3] == unique(scaffold)[e],]
		lines(a_rep[,1:2], lwd=0.8, col="gray29")
	}
}
dev.off()


#########################################################################
#########################################################################
#########################################################################
# Plot Tajima D
#########################################################################
#########################################################################
#########################################################################

# make the plotting polygons
chr_polygons <- list()
for(a in 1:length(chr)) {
	a1 <- windows[x$chr == chr[a]]
	a2 <- a1[length(a1)]
	a1 <- a1[1]
	chr_polygons[[a]] <- rbind(c(a1, tajd_ylim[1]), c(a2, tajd_ylim[1]), c(a2, tajd_ylim[2]), c(a1, tajd_ylim[2]), c(a1, tajd_ylim[1]))
}

# loop for each plot
pdf("_window_TajimaD.pdf", width=6.5, height=6.5)
par(mfrow=c(length(population_order),1))
par(mar=c(2, 2, 0.1, 0.1))
for(b in 1:length(population_order)) {
	# plot the chromosomes
	plot(c(-1,-1), ylim=c(tajd_ylim[1], tajd_ylim[2]), xlim=c(1, length(windows)), xaxt="n", col="white", bty="n", cex.axis=1.1, cex.lab=1.3, ylab="", xlab="")
	odd <- 0
	for(a in 1:length(chr_polygons)) {
		if(odd == 1) {
			polygon(chr_polygons[[a]], col="snow2", border="white")
			odd <- 0	
		} else {
			odd <- 1
		}
	}

	# plot the stat
	# subset this individual or comparison
	b_rep <- data.frame(chromosome=as.character(x$chr), start=as.numeric(x$start), stat=as.numeric(x[,colnames(x) == paste("TajD_", population_order[b], sep="")]), window=as.numeric(rownames(x)))
	# plot
	points(rownames(b_rep), b_rep$stat, pch=19, cex=0.05, col="gray65")
	title(paste(population_names[b], sep=""), adj=0.01, line=-1)
		
	# plot sliding mean line plots
	line_x_axis <- line_rows
	line_y_axis <- c()
	line_scaffold <- c()
	for(e in line_rows) {
		line_y_axis <- c(line_y_axis, mean(b_rep$stat[b_rep$window %in% (b_rep$window[e] - (floor(window_size / 2))):(b_rep$window[e] + (floor(window_size / 2))) & b_rep$chromosome == b_rep$chromosome[e]][1:window_size]))
		line_scaffold <- c(line_scaffold, b_rep$chromosome[e])
	}
	line_plotting <- data.frame(line_x_axis=as.numeric(line_x_axis), line_y_axis=as.numeric(unlist(line_y_axis)), line_scaffold=as.character(unlist(line_scaffold)))
	# plot each scaffold at a time (so lines don't connect between scaffolds)
	for(e in 1:length(unique(scaffold))) {
		a_rep <- line_plotting[line_plotting[,3] == unique(scaffold)[e],]
		lines(a_rep[,1:2], lwd=0.8, col="gray29")
	}
}
dev.off()


