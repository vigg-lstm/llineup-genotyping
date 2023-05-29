library(data.table)
library(ggplot2)
library(stringr)
library(glmmTMB)
library(grid)
library(gridExtra)
library(gtable)

# Load the genotyping data
genotypes <- fread('../../data/ento_geno_plasmo_data.csv', na.strings = c('', NA))

# One thing to note is that there is very little repeat capture within a household. Only 19 households ever 
# produced mosquitoes on more than one occasion, with the max being 2 occasions:
#positive.collections.per.house <- genotypes[, .(collections = length(unique(RND))), by = 'HHID.rnd'][order(collections)]

# When doing it be EA on the other hand, you do get ones with data from all 5 rounds. But most are still just 
# with one round. 
#positive.collections.per.EA <- genotypes[, .(collections = length(unique(RND))), by = c('HSD', 'EA')][order(collections)]

# Recode genotypes
genotypes[, ':='('Vgsc-995F' = str_count(kdr, 'F'),
                 'Vgsc-995S' = str_count(kdr, 'S'),
                 'Cyp4j5-43F' = str_count(Cyp4j5, 'F'),
                 'Coeae1d' = str_count(Coeae1d, 'R'),
                 'Cyp6p4-236M' = str_count(Cyp6p4, 'M'),
                 'ZZB-TE' = str_count(ZZB, 'Z'),
                 'Cyp6aap-Dup1' = str_count(Cyp6aa, 'D'),
				 'Chr2La' = str_count(X2La, 'I')
            )
]

# We will turn the genotypes to haplotypes
gen2hap <- function(n){
	if (is.na(n))
		c(NA, NA)
	else
		c(rep(0, 2-n), rep(1, n))
}

markers <- c(kdrF = 'Vgsc-995F',
             kdrS = 'Vgsc-995S',
             Cyp4j5 = 'Cyp4j5-43F',
             Coeaed1 = 'Coeae1d',
             Cyp6p4 = 'Cyp6p4-236M',
             ZZB = 'ZZB-TE',
             Dup1 = 'Cyp6aap-Dup1', 
             La2 = 'Chr2La')

gam.haplotypes <- genotypes[MolSpecies == 'S', c(lapply(.SD, function(x) as.vector(sapply(x, gen2hap))),
                                                 .(RND = rep(RND, each = 2),
                                                   HSD = rep(HSD, each = 2),
                                                   EA = rep(EA, each = 2),
                                                   HHID.rnd = rep(HHID.rnd, each = 2),
                                                   Wave = rep(Wave, each = 2),
                                                   LLIN.actual = rep(LLIN.actual, each = 2),
                                                   Net.actual = rep(Net.actual, each = 2),
                                                   Location = rep(Location, each = 2),
                                                   wgs.sample.id = rep(wgs.sample.id, each = 2)
                                                 )
                                               ),
                                               .SDcols = markers
]

# There seem to be a lot of NAs in the LLIN.actual column. Due to the clusters which received the wrong
# net allocation. HSDs 27 and 78
gam.haplotypes<- subset(gam.haplotypes, LLIN.actual == "PBO LLIN" | LLIN.actual == "Non-PBO LLIN")

# Plot data by arm
plot.means.and.conf <- function(snp.name, haps) {
	means.and.conf <- haps[, 
                           as.list(unlist(prop.test(sum(.SD[[snp.name]], na.rm = T), 
	                                                sum(!is.na(.SD[[snp.name]]))
	                                      )[c('estimate', 'conf.int')]
	                               )
	                       ),
                           by = c('RND','LLIN.actual')
                      ]
	y.lims <- c(min(means.and.conf[, 3:ncol(means.and.conf)]), max(means.and.conf[, 3:ncol(means.and.conf)]))
	# Empirically, we have observed that the largest range is 0.55. So impose this range on all markers
	# so that we can zoom in on the results while keeping the range constant between markers.
	max.range <- 0.55
	y.range <- y.lims[2] - y.lims[1]
	y.expansion <- (max.range - y.range)/2
	y.lims <- y.lims + c(-y.expansion, y.expansion)
	if (y.lims[1] < 0){
		y.lims[2] <- y.lims[2] - y.lims[1]
		y.lims[1] <- 0
	}
	else if (y.lims[2] > 1){
		y.lims[1] <- y.lims[1] - (y.lims[2] - 1)
		y.lims[2] <- 1
	}
	# Now draw the actual plot
	pp <- ggplot(means.and.conf, aes(x = RND, y = estimate.p, color = LLIN.actual, fill = LLIN.actual)) + 
                 labs(y = 'Allele frequency',title=snp.name ) +
                 geom_polygon(data = rbind(means.and.conf[, .(RND, LLIN.actual, conf = conf.int1)], 
                                           means.and.conf[nrow(means.and.conf):1, .(RND, LLIN.actual, conf = conf.int2)]
                                     ), 
                              aes(x = RND, y = conf),
                              alpha = 0.5,
                              size = 0) +
                 geom_line(size = 1.2) +
	               geom_point(pch = 21, size = 1.8, color = 'grey60') +
	             coord_cartesian(ylim = y.lims) +
                 theme_classic() +
	               theme(legend.title = element_blank())+
                 theme(plot.title = element_text(hjust = 0.5)) 
	pp 
}


combined.plots <- lapply(markers, plot.means.and.conf, haps = gam.haplotypes)
combined.legend = gtable_filter(ggplotGrob(combined.plots[[1]]), "guide-box")

###########
svg('Lynd et al Figure 3_v26May23.svg', width = 12, height = 6)
grid.arrange(
	do.call(arrangeGrob, c(lapply(combined.plots, function(p) p + theme(legend.position='none') + labs(x = '', y = '')), 
	                       nrow = 2,
	                       top = list(textGrob('', gp = gpar(fontface = "bold", cex = 1.3))),
	                       left = list(textGrob('Allele frequency', rot = 90, vjust = 1)),
	                       bottom = list(textGrob('Round', vjust = -1)))
	),
	combined.legend,
    widths=unit.c(unit(1, "npc") - combined.legend$width, combined.legend$width), 
	nrow = 1
)
dev.off()

##############################################################################
# Models being tested:
# H1. Markers will increase in frequency over the course of the intervention
  # Model 1- glmmTMB (marker ~ RND + Location  + LLIN.actual +(1|HSD)) - round as continuous
  # Model 2- glmmTMB (marker ~ RND(1vs5) + Location+ LLIN.actual (BINARY) +(1|HSD))) - baseline vs 25 month
# H2. The increase in cytochrome P450 resistance markers will be lower in An. gambiae s.s. populations 
# from clusters which received PBO-LLINs (Same model as above replicated here in case changed)
  # Model 3- glmmTMB (marker ~ RND + Location + LLIN.actual + LLIN.actual:RND  +(1|HSD))) 
  # Model 4- glmmTMB (marker ~ RND (1vs5) + Location + LLIN.actual + LLIN.actual:RND (1vs5) +(1|HSD))) 
# H3. Dup1/4j5 haplotype will show a greater increase in clusters with Permanet 2 vs Olyset
  # Model 5- glmmTMB (marker ~ RND + Location + P2vsOly + P2vsOly:RND +(1|HSD))) 
  # Model 6- glmmTMB (marker ~ RND (1vs5) + Location + P2vsOly +P2vsOly:RND + (1|HSD))) 

haps <- gam.haplotypes

get.term.details <- function(model, test, model.term = 'RND', test.term = model.term){
	p <- test[['Pr(>Chi)']][which(attributes(test)$row.names == test.term)]
	coefficient <- fixef(model)$cond[model.term]
	c(P = signif(p, 2), coef = signif(coefficient, 2))
}

run.models <- function(marker, haps, verbose = T){
	if (verbose){
		cat('\nRunning models on', marker, '\n')
		cat('\n\tHypo 1_cont:\n')
	}
	
	# Hypothesis 1 
	# model.1 is model for H1 continuous 
	model.1 <- glmmTMB(get(marker) ~ RND + LLIN.actual + Location + (1|HSD), 
	                   data = haps, family = 'binomial')
	test.model.1 <- drop1(model.1, test = 'Chisq')
	if (verbose){
		print(model.1)
		print(test.model.1)
		cat('\n\tHypo 1_ base vs 25:\n')
	}
	# model.2 H1 BAseline vs 25 month 
	model.2 <- glmmTMB(get(marker) ~ as.factor(RND) + LLIN.actual + Location + (1|HSD),
	                   data = haps[RND %in% c(1,5)], family = 'binomial')
	test.model.2 <- drop1(model.2, test = 'Chisq')
	if (verbose){
		print(model.2)
		print(test.model.2)
		cat('\n\tHypo2_cont:\n')
	}
	
	# Hypothesis 2
	# model.3 is model for H2 continuous 
	model.3 <- glmmTMB(get(marker) ~ RND + Location + LLIN.actual + RND:LLIN.actual + (1|HSD),
	                   data = haps[!is.na(LLIN.actual)], family = 'binomial')
	test.model.3 <- drop1(model.3, test = 'Chisq')
	if (verbose){
		print(model.3)
		print(test.model.3)
		cat('\n\tHypo2_ base vs 25:\n')
	}
	
	# model.4 H2 BAseline vs 25 month 
	model.4 <- glmmTMB(get(marker) ~ as.factor(RND) + Location + LLIN.actual + as.factor(RND):LLIN.actual + (1|HSD),
	                   data = haps[RND %in% c(1,5)], family = 'binomial')
	test.model.4 <- drop1(model.4, test = 'Chisq')
	if (verbose){
	  print(model.4)
	  print(test.model.4)
	  cat('\n\tHypo3_cont:\n')
	}
	
	# Hypothesis 3
	# model.5 is model for H3 continuous 
	model.5 <- glmmTMB(get(marker) ~ RND + Net.actual + Location + RND:Net.actual+  (1|HSD),
	                   data = haps[LLIN.actual == 'Non-PBO LLIN'], family = 'binomial')
	test.model.5 <- drop1(model.5, test = 'Chisq')
	if (verbose){
	  print(model.5)
	  print(test.model.5)
	  cat('\n\tHypo3_ base vs 25:\n')
	}
	
	# model.6 H3 BAseline vs 25 month 
	model.6 <- glmmTMB(get(marker) ~ as.factor(RND) + Net.actual + Location + as.factor(RND):Net.actual + (1|HSD),
	                          data = haps[LLIN.actual == 'Non-PBO LLIN' & RND %in% c(1,5)], 
	                          family = 'binomial')
	test.model.6 <- drop1(model.6, test = 'Chisq')
	if (verbose){
	print(model.6)
	print(test.model.6)
	}
	
	model.names <- c('model1', 'model2','model3','model4', 'model5', 'model6')
	
	model.terms <- c('RND', 'as.factor(RND)5',
	                 'RND:LLIN.actualPBO LLIN', 'as.factor(RND)5:LLIN.actualPBO LLIN',
	                 'RND:Net.actualPermaNet-2', 'as.factor(RND)5:Net.actualPermaNet-2') %>%
				   setNames(model.names)
	test.terms <- c('RND', 'as.factor(RND)', 
	                'RND:LLIN.actual', 'as.factor(RND):LLIN.actual',
	                'RND:Net.actual', 'as.factor(RND):Net.actual') %>%
				  setNames(model.names)
	
	all.models <- list(model.1,
	                   model.2,
	                   model.3,
	                   model.4,
	                   model.5,
	                   model.6) %>%
	              setNames(model.names)
	all.tests <- list(test.model.1,
					  test.model.2,
					  test.model.3,
					  test.model.4,
					  test.model.5,
					  test.model.6) %>%
	             setNames(model.names)
	test.summary = t((mapply(get.term.details, all.models, all.tests, model.terms, test.terms)))
	list(all.models = all.models, all.tests = all.tests, test.summary = test.summary)
}

all.glms <- lapply(markers, run.models, gam.haplotypes)

# Print out all of the summaries to the logs
cat('\n\nSummary of results:\n\n')
for (marker in names(markers)){
	cat('\n\t', marker, '\n\n')
	print(all.glms[[marker]])
	cat('\n')
}





