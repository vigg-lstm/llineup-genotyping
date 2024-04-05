library(data.table)
library(ggplot2)
library(stringr)
library(glmmTMB)
library(grid)
library(gridExtra)
library(gtable)
library(bit64)

#######################################################
# Load the genotyping data
ento.data <- fread('../data/ento_geno_plasmo_data.csv', na.strings = c('', 'NA'))

# Keep only arabiensis, funestus and gambiae ss
ento.data.afg <- ento.data[MolSpecies %in% c('AR', 'FUN', 'S'), ]
cat('\nOverall counts of Pfal positive samples:\n')
print(with(ento.data.afg, table(species, Pfal)))

# Number to report in Manuscript: proportion of mosquitoes sporozoite positive:
sporpos <- ento.data[, .(Pfal.rate = round(sum(Pfal) / .N, 4), N = .N), by = c('species')] %>%
           setkey('species')
cat('\nOverall, ', sporpos['An. gambiae', Pfal.rate]*100, '% of An. gambiae and ',
    sporpos['An. funestus', Pfal.rate]*100, '% of An. funestus were positive for ',
    'P. falciparum respectively.\n\n', sep = '')

sporpos.byround <- ento.data[, .(Pfal.total = sum(Pfal),
                                 Pfal.rate = round(sum(Pfal) / .N, 4),
                                 Pvom.total = sum(OVM), 
                                 Povm.rate = round(sum(OVM) / .N, 4),
								 N = .N), 
                               by = c('species', 'RND')] %>%
                   setkey('species')
cat('\nSporozoite positivity split by round:\n')
print(sporpos.byround)

# Since there are no sporozoite positive samples in arabiensis, we can exclude them from some of the plots
ento.data.fg <- ento.data[MolSpecies %in% c('FUN', 'S'), ]

##################################################
split.plasmo.counts <- function(formula, ento.data){
  Fal <- dcast(ento.data, formula, length)
  Fal$N <- Fal[['0']]+Fal[['1']]
  Fal$prop <- (Fal[['1']]/(Fal[['0']]+Fal[['1']]))
  Fal$Perc <- 100*Fal$prop
  prt <- mapply(prop.test, Fal[['1']], Fal[['0']] + Fal[['1']])['conf.int', ]
  Fal$LLim <- 100*sapply(prt, function(x) x[1])
  Fal$ULim <- 100*sapply(prt, function(x) x[2])
#  Fal$ULim <- 100*(Fal$prop + 1.96 * sqrt(Fal$prop * (1-Fal$prop)/Fal$N))
#  Fal$LLim <- 100*(Fal$prop - 1.96 * sqrt(Fal$prop * (1-Fal$prop)/Fal$N))
  #Fal$LLim[Fal$LLim<0] <- 0
  Fal
}

frequency.line.plot <- function(Plasmo, ylab = 'P.falciparum infection prevalence (%)'){
  facet.terms <- setdiff(colnames(Plasmo), c('species', 'RND', '0', '1', 'N', 'prop', 'Perc', 'ULim', 'LLim'))
  facet.formula <- as.formula(paste('~', paste(facet.terms, collapse = ' + ')))
  flp <- ggplot(Plasmo, aes(x=RND, y=Perc, fill = species, color = species)) + 
  geom_line(size = 1.5, position=position_dodge(width=c(0.4)), show.legend = F) +
  geom_point(color = 'grey80', pch = 21, size = 3, position=position_dodge(width=c(0.4))) +
  geom_errorbar(aes(ymin=LLim, ymax=ULim,col=species),width=0.25,cex=0.5, 
                position=position_dodge(width=c(0.4)), show.legend = F)+
  labs(fill="Species") +
  facet_wrap(facet.formula, scales= "free_y") +
  ylab(ylab) +
  xlab('Round') +
  scale_fill_brewer(palette = 'Dark2') +
  scale_color_brewer(palette = 'Dark2') +
  theme_bw() +
  theme(legend.text = element_text(size=15), 
		legend.title = element_text(size = 15),
        strip.text.x = element_text(size = 15, face = 'bold'),
        axis.text = element_text(size=15),
        axis.title = element_text(size=15,face='bold'))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())
  flp
}

Fal.Simp <- split.plasmo.counts(species + RND ~ Pfal, ento.data.fg)
Fal.Loc <- split.plasmo.counts(LLIN.actual + Location + species + RND ~ Pfal, ento.data.fg)
flp.Fal.Loc <- frequency.line.plot(Fal.Loc)
print(flp.Fal.Loc)
ggsave('infection_numbers_location_Fal.png', flp.Fal.Loc, width = 9, height = 7)
Fal.Arm <- split.plasmo.counts(LLIN.actual + species + RND ~ Pfal, ento.data.fg)
flp.Fal.Arm <- frequency.line.plot(Fal.Arm)
print(flp.Fal.Arm)
ggsave('infection_numbers_arm_Fal.png', flp.Fal.Arm, width = 9, height = 7)

###############
# Also for OVM plasmodia

OVM.Arm <- split.plasmo.counts(LLIN.actual+species+RND~OVM, ento.data.fg)
flp.OVM.Arm <- frequency.line.plot(OVM.Arm, ylab = 'P. ovale-vivax-malariae infection prevalence (%)')
print(flp.OVM.Arm)
ggsave('infection_numbers_arm_OVM.png', flp.OVM.Arm, width = 9, height = 7)
###########################################################

library("cowplot")
Plasmodium <- plot_grid(flp.Fal.Arm, flp.OVM.Arm, nrow=2)

pdf("Figure 2_ Lynd et al Plasmodium_11Apr2023.pdf",width=16,height=10)
Plasmodium
dev.off()
# Will remove redundant legend in inkscape

######################
# To produce table 2 in Lynd et al 2023
# Drop arabiensis as always negative
data.Ga.Fun <-  subset(ento.data.afg, MolSpecies!="AR")
data.Ga <-  subset(ento.data.afg, MolSpecies=="S")
data.Fun <-  subset(ento.data.afg, MolSpecies=="FUN")
##########################################

cat('\n  ##############################################')
cat('\n  # Model1: Pfal ~ RND + LLIN.actual + (1|HSD) #')
cat('\n  ##############################################\n')

cat('\n\tgambiae + funestus:\n')
model.1.All <- glmmTMB(Pfal ~ RND + LLIN.actual + species + Location + (1|HSD), 
                       data = data.Ga.Fun,
                       family = 'binomial')
print(drop1(model.1.All, test = 'Chisq'))
print(fixef(model.1.All))
######################
cat('\n\tgambiae alone:\n')
model.1.Ga <- glmmTMB(Pfal ~ RND + LLIN.actual + Location + (1|HSD), 
                      data = data.Ga,
                      family = 'binomial')
print(drop1(model.1.Ga, test = 'Chisq'))
print(fixef(model.1.Ga))
##############################
cat('\n\tfunestus alone:\n')
model.1.Fun <- glmmTMB(Pfal ~ RND + LLIN.actual + Location + (1|HSD), 
                       data = data.Fun,
                       family = 'binomial')
print(drop1(model.1.Fun, test = 'Chisq'))
print(fixef(model.1.Fun))

cat('\n  #############################################')
cat('\n  # Model2: OVM ~ RND + LLIN.actual + (1|HSD) #')
cat('\n  #############################################\n')

cat('\n\tgambiae + funestus:\n')
model.2.All <- glmmTMB(OVM ~ RND + LLIN.actual + species + Location + (1|HSD), 
                       data=data.Ga.Fun,
                       family = 'binomial')
print(drop1(model.2.All, test = 'Chisq'))
print(fixef(model.2.All))
######################
cat('\n\tgambiae alone:\n')
model.2.Ga <- glmmTMB(OVM ~ RND + LLIN.actual + Location + (1|HSD), 
                      data=data.Ga,
                      family = 'binomial')
print(drop1(model.2.Ga, test = 'Chisq'))
print(fixef(model.2.Ga))
##############################
cat('\n\tfunestus alone:\n')
model.2.Fun <- glmmTMB(OVM ~ RND + LLIN.actual + Location + (1|HSD), 
                       data=data.Fun,
                       family = 'binomial')
print(drop1(model.2.Fun, test = 'Chisq'))
print(fixef(model.2.Fun))

cat('\n  ##############################################')
cat('\n  # Model3: Pfal ~ RND * LLIN.actual + (1|HSD) #')
cat('\n  ##############################################\n')

cat('\n\tgambiae + funestus:\n')
model.3.All <- glmmTMB(Pfal ~ RND * LLIN.actual + species + Location + (1|HSD), 
                       data=data.Ga.Fun,
                       family = 'binomial')
print(drop1(model.3.All, test = 'Chisq'))
print(fixef(model.3.All))
######################
cat('\n\tgambiae alone:\n')
model.3.Ga <- glmmTMB(Pfal ~ RND * LLIN.actual + Location + (1|HSD), 
                      data=data.Ga,
                      family = 'binomial')
print(drop1(model.3.Ga, test = 'Chisq'))
print(fixef(model.3.Ga))
##############################
cat('\n\tfunestus alone:\n')
model.3.Fun <- glmmTMB(Pfal ~ RND * LLIN.actual + Location + (1|HSD), 
                       data=data.Fun,
                       family = 'binomial')
print(drop1(model.3.Fun, test = 'Chisq'))
print(fixef(model.3.Fun))

cat('\n  #############################################')
cat('\n  # Model4: OVM ~ RND * LLIN.actual + (1|HSD) #')
cat('\n  #############################################\n')

cat('\n\tgambiae + funestus:\n')
model.4.All <- glmmTMB(OVM ~ RND * LLIN.actual + species + Location + (1|HSD), 
                       data=data.Ga.Fun,
                       family = 'binomial')
print(drop1(model.4.All, test = 'Chisq'))
print(fixef(model.4.All))
######################
cat('\n\tgambiae alone:\n')
model.4.Ga <- glmmTMB(OVM ~ RND * LLIN.actual + Location + (1|HSD), 
                      data=data.Ga,
                      family = 'binomial')
print(drop1(model.4.Ga, test = 'Chisq'))
print(fixef(model.4.Ga))
##############################
cat('\n\tfunestus alone:\n')
model.4.Fun <- glmmTMB(OVM ~ RND * LLIN.actual + Location + (1|HSD), 
                       data=data.Fun,
                       family = 'binomial')
print(drop1(model.4.Fun, test = 'Chisq'))
print(fixef(model.4.Fun))
########################################################################

#######################################################################
# Table 3. Associations between resistance markers and infection status
#######################################################################
Inf.Gam <- ento.data[MolSpecies=="S"]
# Recode the 2La 
Inf.Gam[, X2La := gsub('I', 'a', gsub('S', '+', unlist(X2La)))]

cat('\n  #####################################################################')
cat('\n  # Fisher tests of association of P. falciparum with genetic markers #')
cat('\n  #####################################################################\n')

cat('\n\tCyp6p4\n')
print(with(Inf.Gam, table(Cyp6p4,Pfal)))
print(with(Inf.Gam, fisher.test(Cyp6p4, Pfal)))

cat('\n\tCoeae1d\n')
print(with(Inf.Gam, table(Coeae1d,Pfal)))
print(with(Inf.Gam, fisher.test(Coeae1d,Pfal)))

cat('\n\tCyp4j5\n')
print(with(Inf.Gam, table(Cyp4j5,Pfal)))
print(with(Inf.Gam, fisher.test(Cyp4j5,Pfal)))

cat('\n\tkdr\n')
print(with(Inf.Gam, table(kdr,Pfal)))
print(with(Inf.Gam, fisher.test(kdr,Pfal)))

cat('\n\t2La\n')
print(with(Inf.Gam, table(X2La,Pfal)))
print(with(Inf.Gam, fisher.test(X2La,Pfal)))

cat('\n  #######################################################################')
cat('\n  # Fisher tests of association of other Plasmodia with genetic markers #')
cat('\n  #######################################################################\n')

cat('\n\tCyp6p4\n')
print(with(Inf.Gam, table(Cyp6p4,OVM)))
print(with(Inf.Gam, fisher.test(Cyp6p4, OVM)))

cat('\n\tCoeae1d\n')
print(with(Inf.Gam, table(Coeae1d,OVM)))
print(with(Inf.Gam, fisher.test(Coeae1d,OVM)))

cat('\n\tCyp4j5\n')
print(with(Inf.Gam, table(Cyp4j5,OVM)))
print(with(Inf.Gam, fisher.test(Cyp4j5,OVM)))

cat('\n\tkdr\n')
print(with(Inf.Gam, table(kdr,OVM)))
print(with(Inf.Gam, fisher.test(kdr,OVM)))

cat('\n\t2La\n')
print(with(Inf.Gam, table(X2La,OVM)))
print(with(Inf.Gam, fisher.test(X2La,OVM)))
####################

# Make a figure. 
change.col <- function(color, light.change){
	col.rgb <- col2rgb(color)/255
	if (light.change < 1){
		red <- 1-(1-col.rgb[1])*light.change
		blue <- 1-(1-col.rgb[2])*light.change
		green <- 1-(1-col.rgb[3])*light.change
	}
	else {
		red <- col.rgb[1]*(2-light.change)
		blue <- col.rgb[2]*(2-light.change)
		green <- col.rgb[3]*(2-light.change)
	}
	rgb(cbind(red, blue, green))
}

plasmo.table <- function(marker, infection.table, allele.order = NULL,
						 include.plot = T, cel.col = 'cadetblue', title = marker, cex = 0.8,
						 include.pvalue = T){
	infection.rates <- infection.table[, .(N = .N, 
								         Pf.num = sum(Pfal),
								         Pf.rate = round(sum(Pfal)/.N, 3),
								         Povm.num = sum(OVM),
								         Povm.rate = round(sum(OVM)/.N, 3)),
                                       by = marker] %>%
	                   .[!is.na(get(marker))]
	if (!is.null(allele.order)){
		infection.rates$order <- setNames(1:length(allele.order), allele.order)[infection.rates[[marker]]]
		infection.rates <- infection.rates[order(order)][, order := NULL]
	}
	if (include.plot){
		table.for.plotting <- infection.rates[, .(genotype = get(marker),
												  N = N,
												  'P. falciparum' = paste(Pf.num, ' (', Pf.rate*100, '%)', sep = ''),
		                                          'P. ovm' = paste(Povm.num, ' (', Povm.rate*100, '%)', sep = '')),
										      ]
		par(mar = c(1, 0.1, 3, 0.1), cex = cex)
		plot(c(0,ncol(table.for.plotting)), -c(0,nrow(table.for.plotting)), 
			 type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', 
			 xlab = '', ylab = '')
		dark.col <- change.col(cel.col, 1.4)
		title(title, line = 1.5, col.main = dark.col)
		x.matrix <- matrix(1:ncol(table.for.plotting), nrow(table.for.plotting), ncol(table.for.plotting), byrow = T)
		y.matrix <- matrix(1:nrow(table.for.plotting), nrow(table.for.plotting), ncol(table.for.plotting))
		rect(x.matrix - 1, 1 - y.matrix, x.matrix, - y.matrix, col = cel.col, border = 'white')
		text(x.matrix - 0.5, 0.5 - y.matrix, as.matrix(table.for.plotting), col = 'white', font = 2)
		text(x.matrix[1,] - 0.5, (par('usr')[4] - par('usr')[3])/100, adj = c(0.5,-0.4), col = dark.col, colnames(table.for.plotting), xpd = NA, font = 2)
		if (include.pvalue){
			pfal.fisher <- signif(fisher.test(infection.rates[, .(Pf.num, N-Pf.num)])$p.value, 2)
			povm.fisher <- signif(fisher.test(infection.rates[, .(Povm.num, N-Povm.num)])$p.value, 2)
			mtext('Fisher test P =', 1, line = 0, at = 2, adj = 1.05, col = dark.col, font = 2, cex = cex)
			mtext(pfal.fisher, 1, line = 0, at = 2.5, adj = 0.5, col = dark.col, font = 2, cex = cex)
			mtext(povm.fisher, 1, line = 0, at = 3.5, adj = 0.5, col = dark.col, font = 2, cex = cex)
		}
	}
	invisible(infection.rates)
}

svg('kdr_plasmo_assoc.svg', width = 4, height = 3)
plasmo.table('kdr', Inf.Gam, title = 'Vgsc-995', allele.order = c('LL', 'LS', 'FL', 'FS', 'SS', 'FF'))
dev.off()
#
svg('Cyp6p4_plasmo_assoc.svg', width = 4, height = 1.75)
plasmo.table('Cyp6p4', Inf.Gam, cel.col = 'indianred1', title = 'Cyp6p4-236M', allele.order = c('II', 'IM', 'MM'))
dev.off()
#
svg('Cyp4j5_plasmo_assoc.svg', width = 4, height = 1.75)
plasmo.table('Cyp4j5', Inf.Gam, cel.col = 'thistle4', title = 'Cyp4j5-43F', allele.order = c('LL', 'LF', 'FF'))
dev.off()
#
svg('Coeae1d_plasmo_assoc.svg', width = 4, height = 1.75)
plasmo.table('Coeae1d', Inf.Gam, cel.col = 'palegreen3', allele.order = c('SS', 'SR', 'RR'))
dev.off()
#
svg('2La_plasmo_assoc.svg', width = 4, height = 1.75)
plasmo.table('X2La', Inf.Gam, cel.col = 'gold3', title = '2La inversion', allele.order = c('++', '+a', 'aa'))
dev.off()

# Output a table summary of the P-values of interest
p.values.of.interest <- c(model1_Gam = drop1(model.1.Ga, test = 'Chisq')['RND', 'Pr(>Chi)'],
                          model1_Fun = drop1(model.1.Fun, test = 'Chisq')['RND', 'Pr(>Chi)'],
                          model2_Gam = drop1(model.2.Ga, test = 'Chisq')['RND', 'Pr(>Chi)'],
                          model2_Fun = drop1(model.2.Fun, test = 'Chisq')['RND', 'Pr(>Chi)'])

p.value.table <- data.table(model = names(p.values.of.interest),
                            P = p.values.of.interest,
                            FDR = p.adjust(p.values.of.interest, method = 'BH'))

fwrite(p.value.table, 'sporozoite_P_values_summary.csv', sep = '\t')

