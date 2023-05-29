library(data.table)
library(ggplot2)
library(stringr)

# Load the data
all.data <- fread('PBO Ento ALL molecular data V5 31-03-23 merged.csv', na.strings = '')

# Give the net types their actual names
net.code.conversion <- c('PermaNet-2', 'PermaNet-3', 'Olyset', 'OlysetPlus')
all.data$Net.intend <- net.code.conversion[all.data$Net.intend]
all.data$Net.actual <- net.code.conversion[all.data$Net.actual]
llin.code.conversion <- c('Non-PBO LLIN', 'PBO LLIN')
all.data$LLIN.intend <- llin.code.conversion[all.data$LLIN.intend]
all.data$LLIN.actual <- llin.code.conversion[all.data$LLIN.actual]
all.data[, wgs.sample.id := paste(HHID, str_pad(ID, 3, pad = '0'), sep = '')]

# Fix the species names
species.conversion <- c('AO'  = 'An. sp',
                        'AR'  = 'An. arabiensis',
                        'FUN' = 'An. funestus',
                        'LEE' = 'An. leesoni',
                        'M'   = 'An. coluzzii',
                        'MOU' = 'An. moucheti',
                        'PAR' = 'An. parensis',
                        'RIV' = 'An. rivulorum',
                        'S'   = 'An. gambiae')
all.data[, species := ..species.conversion[MolSpecies]]

# Fix some of the gene names
gene.conversion <- c('Cyp4J5' = 'Cyp4j5',
                     'Coeaed1' = 'Coeae1d',
                     'Cyp6P4' = 'Cyp6p4',
                     'X6AA' = 'Cyp6aa',
                     'X6AA_Dup' = 'Cyp6aa_Dup')
setnames(all.data, names(gene.conversion), gene.conversion)

# Remove rows with no clear LLIN or species data
all.data <- all.data[!is.na(LLIN.actual), ] %>%
            .[!is.na(species), ]
fwrite(all.data, file="ento_geno_plasmo_data.csv", sep = ',')

# Plot sample sizes
num.mozzies <- all.data[MolSpecies %in% c('AR', 'FUN', 'S'), 
                        aggregate(wgs.sample.id, 
                                  by = list(species = species, 
                                            LLIN.actual = LLIN.actual, 
                                            Location = Location, 
                                            RND = RND), 
                                  length, 
                                  drop = F
                        )
]

frequency.line.plot <- ggplot(num.mozzies, aes(x=RND, y=x, fill = species, color = species)) + 
                       geom_line(size = 1.2) +
                       geom_point(color = 'grey80', pch = 21, size = 2) +
                       geom_text(aes(label = x, 
                                     y = c('An. arabiensis' = 650, 
                                           'An. funestus' = 700,
                                           'An. gambiae' = 750)[species],
                                     fontface = 'bold'),
                                 show.legend = F) +
                       facet_wrap(~ Location + LLIN.actual) +
			           ylab('Number of mosquitoes') +
			           scale_fill_brewer(palette = 'Dark2') +
			           scale_color_brewer(palette = 'Dark2') +
			           theme_bw() +
                       theme(panel.grid.major.x = element_blank(),
                             panel.grid.minor.x = element_blank(),
                             panel.grid.minor.y = element_blank())
#print(frequency.line.plot)
ggsave('sample_numbers.png', frequency.line.plot, width = 9, height = 7)

# Sample sizes for manuscript Table 1:
sample.sizes <- data.table(num.mozzies)[, sum(x), by = c('species', 'RND')] %>%
                dcast(species ~ RND) %>%
                setkey('species') %>%
				.[c('An. gambiae', 'An. arabiensis', 'An. funestus')]

cat('\nSample sizes of three main species:\n')
print(sample.sizes)

other.mozzies <- all.data[MolSpecies %in% c('M', 'LEE', 'PAR', 'RIV', 'MOU'), 
                          aggregate(wgs.sample.id, 
                                    by = list(RND = RND), 
                                    length, 
                                    drop = F
                        )
]
cat('\nSample sizes of combined coluzzii, parensis, leesoni, rivulorum and mouchetti:\n')
print(other.mozzies)


