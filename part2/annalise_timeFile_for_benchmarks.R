

infile <- read.csv("timeFile.csv", header = FALSE, stringsAsFactors = FALSE)

library(reshape2)
library(ggplot2)
library(gridExtra)


colnames(infile) <- c("N_in", "process", "time")
#connvert N_in to actual values
infile$N_in <- as.integer(substr(infile$N_in, 19, length(infile$N_in)-14))

##blast hits vs queries
#create a new df with the N_in and N_blst_hits 
blastHits <- infile[infile$process=="blastHits",][,-2]
colnames(blastHits) <- c("N_in", "N_blst")


#creat on theme for ggplot to use thoughout
phyloPiTheme <- theme_bw(base_size = 14) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


#plot and fit N_blst against N_in, force through 0
fit <-  lm(blastHits$N_blst ~ blastHits$N_in-1)
coefficients(fit)
summary(fit)

nBlstPlot <- ggplot(blastHits, aes(x = N_in, y = N_blst)) +
  geom_point(shape = 1) + 
  geom_smooth(method=lm, se = FALSE, colour='black',formula=y~x-1, size = 0.25) + 
  labs(x="Number of input sequences", y="Number of sequences retrieved by BLAST") +
  annotate("text", x=41, y=72, label = "y == 4.628 * x", parse=T)+
  annotate("text", x=40, y=60, label = "R^2 == 0.998", parse=T) +
  #annotate("label", x = 0, y = 250)
  phyloPiTheme
nBlstPlot
ggsave(file="blast hits vs queries.svg", plot=nBlstPlot)



#parse the dataframe for blast
blast <- infile[infile$process=="blast",][-2]
blast <- cbind(blastHits, blast)[-3]

fitBlastT <-  lm(blast$time ~ blast$N_in-1)
coefficients(fitBlastT)
summary(fitBlastT)

blastT <- ggplot(blast, aes(x = N_in, y = time)) + 
  geom_point(shape = 1) +
  geom_smooth(method=lm, se = FALSE, colour='black',formula=y~x-1, size = 0.25) + 
  labs(x="Number of input sequences", y="time (s)") +
  annotate("text", x=40, y=200, label = "y == 11.02 * x", parse=T)+
  annotate("text", x=40, y=170, label = "R^2 == 1", parse=T) +
  phyloPiTheme
blastT
ggsave(file="blast time.svg", plot = blastT)



#parse the dataframe for mafft
mafft <- infile[infile$process=="mafftTime",][-2]
mafft <- cbind(blastHits, mafft)[-3]

#plot and fit time against N-blst, mafft
t <- mafft$time
N <- mafft$N_blst
fit <- nls(t~b * N**a, start = list(a=2,b=0.1))
cor(t,predict(fit))

summary(fit)
fit2 <- nls(t~b * N**2, start = list(b=0.1))
summary(fit2)
anova(fit, fit2)

mafftT <- ggplot(mafft, aes(x=N_blst, y=time))+
  geom_point(shape = 1) +
  geom_smooth(method="nls",
              formula = y ~ b * x**a, 
              method.args=list(start = c(a = 2, b = 0.1)), 
              se = FALSE, colour='black', size = 0.25) +
  labs(x="Number of sequences in alignment", y="time (s)") +
  annotate("text", x=190, y=1800, label = "y == 0.153 * x^1.92", parse=T)+
  phyloPiTheme 
  
mafftT
ggsave(file="mafftTime.svg", plot=mafftT)



y <- mafft$time
x <- mafft$N_blst
fitMafft <- nls(y ~ b * x**a, start = list(a = 2, b = 0.1))
plot(x=mafft$N_blst, y=mafft$time)
lines(x, predict(fitMafft))
plot(x, resid(fitMafft))

hist(resid(fitMafft), breaks = 25)
abline(h=0)
qqnorm(resid(fitMafft))
qqline(resid(fitMafft))



#plot for fastTree
#for number of queries in input
fastTreeNq <- infile[infile$process == "fasttreeTime",][-2]
#add blastHits here
fastTreeWithBlstHits <- cbind(fastTreeNq,blastHits)[-3]



##an of time vs input query
fitFastTreeNq <-  lm(fastTreeNq$time ~ blast$N_in-1)
coefficients(fitFastTreeNq)
summary(fitFastTreeNq)

fastTreeNqP <- ggplot(fastTreeNq, aes(x = N_in, y = time)) + 
  geom_point(shape = 1) +
  geom_smooth(method=lm, se = FALSE, colour='black',formula=y~x-1, size = 0.25) + 
  labs(x="Number of query sequences", y="time (s)") +
  annotate("text", x=40, y=70, label = "y == 3.02 * x", parse=T)+
  annotate("text", x=40, y=60, label = "R^2 == .994", parse=T)+
  phyloPiTheme
fastTreeNqP
ggsave(file="fastTree query sequences.svg", plot = fastTreeNqP)

##an of time vs blast results, thus number of samples in analysis
fitFastTreeWithBlstHits <-  lm(fastTreeWithBlstHits$time ~ fastTreeWithBlstHits$N_blst -1)
coefficients(fitFastTreeWithBlstHits)
summary(fitFastTreeWithBlstHits)

fastTreeWithBlastHitsP <- ggplot(fastTreeWithBlstHits, aes(x = N_blst, y = time)) + 
  geom_point(shape = 1) +
  geom_smooth(method=lm, se = FALSE, colour='black',formula=y~x-1, size = 0.25) + 
  labs(x="Number of sequences in alingment", y="time (s)") +
  annotate("text", x=170, y=60, label = "y == 0.652 * x", parse=T)+
  annotate("text", x=170, y=50, label = "R^2 == .993", parse=T)+
  phyloPiTheme
fastTreeWithBlastHitsP
ggsave(file="fastTree total sequences.svg", plot = fastTreeWithBlastHitsP)



plots = list(nBlstPlot + ggtitle("(a)"),
             blastT + ggtitle("(b)"), 
             mafftT + ggtitle("(c)"), 
             fastTreeWithBlastHitsP + ggtitle("(d)"))
figure <- grid.arrange(grobs= lapply(plots, "+", theme(plot.margin=margin(10,25,10,25))))
figure

ggsave(file = "benchFig.svg",plot = figure, width = 7.66, height = 10)
ggsave(file = "benchFig.png",plot = figure, width = 7.66, height = 10)


