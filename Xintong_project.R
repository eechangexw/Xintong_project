#Calculate Acoustic Index(A total of 27 files, repeat the following steps for each file)
setwd("E:/RBird_001")
files <- list.files(path = getwd(), pattern = "wav$", ignore.case = T )
df <- data.frame(matrix(ncol = length(names(z)), nrow = 0))
colnames(df) <- names(z)
results <- list() # Storage for all acoustic index results in WAV files
for (file in files) {
  wav <- readWave(file, from = 0, to = 60, units = "seconds", header = FALSE, toWaveMC = NULL)
  full_path <- file.path(getwd(), file)
  
  # calculate
  x <- soundscapespec(wav, plot = TRUE, main = file) 
  v <- as.vector(x[, 2]) 
  f <-48000
  
  AEI <- acoustic_evenness(wav, max_freq = 10000)
  AEI.L <- AEI$aei_left
  ADI <- acoustic_diversity(wav, max_freq = 10000) 
  ADI.L <- ADI$adi_left
  ACI <- acoustic_complexity(wav, min_freq = 2000, max_freq = 11000)
  ACI.L <- ACI$AciTotAll_left
  BIO <- bioacoustic_index(wav, min_freq = 2000, max_freq = 8000)
  BIO.L <-BIO$left_area
  NDSI <- ndsi(wav)
  NDSI.L<-(NDSI$ndsi_left)
  TE <- H(wav)
  envorni <- env(wav, f = 22050, plot = FALSE)
  Ht <- th(envorni)
  speca <- spec(wav, f = f, plot = FALSE) 
  Hf <- sh(speca) 
  MAE <- M(wav)
  spec <- meanspec(wav, plot = FALSE) 
  peaks <- fpeaks(spec, plot = FALSE)
  NP <- length(peaks) / 2
  
  # add result into the list
  finf <- data.table(file.info(dir(getwd()), extra_cols = F))
  finft <- data.table(finf$mtime[file = file])
  z <- list(AEI = AEI.L, ADI = ADI.L, ACI = ACI.L, BIO = BIO.L, NDSI = NDSI.L, 
            TE = TE, Ht = Ht, Hf = Hf,
            MAE = MAE, NP = NP, DateTime = finft)
  df <- rbind(df, data.frame(z, row.names = make.names(rep(file, length(z[[1]])), unique = TRUE)))
}

library(openxlsx)
write.xlsx(df, file = "Bird001.xlsx", rowNames = TRUE)

#Then, I integrated 27 tables into a large table and named "BirdQ"



#PCA analysis
setwd("C:/Users/xinto/Desktop/Data output")
library(vegan)
data <- read.csv("BirdQ.csv", header=TRUE, sep=",")
acoustic_data <- data[, c("AEI", "ADI", "ACI", "BIO", "NDSI", "TE", "Ht", "Hf", "MAE", "NP")]
acoustic_data_scaled <- scale(acoustic_data)

# pca
pca_result <- rda(acoustic_data)
pca <- rda(acoustic_data_scaled)

# read scores
scores <- scores(pca)
loadings <- scores(pca)$species
# result
df <- summary(pca_result)
library(openxlsx)
pca_result <- rda(acoustic_data_scaled)
pca_scores <- scores(pca_result) 
pca_scores_df <- as.data.frame(pca_scores$sites)  
# write.xlxs
write.xlsx(pca_scores_df, file = "pca_scores.xlsx", rowNames = FALSE)

# Then I added pc1&pc2 value into BirdQ.csv




#GAMM model
setwd("C:/Users/xinto/Desktop/Data output")
data <- read.csv('BirdQ.csv')
library(mgcv)

# Basic
Forest_DIS <- as.factor(data$Forest_DIS)
formula <- PC1 ~ s(Time) + s(Time, Forest_DIS) 
model1 <- gamm(formula, data = data, random = list(Sample_point = ~1),
               method = "REML", family = gaussian(link = "identity"))
summary(model1$gam)
summary(model1$lme)
formula2 <- PC2 ~ s(Time) + s(Time, Forest_DIS) 
model2 <- gamm(formula2, data = data, random = list(Sample_point = ~1),
               method = "REML", K=K, family = gaussian(link = "identity"))
summary(model2$gam)
summary(model2$lme)
#Best K
aic_list <- list()
for (k in 1:10) {
  model1 <- gamm(formula , data = data,random = list(Sample_point = ~1),
                 method = "REML", family = gaussian(link = "identity"))
  aic <- AIC(gamm_model$lme)
  aic_list[[as.character(k)]] <- aic
}

min_aic_model <- which.min(unlist(aic_list))

#Date
Date <- data$Date
model1 <- gamm(data$PC2 ~ s(data$Time) + s(data$Time, data$Forest_DIS) + data$Date, 
               data = data, random = list(Sample_point = ~1),
               method = "REML", family = gaussian(link = "identity"))
summary(model1$gam)
AIC(model1$lme)
summary(model1$lme) 

# test spatial correlation
library(gstat) 
library(geoR) 
library(sp)
coordinates(data) <- ~ Long + Lat
distances <- dist(coordinates(data))

# Create a variogram for PC1
variogram_pc1 <- variogram(PC1 ~ 1, data = data, distances = distances)

# Create a variogram for PC2
variogram_pc2 <- variogram(PC2 ~ 1, data = data, distances = distances)

# Plot the variograms
plot(variogram_pc1, main = "Variogram for PC1")
plot(variogram_pc2, main = "Variogram for PC2")

# ar1
model_ar1 <- gamm(formula, data = data, correlation = corAR1(form = ~ 1 | Time + Sample_point),
                  random = list(Sample_point = ~1))
AIC(model_ar1$lme)
model_ar1$gam
summary(model_ar1$gam)
model2_ar1 <- gamm(formula2, data = data, correlation = corAR1(form = ~ 1 | Time + Sample_point), 
                   random = list(Sample_point = ~1))
model2_ar1$lme
AIC(model2_ar1$lme)
model2_ar1$gam
summary(model2_ar1$gam)
#ARMA
model_corARMA <- gamm(formula , data = data, 
                      correlation = corARMA(p=2, q=2, form = ~ 1 | Time + Sample_point), 
                      random = list(Sample_point = ~1))
model_corARMA$lme
AIC(model_corARMA$lme)
summary(model_corARMA$gam)
model2_corARMA <- gamm(formula2, data = data, 
                       correlation = corARMA(p=2, q=2, form = ~ 1 | Time + Sample_point), 
                       random = list(Sample_point = ~1))
model2_corARMA$lme
AIC(model2_corARMA$lme)
summary(model2_corARMA$gam)
#corCAR
model_corCAR <- gamm(formula , data = data, 
                     correlation = corCAR1(form = ~ 1 | Time + Sample_point), 
                     random = list(Sample_point = ~1))
model_corCAR$lme
AIC(model_corCAR$lme)
summary(model_corCAR$gam)
model2_corCAR <- gamm(formula2, data = data, 
                      correlation = corCAR1(form = ~ 1 | Time + Sample_point), 
                      random = list(Sample_point = ~1))
model2_corCAR$lme
model2_corCAR$gam
summary(model_corcar2$gam)

#Final
model_Final <- gamm(data$PC1 ~ s(data$Time) + s(data$Time, data$Forest_DIS) + data$Date , 
                    data = data, correlation = corARMA(p=2, q=2, form = ~ 1 | Time + Sample_point), 
                    random = list(Sample_point = ~1))
model_Final$lme
AIC(model_Final$lme)
summary(model_Final$gam)
model2_Final <- gamm(data$PC2 ~ s(data$Time) + s(data$Time, data$Forest_DIS) + data$Date , 
                     data = data, correlation = corARMA(p=2, q=2, form = ~ 1 | Time + Sample_point), 
                     random = list(Sample_point = ~1))
model2_Final$lme
AIC(model2_Final$lme)
summary(model2_Final$gam)

par(mfrow = c(1, 1)) 
# Temporal correlations
residuals <- resid(model1$gam)
acf_values <- acf(residuals, plot = FALSE)
plot(acf_values, main = "PC1 GAMM before")
residuals_Final <- resid(model_Final$gam)
acf_values_Final <- acf(residuals_Final, plot = FALSE)
plot(acf_values_Final, main = "PC1 GAMM after")
residuals2 <- resid(model2$gam)
acf_values2 <- acf(residuals2, plot = FALSE)
plot(acf_values2, main = "PC2 GAMM before")
residuals2_Final <- resid(model2_Final$gam)
acf_values2_Final <- acf(residuals2_Final,  plot = FALSE)
plot(acf_values2_Final, main = "PC2 GAMM after")

#smooth curve
a_data <- subset(data, Distance_Cata == "1")
b_data <- subset(data, Distance_Cata == "2")
c_data <- subset(data, Distance_Cata == "3")
d_data <- subset(data, Distance_Cata == "4")
pc1_a_data <- gamm(PC1 ~ s(Time) + Date, data = a_data, correlation = corARMA(p=2, q=2, form = ~ 1 | Time), 
                   random = list(Sample_point = ~1))
pc2_a_data <- gamm(PC2 ~ s(Time) + Date, data = a_data, correlation = corARMA(p=2, q=2, form = ~ 1 | Time), 
                   random = list(Sample_point = ~1))
pc1_b_data <- gamm(PC1 ~ s(Time) + Date, data = b_data, correlation = corARMA(p=2, q=2, form = ~ 1 | Time), 
                   random = list(Sample_point = ~1))
pc2_b_data <- gamm(PC2 ~ s(Time) + Date, data = b_data, correlation = corARMA(p=2, q=2, form = ~ 1 | Time), 
                   random = list(Sample_point = ~1))
pc1_c_data <- gamm(PC1 ~ s(Time) + Date, data = c_data, correlation = corARMA(p=2, q=2, form = ~ 1 | Time), 
                   random = list(Sample_point = ~1))
pc2_c_data <- gamm(PC2 ~ s(Time) + Date, data = c_data, correlation = corARMA(p=2, q=2, form = ~ 1 | Time), 
                   random = list(Sample_point = ~1))
pc1_d_data <- gamm(PC1 ~ s(Time) + Date, data = d_data, correlation = corARMA(p=2, q=2, form = ~ 1 | Time), 
                   random = list(Sample_point = ~1))
pc2_d_data <- gamm(PC2 ~ s(Time) + Date, data = d_data, correlation = corARMA(p=2, q=2, form = ~ 1 | Time), 
                   random = list(Sample_point = ~1))
summary(pc1_a_data$lme)
summary(pc1_a_data$gam)
summary(pc1_b_data$lme)
summary(pc1_b_data$gam)
summary(pc1_c_data$lme)
summary(pc1_c_data$gam)
summary(pc1_d_data$lme)
summary(pc1_d_data$gam)
summary(pc2_a_data$lme)
summary(pc2_a_data$gam)
summary(pc2_b_data$lme)
summary(pc2_b_data$gam)
summary(pc2_c_data$lme)
summary(pc2_c_data$gam)
summary(pc2_d_data$lme)
summary(pc2_d_data$gam)
library(ggplot2)




plot1_a_data<- ggplot(a_data, aes(x = Time, y = -1 * PC1)) +
  geom_smooth(method = "gam", se = TRUE, color = "black", size=1.5)  +
  geom_point(size = 1) +
  labs(title = "In the forest", x = "Time(h)", y = "PC1") +
  annotate("text", x = max(23), y = max(4), label = paste("Edf=", 7.87, ", F=", 20.3, ", P<", 0.001), hjust = 1, vjust = 0.5, size = 8) +
  theme_minimal()+
  theme(axis.line = element_line(color = "black"),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 29),
        panel.grid = element_blank())
plot1_a_data
plot1_b_data<- ggplot(b_data, aes(x = Time, y = -1 * PC1)) +
  geom_smooth(method = "gam", se = TRUE, color = "black", size = 1.5)  +
  geom_point(size = 1) +
  labs(title = "0-100m", x = "Time(h)", y = "PC1") +
  annotate("text", x = max(23), y = max(4), label = paste("Edf=", 7.35, ", F=", 14.29, ", P<", 0.001), hjust = 1, vjust = 0.5, size = 8) +
  theme_minimal()+
  theme(axis.line = element_line(color = "black"),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 29),
        panel.grid = element_blank())
plot1_b_data
plot1_c_data<- ggplot(c_data, aes(x = Time, y = -1 * PC1)) +
  geom_smooth(method = "gam", se = TRUE, color = "black", size = 1.5)  +
  geom_point(size = 1) +
  labs(title = "100-250m", x = "Time(h)", y = "PC1") +
  annotate("text", x = max(23), y = max(4), label = paste("Edf=", 7.38, ", F=", 17.81, ", P<", 0.001), hjust = 1, vjust = 0.5, size = 8) +
  theme_minimal()+
  theme(axis.line = element_line(color = "black"),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 29),
        panel.grid = element_blank())
plot1_c_data
plot1_d_data<- ggplot(d_data, aes(x = Time, y = -1 * PC1)) +
  geom_smooth(method = "gam", se = TRUE, color = "black", size = 1.5)  +
  geom_point(size = 1) +
  labs(title = ">250m", x = "Time(h)", y = "PC1") +
  annotate("text", x = max(23), y = max(4), label = paste("Edf=", 5.85, ", F=", 16.73, ", P<", 0.001), hjust = 1, vjust = 0.5, size = 8) +
  theme_minimal()+
  theme(axis.line = element_line(color = "black"),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 29),
        panel.grid = element_blank())
plot1_d_data

combined_plot1 <- ggplot() +
  geom_smooth(data = a_data, aes(x = Time, y = -1 * PC1, color = "In the forest"), method = "gam", se = TRUE, linetype = "twodash", size = 1.5) +
  geom_smooth(data = b_data, aes(x = Time, y = -1 * PC1, color = "0-100"), method = "gam", se = TRUE, linetype = "twodash", size = 1.5) +
  geom_smooth(data = c_data, aes(x = Time, y = -1 * PC1, color = "100-250"), method = "gam", se = TRUE, linetype = "twodash", size = 1.5) +
  geom_smooth(data = d_data, aes(x = Time, y = -1 * PC1, color = ">250"), method = "gam", se = TRUE, linetype = "twodash", size = 1.5) +
  labs(x = "Time(h)", y = "PC1") +
  scale_color_manual(
    values = c("In the forest" = "cyan3", "0-100" = "orange3", "100-250" = "green3", ">250" = "magenta4"), 
    labels = c("In the forest", "0-100", "100-250", ">250"),
    breaks = c("In the forest", "0-100", "100-250", ">250")
  ) +
  guides(color = guide_legend(title = "Forest distance(m)")) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    plot.margin = margin(0, 2, 0, 0, "cm"),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 25),
    axis.title = element_text(size = 25),
    plot.title = element_text(size = 29),
    legend.title = element_text(size = 25),  
    legend.text = element_text(size = 24),
    legend.position = "right",
    panel.grid = element_blank())

combined_plot1
library(gridExtra)
common_ylim <- c(-2, 4)  

plot1_a_data <- plot1_a_data + coord_cartesian(ylim = common_ylim)
plot1_b_data <- plot1_b_data + coord_cartesian(ylim = common_ylim)
plot1_c_data <- plot1_c_data + coord_cartesian(ylim = common_ylim)
plot1_d_data <- plot1_d_data + coord_cartesian(ylim = common_ylim)
row1 <- grid.arrange(plot1_a_data, plot1_b_data, nrow = 1)
row2 <- grid.arrange(plot1_c_data, plot1_d_data, nrow = 1)




plot2_a_data<- ggplot(a_data, aes(x = Time, y = -1 * PC2)) +
  geom_smooth(method = "gam", se = TRUE, color = "black", size = 1.5)  +
  geom_point(size = 1) +
  labs(title = "In the forest", x = "Time(h)", y = "PC2") +
  annotate("text", x = max(23), y = max(2), label = paste("Edf=", 4.25, ", F=", 22.66, ", P<", 0.001), hjust = 1, vjust = 0.5, size = 8) +
  theme_minimal()+
  theme(axis.line = element_line(color = "black"),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 29),
        panel.grid = element_blank())
plot2_a_data
plot2_b_data<- ggplot(b_data, aes(x = Time, y = -1 * PC2)) +
  geom_smooth(method = "gam", se = TRUE, color = "black", size = 1.5)  +
  geom_point(size = 1) +
  labs(title = "0-100m", x = "Time(h)", y = "PC2") +
  annotate("text", x = max(23), y = max(2), label = paste("Edf=", 6.02, ", F=", 15.44, ", P<", 0.001), hjust = 1, vjust = 0.5, size = 8) +
  theme_minimal()+
  theme(axis.line = element_line(color = "black"),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 29),
        panel.grid = element_blank())
plot2_b_data
plot2_c_data<- ggplot(c_data, aes(x = Time, y = -1 * PC2)) +
  geom_smooth(method = "gam", se = TRUE, color = "black", size = 1.5)  +
  geom_point(size = 1) +
  labs(title = "100-250m", x = "Time(h)", y = "PC2") +
  annotate("text", x = max(23), y = max(2), label = paste("Edf=", 5.61, ", F=", 12.12, ", P<", 0.001), hjust = 1, vjust = 0.5, size = 8) +
  theme_minimal()+
  theme(axis.line = element_line(color = "black"),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 29),
        panel.grid = element_blank())
plot2_c_data
plot2_d_data<- ggplot(d_data, aes(x = Time, y = -1 * PC2)) +
  geom_smooth(method = "gam", se = TRUE, color = "black", size = 1.5)  +
  geom_point(size = 1) +
  labs(title = ">250m", x = "Time(h)", y = "PC2") +
  annotate("text", x = max(23), y = max(2), label = paste("Edf=", 6.03, ", F=", 14.25, ", P<", 0.001), hjust = 1, vjust = 0.5, size = 8) +
  theme_minimal()+
  theme(axis.line = element_line(color = "black"),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 29),
        panel.grid = element_blank())
plot2_d_data
combined_plot2 <- ggplot() +
  geom_smooth(data = a_data, aes(x = Time, y = -1 * PC2, color = "In the forest"), method = "gam", se = TRUE, linetype = "twodash", size = 1.5) +
  geom_smooth(data = b_data, aes(x = Time, y = -1 * PC2, color = "0-100"), method = "gam", se = TRUE, linetype = "twodash", size = 1.5) +
  geom_smooth(data = c_data, aes(x = Time, y = -1 * PC2, color = "100-250"), method = "gam", se = TRUE, linetype = "twodash", size = 1.5) +
  geom_smooth(data = d_data, aes(x = Time, y = -1 * PC2, color = ">250"), method = "gam", se = TRUE, linetype = "twodash", size = 1.5) +
  labs(x = "Time(h)", y = "PC2") +
  scale_color_manual(
    values = c("In the forest" = "cyan3", "0-100" = "orange3", "100-250" = "green3", ">250" = "magenta4"), 
    labels = c("In the forest", "0-100", "100-250", ">250"),
    breaks = c("In the forest", "0-100", "100-250", ">250")
  ) +
  guides(color = guide_legend(title = "Forest distance(m)")) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    plot.margin = margin(0, 2, 0, 0, "cm"),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 25),
    axis.title = element_text(size = 25),
    plot.title = element_text(size = 29),
    legend.title = element_text(size = 25),  
    legend.text = element_text(size = 24),
    legend.position = "right",
    panel.grid = element_blank())

combined_plot2
common_ylim <- c(-1, 2) 
plot2_a_data <- plot2_a_data + coord_cartesian(ylim = common_ylim)
plot2_b_data <- plot2_b_data + coord_cartesian(ylim = common_ylim)
plot2_c_data <- plot2_c_data + coord_cartesian(ylim = common_ylim)
plot2_d_data <- plot2_d_data + coord_cartesian(ylim = common_ylim)
row3 <- grid.arrange(plot2_a_data, plot2_b_data, nrow = 1)

row4 <- grid.arrange(plot2_c_data, plot2_d_data, nrow = 1)

Dawn_data <- subset(data, Time == "6")
Dusk_data <- subset(data, Time == "18")
midan_data <- subset(data, Time == "12")
midnight_data <- subset(data, Time == "0")

#lm
pc1_Dawn <- lm(PC1 ~ Forest_DIS, data = Dawn_data)
pc1_Dusk <- lm(PC1 ~ Forest_DIS, data = Dusk_data)
pc1_midan <- lm(PC1 ~ Forest_DIS, data = midan_data)
pc1_midnight <- lm(PC1 ~ Forest_DIS, data = midnight_data)
summary(pc1_Dawn)
summary(pc1_Dusk)
summary(pc1_midan)
summary(pc1_midnight)

plot1_dawn<- ggplot(Dawn_data, aes(x = Forest_DIS, y = -1*PC1)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", level = 0.94)  +
  geom_point(size=1) +
  labs(x = " ", y = "PC1") +
  annotate("text", x = max(600), y = max(2), label = paste("P<", 0.001), hjust = 1, vjust = 0.5, size = 8) +
  theme_minimal()+
  theme(axis.line = element_line(color = "black"),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 29),
        panel.grid = element_blank())
plot1_dawn


plot1_dusk<- ggplot(Dusk_data, aes(x = Forest_DIS, y = -1 * PC1)) +
  geom_smooth(method = "lm", se = TRUE, color = "black",level=0.93)  +
  geom_point(size=1) +
  labs(x = " ", y = "PC1") +
  annotate("text", x = max(600), y = max(2), label = paste("P<", 0.01), hjust = 1, vjust = 0.5, size = 8) +
  theme_minimal()+
  theme(axis.line = element_line(color = "black"),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 29),
        panel.grid = element_blank())
plot1_dusk
plot1_midan<- ggplot(midan_data, aes(x = Forest_DIS, y = -1 * PC1)) +
  geom_smooth(method = "lm", se = TRUE, color = "black")  +
  geom_point(size=1) +
  labs(x = " ", y = "PC1") +
  annotate("text", x = max(600), y = max(2), label = paste("P=", 0.87), hjust = 1, vjust = 0.5, size = 8) +
  theme_minimal()+
  theme(axis.line = element_line(color = "black"),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 29),
        panel.grid = element_blank())
plot1_midan
plot1_midnight<- ggplot(midnight_data, aes(x = Forest_DIS, y = -1 * PC1)) +
  geom_smooth(method = "lm", se = TRUE, color = "black")  +
  geom_point(size=1) +
  labs(x = " ", y = "PC1") +
  annotate("text", x = max(600), y = max(2), label = paste("P<", 0.1), hjust = 1, vjust = 0.5, size = 8) +
  theme_minimal()+
  theme(axis.line = element_line(color = "black"),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 29),
        panel.grid = element_blank())
plot1_midnight

plot2_dawn<- ggplot(Dawn_data, aes(x = Forest_DIS, y = -1 * PC2)) +
  geom_smooth(method = "lm", se = TRUE, color = "black")  +
  geom_point(size=1) +
  labs(x = " ", y = "PC2") +
  annotate("text", x = max(600), y = max(1), label = paste("P=", 0.51), hjust = 1, vjust = 0.5, size = 8) +
  theme_minimal()+
  theme(axis.line = element_line(color = "black"),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 29),
        panel.grid = element_blank())+
  scale_y_continuous(labels = function(x) gsub("\\.0$", "", format(x, nsmall = 1)))
plot2_dawn

plot2_dusk<- ggplot(Dusk_data, aes(x = Forest_DIS, y = -1 * PC2)) +
  geom_smooth(method = "lm", se = TRUE, color = "black")  +
  geom_point(size=1) +
  labs(x = " ", y = "PC2") +
  annotate("text", x = max(600), y = max(1), label = paste("P<", 0.05), hjust = 1, vjust = 0.5, size = 8) +
  theme_minimal()+
  theme(axis.line = element_line(color = "black"),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 29),
        panel.grid = element_blank())+
  scale_y_continuous(labels = function(x) gsub("\\.0$", "", format(x, nsmall = 1)))
plot2_dusk

plot2_midan<- ggplot(midan_data, aes(x = Forest_DIS, y = -1 * PC2)) +
  geom_smooth(method = "lm", se = TRUE, color = "black",level=0.83)  +
  geom_point(size=1) +
  labs(x = " ", y = "PC2") +
  annotate("text", x = max(600), y = max(1), label = paste("P=", 0.79), hjust = 1, vjust = 0.5, size = 8) +
  theme_minimal()+
  theme(axis.line = element_line(color = "black"),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 29),
        panel.grid = element_blank())+
  scale_y_continuous(labels = function(x) gsub("\\.0$", "", format(x, nsmall = 1)))
plot2_midan

plot2_midnight<- ggplot(midnight_data, aes(x = Forest_DIS, y = -1 * PC2)) +
  geom_smooth(method = "lm", se = TRUE, color = "black")  +
  geom_point(size=1) +
  labs(x = " ", y = "PC2") +
  annotate("text", x = max(600), y = max(1), label = paste("P<", 0.001), hjust = 1, vjust = 0.5, size = 8) +
  theme_minimal()+
  theme(axis.line = element_line(color = "black"),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 29),
        panel.grid = element_blank())+
  scale_y_continuous(labels = function(x) gsub("\\.0$", "", format(x, nsmall = 1)))
plot2_midnight
common_ylim <- c(-2, 2)

plot1_dawn_data <- plot1_dawn + coord_cartesian(ylim = common_ylim)
plot1_dusk_data <- plot1_dusk + coord_cartesian(ylim = common_ylim)
plot1_midan_data <- plot1_midan + coord_cartesian(ylim = common_ylim)
plot1_midnight_data <- plot1_midnight + coord_cartesian(ylim = common_ylim)

common_ylim <- c(-1, 1) 

plot2_dawn_data <- plot2_dawn + coord_cartesian(ylim = common_ylim)
plot2_dusk_data <- plot2_dusk + coord_cartesian(ylim = common_ylim)
plot2_midan_data <- plot2_midan + coord_cartesian(ylim = common_ylim)
plot2_midnight_data <- plot2_midnight + coord_cartesian(ylim = common_ylim)
row1 <- grid.arrange(plot1_dawn_data, plot2_dawn_data, nrow = 1)
row2 <- grid.arrange(plot1_dusk_data, plot2_dusk_data, nrow = 1)
row3 <- grid.arrange(plot1_midan_data, plot2_midan_data, nrow = 1)
row4 <- grid.arrange(plot1_midnight_data, plot2_midnight_data, nrow = 1)


#PCoA analysis 
setwd("C:/Users/xinto/Desktop/Data output") 
data <- read.csv("Point_counts.csv") 
richness <- table(data$Site)
library(sp)      # 用于空间操作

# 计算物种丰富度
richness <- tapply(data$Species, data$Site, FUN = function(x) length(unique(x)))

# 提取距离森林数据
forest_distance <- data$Distance

# 创建数据框
df <- data.frame(Richness = richness, DistanceToForest = forest_distance)
Rich <- data1$Richness
DIS <- data$NEAR_DIS 
# 执行相关性分析
correlation <- cor(Rich, DIS, use = "complete.obs")

# 绘制回归线图
ggplot(df, aes(x = DistanceToForest, y = Richness)) +
  geom_point() +   
  geom_smooth(method = "lm", se = FALSE) +  
  labs(x = "Distance to Forest", y = "Species Richness") +   
  ggtitle("Correlation between Species Richness and Distance to Forest") 

# 计算物种丰富度
richness <- tapply(data$Species, data$Site, FUN = function(x) length(unique(x)))

# Compute mean species occurrence for each site
site_occurrence <- data %>% 
  group_by(Site) %>% 
  summarize(Species = n_distinct(Species))

# Convert data to presence-absence matrix
presence_absence <- xtabs(~Site + Species, data = data)

# Calculate Sorensen dissimilarity index
sorensen_index <- vegdist(presence_absence, method = "bray")

# Perform PCoA
pcoa <- cmdscale(sorensen_index, k = 2)
str(data)
# Combine PCoA results with site occurrence data
pcoa_plot <- data.frame(Site = as.character(site_occurrence$Site), PC1 = pcoa[, 1], PC2 = pcoa[, 2])

pcoa_plot$Category <- data$Cata

# Plot PCoA results
labels <- c("a" = "In the forest", "b" = "0-100", "c" = "100-250", "d" = ">250", size = 8)
colors <- c("a" = "red", "b" = "blue", "c" = "green", "d" = "orange")


PCoA <- ggplot(data, aes(x = PC1, y = PC2, color = Cata)) +
  geom_point() +
  labs(title = "PCoA of Community Composition(Sorensen dissimilarity index)", color = "Distance from points to forest(m)", size = 20) +
  scale_color_manual(values = colors, labels = labels)+
  theme(axis.title.x = element_text(size = 20),  
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 19)) 


PCoA <- ggplot(data, aes(x = PC1, y = PC2, color = Site)) +
  geom_point() +
  labs(title = "PCoA of Community Composition(Sorensen dissimilarity index)", color = "Site") 

PCoA

site_diversity <- diversity(table(data$Site, data$Species), index = "shannon")


Shannon <- ggplot(data, aes(x = Distance, y = Shannon)) +
  geom_point() +
  labs(x = " ",
       y = "Shannon Index")+
  annotate("text", x = 600, y = 3.5, label = expression(paste(R^2, " = 0.24, ", P < 0.001)), hjust = 1, vjust = 1, size = 6) +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(size = 18),
        panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.line = element_line(color = "black", size = 0.5))+
  geom_smooth(method = "lm", se = TRUE, color = "black", level = 0.95)

Shannon


model <- lm(Richness ~ Distance, data = data)
model <- lm(Shannon ~ Distance, data = data)
model <- lm(PC1 ~ Distance, data = data)
model <- lm(PC2 ~ Distance, data = data)


Richness <- ggplot(data, aes(x = Distance, y = Richness)) +
  geom_point() +
  labs(x = " ",
       y = "Richness")+
  annotate("text", x = 600, y = 40, label = expression(paste(R^2, " = 0.17, ", P < 0.001)), hjust = 1, vjust = 1, size = 6) +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(size = 18),
        panel.background = element_blank(),
        panel.grid = element_blank(),  # 移除背景的网格线
        axis.line = element_line(color = "black", size = 0.5))+
  geom_smooth(method = "lm", se = TRUE, color = "black", level = 0.95) # 添加轴线


Richness


PC2 <- ggplot(data, aes(x = Distance, y = PC2)) +
  geom_point() +
  labs(x = " ", y = "PC2") +
  annotate("text", x = 600, y = 0.48,
           label = expression(paste(R^2, " = 0.23, ", P < 0.001)),
           hjust = 1, vjust = 1, size = 6) +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 18),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  ) +
  geom_smooth(method = "lm", se = TRUE, color = "black", level = 0.5)

PC2

PC1 <- ggplot(data, aes(x = Distance, y = PC1)) +
  geom_point() +
  labs(x = " ", y = "PC1") +
  annotate("text", x = 600, y = 0.50, label = paste("P <", 0.001), hjust = 1, vjust = 1, size = 6) +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(size = 18),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black", size = 0.5)) +
  geom_smooth(method = "gam", se = TRUE, color = "black", level = 0.95) 

PC1 
row <- grid.arrange(Shannon, Richness, nrow = 1)

model <- segmented(lm(PC1 ~ Distance, data = data), seg.Z = ~ Distance, psi = 240)
summary(model)
plot(model, xlab = "Distance to forest(m)", ylab = "PC1", main = "Piecewise Regression", color = "black")
plot(data$Distance, data$PC1,  xlab = "Distance", ylab = "PC1", main = "Piecewise Regression")
lines(model, col = "red")
library(segmented)

model <- segmented(lm(PC1 ~ Distance, data = data), seg.Z = ~ Distance, psi = 200)
summary(model)

points(data$Distance, data$PC1)
points(data$Distance, predict(model), col = "red", type = "l")

data_segment1 <- subset(data, Distance <= 239)
data_segment2 <- subset(data, Distance > 239)

PC1 <- ggplot(data, aes(x = Distance, y = PC1)) +
  geom_point() +
  labs(x = " ", y = "PC1") +
  annotate("text", x = 600, y = 0.50, label = expression(paste(R^2, " = 0.71, ", P < 0.001)), hjust = 1, vjust = 1, size = 6) +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 18),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

PC1 <- PC1 +
  geom_smooth(data = data_segment1, method = "lm", formula = y ~ x, se = TRUE, color = "black", level = 0.9999) +
  geom_smooth(data = data_segment2, method = "lm", formula = y ~ x, se = TRUE, color = "black", level = 0.9999)

PC1