
remove(list=ls())
library(writexl)
library(ggplot2)
library(tidyverse)

# Reading in data
y_data <- read.csv("CCLE_Ydata.csv")
y_data <- y_data[, colSums(apply(y_data, 2, is.na))==0]
nfolds <- 10
method_num <- 10

# Reading in the results
results <- readRDS("Drug_results_new.rds")

mspe <- array(unlist(results$mspe), dim = c(nfolds, ncol(y_data), method_num))
mad <- array(unlist(results$mad), dim = c(nfolds, ncol(y_data), method_num))
pred_data <- array(unlist(results$pred_data), dim = c(nrow(y_data), ncol(y_data), method_num)) 
LR_test <- unlist(results$lr_test)
ecp <-  array(unlist(results$ecp), dim = c(nrow(y_data), ncol(y_data), 4))
PI_len  <- array(unlist(results$PI_len), dim = c(nrow(y_data), ncol(y_data), 4))

##### Processing Results #####
sd_mspe_outcome <- do.call(cbind, (lapply(apply(mspe, 3, list), function(x){apply(x[[1]], 2, sd)})))/sqrt(nfolds)
sd_mad_outcome <- do.call(cbind, (lapply(apply(mad, 3, list), function(x){apply(x[[1]], 2, sd)})))/sqrt(nfolds)

test_array <- array( as.matrix(results$test_data), dim = dim(pred_data) )
squared_errors <- (test_array -pred_data)^2
mspe_outcome <- apply(squared_errors, 3, colMeans, na.rm=T)

abs_errors <- abs(test_array-pred_data)
mads_outcome <- apply(abs_errors, 3, function(x){apply(x, 2, median, na.rm=T)})

methods <- c("PROBE", "MCP", "SCAD", "LASSO", "ALASSO", "EBREG", "VARBVS", "SSLASSO", "SPARSEVB", "PROBE.1")
drugs <- names(y_data)
colnames(sd_mspe_outcome) = colnames(sd_mad_outcome) = colnames(mspe_outcome) = colnames(mads_outcome)  <- methods
rownames(sd_mspe_outcome) = rownames(sd_mad_outcome) = rownames(mspe_outcome) = rownames(mads_outcome) <-  drugs

#Proportion of test PIs that contain the test observation for PROBE
ecp_outcome <- apply(ecp, 3, colMeans, na.rm=T)
PILen_outcome <- apply(PI_len, 3, colMeans, na.rm=T)

#Mean PI_length by fold
fold_num <- results$fold_num
mean_PI <- array(0, dim = c(ncol(y_data),nfolds,4))
for(i in 1:8){
  for(k in 1:nfolds){
    mean_PI[i,k,] <- apply( PI_len[fold_num == k,i,], 2, mean)
  }
}

#SD of PI_length by fold
sd_PI_length <- do.call(cbind, (lapply(apply(mean_PI, 3, list), function(x){apply(x[[1]], 1, sd)})))

Method_PI <- c("PROBE", "Conform Split", "Conform Jackknife", "EBREG")
colnames(PILen_outcome) = colnames(ecp_outcome) = colnames(sd_PI_length)  <- Method_PI

#Save excel files with output
sheets <- list("MPSE" = data.frame(drugs,mspe_outcome), 
               "MPSE SD" = data.frame(drugs,sd_mspe_outcome), 
               "MAD" = data.frame(drugs,mads_outcome), 
               "MAD SD" = data.frame(drugs,sd_mad_outcome), 
               "ECP PI" =data.frame(drugs, ecp_outcome), 
               "PI Length" =data.frame(drugs, PILen_outcome), 
               "SD PI Length"= data.frame(drugs, sd_PI_length))

##### Exporting results to excel file  #####
write_xlsx(sheets, "Table Results.xlsx")


## Best MPSE and MAD
best_MPSE <- methods[ apply( mspe_outcome, 1, which.min)]
best_MAD <- methods[ apply( mads_outcome, 1, which.min)]
data.frame(drugs, best_MPSE, best_MAD)





##### Create Plots of MPSE and MAD.  #####


cols4methods <- c("black","firebrick1" ,  "darkgoldenrod4" , "maroon2",         "cornflowerblue", "chartreuse4" ,            "palevioletred2",  "lightsteelblue4"   ,   "blue1", "coral")
# cols4methods <- c("black","firebrick1" ,  "darkgoldenrod4" , "maroon2",         "cornflowerblue" ,            "palevioletred2",  "lightsteelblue4"   ,   "blue1", "coral")
cols4methods_nosvb <- cols4methods[-c(7)]

shape_vals <- c(15:24)
shape_vals_nosvb <- shape_vals[-c(7)]

pos_wid <- 0.88
pos_siz <- 3

#Plot for MPSE
mspe_data <- data.frame(cbind(unlist(data.frame(mspe_outcome)), unlist(data.frame(sd_mspe_outcome)), rep(methods, each=8), rep(drugs, length(methods))))
names(mspe_data) <- c("MPSE", "SD", "Method", "Drug")
mspe_data$Drug <- gsub("_ActArea", "", mspe_data$Drug)
mspe_data$Drug <- gsub("\\.", "-", mspe_data$Drug)
mspe_data[,1] <- as.numeric(mspe_data[,1])
mspe_data[,2] <- as.numeric(mspe_data[,2])
mspe_data$lower <- mspe_data$MPSE-mspe_data$SD
mspe_data$upper <- mspe_data$MPSE+mspe_data$SD


ggplot(mspe_data, aes(x=Drug, y=MPSE, color = Method)) + 
  scale_color_manual(values=cols4methods) +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = pos_wid, 
                position=position_dodge(width=pos_wid)) +
  geom_point(aes(shape=Method), size=pos_siz, position=position_dodge(width=pos_wid)) + 
  scale_shape_manual(values=shape_vals) +
theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(angle = 20,hjust=1,vjust=0.98,size=12), 
        axis.text.y=element_text(size=12), 
        axis.title = element_text(color="#666666",
                                  face="bold", size=18),
        strip.text.x = element_text(size=12),
        strip.text.y = element_text(size=12))
ggsave("Example1_MSPE.pdf", dpi = 600, width = 15*0.9, height = 6*0.9)




#Plot for MAD
mad_data <- data.frame(cbind(unlist(data.frame(mads_outcome)), unlist(data.frame(sd_mad_outcome)), rep(methods, each=8), rep(drugs, length(methods))))
names(mad_data) <- c("MAD", "SD", "Method", "Drug")
mad_data$Drug <- gsub("_ActArea", "", mad_data$Drug)
mad_data$Drug <- gsub("\\.", "-", mad_data$Drug)
mad_data[,1] <- as.numeric(mad_data[,1])
mad_data[,2] <- as.numeric(mad_data[,2])
mad_data$lower <- mad_data$MAD-mad_data$SD
mad_data$upper <- mad_data$MAD+mad_data$SD


ggplot(mad_data, aes(x=Drug, y=MAD, color = Method)) + 
  scale_color_manual(values=cols4methods) +
  theme_bw() + 
  theme_light() +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = pos_wid, 
                position=position_dodge(width=pos_wid)) +
  geom_point(aes(shape=Method), size=pos_siz, position=position_dodge(width=pos_wid)) + 
  scale_shape_manual(values=shape_vals) +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(angle = 20,hjust=1,vjust=0.98,size=12), 
        axis.text.y=element_text(size=12), 
        axis.title = element_text(color="#666666",
                                  face="bold", size=18),
        strip.text.x = element_text(size=12),
        strip.text.y = element_text(size=12))
ggsave("Example1_MAD.pdf", dpi = 600, width = 15*0.9, height = 6*0.9)



### Redo without SparseVB
mspe_data <- mspe_data[mspe_data$Method != 'SPARSEVB',]


ggplot(mspe_data, aes(x=Drug, y=MPSE, color = Method)) + 
  scale_color_manual(values=cols4methods_nosvb) +
  theme_bw() + 
  theme_light() +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = pos_wid, 
                position=position_dodge(width=pos_wid)) +
  geom_point(aes(shape=Method), size=pos_siz, position=position_dodge(width=pos_wid)) + 
  scale_shape_manual(values=shape_vals_nosvb) +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(angle = 20,hjust=1,vjust=0.98,size=12), 
        axis.text.y=element_text(size=12), 
        axis.title = element_text(color="#666666",
                                  face="bold", size=18),
        strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14))
ggsave("Example1_MSPE_nosparsevb.pdf", dpi = 600, width = 15*0.9, height = 6*0.9)


#Plot for MAD
mad_data <- mad_data[mad_data$Method != 'SPARSEVB',]
ggplot(mad_data, aes(x=Drug, y=MAD, color = Method)) + 
  scale_color_manual(values=cols4methods_nosvb) +
  theme_bw() + 
  theme_light() +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = pos_wid, 
                position=position_dodge(width=pos_wid)) +
  geom_point(aes(shape=Method), size=pos_siz, position=position_dodge(width=pos_wid)) + 
  scale_shape_manual(values=shape_vals_nosvb) +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(angle = 20,hjust=1,vjust=0.98,size=12), 
        axis.text.y=element_text(size=12), 
        axis.title = element_text(color="#666666",
                                  face="bold", size=18),
        strip.text.x = element_text(size=12),
        strip.text.y = element_text(size=12))
ggsave("Example1_MAD_nosparsevb.pdf", dpi = 600, width = 15*0.9, height = 6*0.9)




#Plot for MSPE and MAD

mspe_data <- mspe_data[mspe_data$Method != 'SPARSEVB',]
mad_data <- mad_data[mad_data$Method != 'SPARSEVB',]

mspe_data <- mspe_data %>% 
  rename(Error = MPSE) %>% 
  mutate(Measure = "MSPE")

mad_data <- mad_data %>% 
  rename(Error = MAD) %>% 
  mutate(Measure = "MAD")

err_data <- mspe_data %>% 
  bind_rows(mad_data)  %>%
  mutate(across(Measure, factor, levels=c("MSPE","MAD"))) 

ggplot(err_data, aes(x=Drug, y=Error, color = Method)) + 
  scale_color_manual(values=cols4methods_nosvb) +
  theme_bw() + 
  theme_light() +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = pos_wid, 
                position=position_dodge(width=pos_wid)) +
  geom_point(aes(shape=Method), size=pos_siz, position=position_dodge(width=pos_wid)) + 
  scale_shape_manual(values=shape_vals_nosvb) +
  facet_wrap(~Measure,ncol = 1, scales = "free") +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(angle = 20,hjust=1,vjust=0.98,size=12), 
        axis.text.y =element_text(size=12), 
        axis.title.y = element_text(size=0),
        axis.title.x = element_text(color="#666666",
                                  face="bold", size=18),
        strip.text.x = element_text(color="#666666",
                                    face="bold", size=18))
ggsave("Example1_MADMSPE_nosparsevb.pdf", dpi = 600, width = 15*0.9, height = 12*0.9)


