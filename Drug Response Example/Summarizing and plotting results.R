
remove(list=ls())
library(writexl)
library(ggplot2)

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
sd_mspe_outcome <- do.call(cbind, (lapply(apply(mspe, 3, list), function(x){apply(x[[1]], 2, sd)})))
sd_mad_outcome <- do.call(cbind, (lapply(apply(mad, 3, list), function(x){apply(x[[1]], 2, sd)})))

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
  theme_bw() + 
  theme_light() +
  theme(axis.text.x = element_text(angle = 30,hjust=1,vjust=0.98,size=15), axis.text.y=element_text(size=15), axis.title=element_text(size=17)) +
  theme(legend.text=element_text(size=15), legend.title=element_text(size=17) )
ggsave("Example1_MSPE.pdf", dpi = 600, width = 15*0.9, height = 6*0.9)
dev.off()



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
  theme(axis.text.x = element_text(angle = 30,hjust=1,vjust=0.98,size=15), axis.text.y=element_text(size=15), axis.title=element_text(size=17)) +
  theme(legend.text=element_text(size=15), legend.title=element_text(size=17) )
ggsave("Example1_MAD.pdf", dpi = 600, width = 15*0.9, height = 6*0.9)
dev.off()


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
  
  theme(axis.text.x = element_text(angle = 30,hjust=1,vjust=0.98,size=15), axis.text.y=element_text(size=15), axis.title=element_text(size=17)) +
  theme(legend.text=element_text(size=15), legend.title=element_text(size=17) )
ggsave("Example1_MSPE_nosparsevb.pdf", dpi = 600, width = 15*0.9, height = 6*0.9)
dev.off()

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
  theme(axis.text.x = element_text(angle = 30,hjust=1,vjust=0.98,size=15), axis.text.y=element_text(size=15), axis.title=element_text(size=17)) +
  theme(legend.text=element_text(size=15), legend.title=element_text(size=17) )
ggsave("Example1_MAD_nosparsevb.pdf", dpi = 600, width = 15*0.9, height = 6*0.9)
dev.off()







##### Create Plots of PI length by fold #####

drugs <- gsub("_ActArea", "", drugs)
drugs <- gsub("\\.", "-", drugs)

## Adding Jackknife data to pred_data
conf_jack_w_cv <- readRDS("Conform_lasso_jack.rds")
EBreg_unit <- readRDS("EBreg_unif.rds")
pred_data2 <- array(unlist(results$pred_data), dim = c(nrow(y_data), ncol(y_data), 6)) 
pred_data <- array(0, dim = c(nrow(y_data), ncol(y_data), 8)) 
pred_data[,,1:6] <- pred_data2
pred_data[,,7] <- as.matrix(conf_jack_w_cv$Pred)
pred_data[,,8] <- as.matrix(EBreg_unit$Pred)

## Creating a data.frame with the results
drug_ind <- which.max(drugs == "PD-0325901")
plot_mat <- data.frame(ID = rep(1:nrow(y_data),4), Fold = as.factor(rep(fold_num,4)), 
                       Method = factor( rep( c("PROBE", "Conformal Split", 
                       "Conformal Jackknife", "EBREG"), each = nrow(y_data)),
                                        levels = c("PROBE", "Conformal Split", 
                                                   "Conformal Jackknife", "EBREG")),
                       Outcome = c(drugs[drug_ind]),  Y = y_data[,drug_ind], Pred = 
                         array( pred_data[, drug_ind, c(1,6,7, 8)]), PI_Length = 
                         array(PI_len[,drug_ind,]) )



# Plot of PI length by fold
limits <- aes(y = PI_Length, x = ID, ymax = PI_Length, ymin=0*PI_Length, color = Fold) 

p <- ggplot(data=plot_mat, aes(y = PI_Length, x = ID, color = Fold)) + 
  geom_point(size = 0.7) + facet_wrap(~Method, nrow =4,scales="fixed") 

p <- p +  labs(x="", y=paste("Prediction Interval Length"),size=60) + theme_minimal() + 
  theme(axis.text=element_text(size=15), axis.title = element_text(size=17), axis.text.x=element_text(size=0)) + 
  scale_colour_hue(l = 50, c = 100)+ 
  theme(strip.text.x = element_text(size = 17, colour = "black"))+
  theme(legend.text=element_text(size=15), legend.title=element_text(size=17) )

print(p)
ggsave("PI_length_fold.pdf", dpi = 600, width = 14, height = 5)
dev.off()


# Plot of PI length by predicted value

p <- ggplot(subset(plot_mat), aes(x = Pred, y = PI_Length, color = Fold)) + 
  geom_point() + facet_wrap(~Method, nrow =1,scales="fixed") 


p <- p +  labs(x="Predicted Value", y=paste("Prediction Interval Length"),size=60) + 
  theme_minimal() + theme(axis.text=element_text(size=14), axis.title = element_text(size=17)) +
  scale_colour_hue(l = 50, c = 100)+ 
  theme(strip.text.x = element_text(size = 14, colour = "black"))+
  theme(legend.text=element_text(size=15), legend.title=element_text(size=17) )
print(p)

ggsave("PI_length_by_pred.pdf", dpi = 600, width = 14, height = 5)
dev.off()

# Plot of PI's by predicted value by fold for the first 4 folds
limits <- aes(y = Pred, x = Pred, ymax = Pred + PI_Length/2, ymin = 
                Pred - PI_Length/2, color = Fold) 

p <- ggplot(subset(plot_mat, Fold %in% c(1:4)), aes(y = Pred, x = ID)) + 
  geom_pointrange(limits,fatten=0.1,size = 0.7) + facet_wrap(~Method, nrow =1,scales="fixed") 


cols4methods <- c("firebrick2", "dodgerblue1", "snow3","green4", "slateblue1")

p <- p +  labs(x="Predicted Value", y=paste("Prediction Interval"),size=60) + 
  theme_minimal() + theme(axis.text=element_text(size=14), axis.title = element_text(size=17)) +
  scale_color_manual(values=cols4methods) +
  theme(strip.text.x = element_text(size = 14, colour = "black"))+
  theme(legend.text=element_text(size=15), legend.title=element_text(size=17) )
print(p)

ggsave("PI_length_pred.pdf", dpi = 600, width = 14, height = 5)
dev.off()



