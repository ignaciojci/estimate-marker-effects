library(dplyr)
library(readxl)
library(gaston)
library(BGLR)

setwd("~/projects/PAS2444/jignacio/2025/gs-predictions")

# Get BLUES
cohorts <- paste0("OH",13:19)

blues_file <- "data/ALL BLUES BY COHORT.xlsx" # file from Clay

df_list <- lapply(cohorts, function(cohort){
  read_excel(blues_file, sheet = cohort, na = c("",".","NA")) %>%
    mutate(COHORT = cohort, .before=NAME)
})
names(df_list) <- cohorts
str(df_list)

blues_df <- do.call(rbind, df_list)

# Get marker data
vcf_data <- read.vcf("data/06_oh13-22_lines_imp.vcf.gz", convert.chr = F)
sample_names <- vcf_data@ped$famid

# Match sample names
match_idx <- match(blues_df$NAME, sample_names)
has_marker_data <- which(!is.na(match_idx))
has_marker_data_match_idx <- match_idx[has_marker_data]

blues_df_with_marker <- blues_df[has_marker_data,]
table(blues_df_with_marker$COHORT)
marker_data <- as.matrix(vcf_data)[has_marker_data_match_idx,] 
length(which(duplicated(rownames(marker_data)))) # no duplicate names

# LD Prune
plink <- "~/softwares/plink_1.9/plink" # specify path to plink

X.tmp <- as.matrix(marker_data)
ifelse(!dir.exists("tmp"),dir.create("tmp"),print("dir exists"))
ldfile <- "tmp/ld.filt"
write.bed.matrix(vcf_data,ldfile)

system(paste(plink,"--bfile",ldfile,"--allow-extra-chr","--indep-pairwise 250 10 0.4","--out",ldfile))
ldfilt <- read.table(paste0(ldfile,".prune.in"))
X <- X.tmp[,ldfilt$V1]

# Read list of markers of interest
target_markers_file <- "data/target_markers.csv" ### CHANGE THIS FILE TO CONTAIN THE ACTUAL TARGET MARKERS
#write.table(data.frame(colnames(X)[sample(1:ncol(X), 89)]), "data/target_markers.csv", sep=",", row.names = F, quote = F, col.names = F)
target_markers <- unlist(read.csv(target_markers_file, header=F)[,1])
non_target_markers <- colnames(X)[!colnames(X) %in% target_markers]

# Run parameters
nIter = 1500 # change to 15,000 when running final analysis
burnIn = 500 # change to 3,000 when running final analysis

out_matrix <- matrix(NA, nrow=3, ncol=5)
colnames(out_matrix) <- c("test", "target mean effect", "non-target mean effect", "t-test p-value", "accuracy")

# 1.	All 1,173 markers as random
test_num = 1
y = blues_df_with_marker$YLD
ETA = list(list(X=X, model="BRR"))
m <- BGLR(y = y, ETA = ETA)

b <- m$ETA[[1]]$b # marker effects
b_target <- b[target_markers]
b_non_target <- b[non_target_markers]

get_bglr_stats <- function(b_target, b_non_target = NA, i, mat = out_matrix){
  data <- data.frame(
    values = c(b_target, b_non_target),
    group = c(rep("target",length(b_target)), rep("non-target",length(b_non_target)))
  )
  
  boxplot(values ~ group, data = data,
          main = paste("Marker effects per marker group of scenario ", i),
          ylab = "marker effect",
          col = c("lightblue", "lightgreen"))
  
  mat[i,1] <- i
  mat[i,2] <- mean(b_target)
  mat[i,3] <- mean(b_non_target)
  tryCatch(
    mat[i,4] <- t.test(b_target, b_non_target)$p.value,
    error = function(e) e
  )
  mat[i,5] <- cor(y, m$yHat, use="complete.obs")
  return(mat)
}

out_matrix <- get_bglr_stats(b_target, b_non_target, test_num)

# 2.	89 consistent markers as fixed, perhaps weighted, all other random (need to use BGLR?)
test_num = 2
y = blues_df_with_marker$YLD
X1 = X[,target_markers]
X2 = X[,non_target_markers]

ETA = list(list(X=X1, model="FIXED"),
           list(X=X2, model="BRR"))
m <- BGLR(y = y, ETA = ETA)

b_target <- m$ETA[[1]]$b
b_non_target <- m$ETA[[2]]$b

out_matrix <- get_bglr_stats(b_target, b_non_target, test_num)

# 3.	Just the 89
test_num = 3
y = blues_df_with_marker$YLD
X1 = X[,target_markers]

ETA = list(list(X=X1, model="FIXED"))
m <- BGLR(y = y, ETA = ETA)

b_target <- m$ETA[[1]]$b
b_non_target <- NA

out_matrix <- get_bglr_stats(b_target, b_non_target, test_num)

out_matrix