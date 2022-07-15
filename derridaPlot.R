setwd("C:\\Users\\Pankaz-Lab-PC1\\Desktop\\baby\\boolean analysis")

library(BoolNet)
brules = loadNetwork("Total Factors Combined.txt")

normal= plotSequence(brules, startState = c("DJ1" = 1, "PINK1" = 1, "SNCA" = 1, "Ca" = 0,  "ROS" = 0, "PTP" = 0, "SNCAagg" = 0, "MP" = 1, "ATP" = 1, "DJBP" = 1, "NCX" = 1, "NOX" = 0, "SOD" = 1, "LAMP2A" = 1, "FFASB" = 1)
                       , includeAttractorStates = "all")

disease= plotSequence(brules, startState = c("DJ1" = 1, "PINK1" = 0, "SNCA" = 1, "Ca" = 0,  "ROS" = 0, "PTP" = 0, "SNCAagg" = 0, "MP" = 1, "ATP" = 1, "DJBP" = 1, "NCX" = 1, "NOX" = 0, "SOD" = 1, "LAMP2A" = 1, "FFASB" = 1)
                  , includeAttractorStates = "all")


flip = function(x){
  return(ifelse(x==0,1,0))
}

i = c(1,2)
ini1 = sample(c(0,1),29, rep = T)
ini2 = ini1
ini2[i] = flip(ini2[i])

hammingDist = data.frame(c())
for (h in 1:27) {
  hammingTable = data.frame(c())
  for (j in 1:500) {
    x = c(1:29)
    x = x[-i]
    flipPos = sample(x,h)
    inp1 = ini1
    inp2 = ini2
    inp1[flipPos] = flip(inp1[flipPos])
    inp2[flipPos] = flip(inp2[flipPos])
    
    results1= plotSequence(brules, startState = inp1, includeAttractorStates = "all")
    op1 = rev(as.vector(results1[,2]))
    
    results2= plotSequence(brules, startState = inp2, includeAttractorStates = "all")
    op2 = rev(as.vector(results2[,2]))
    
    hd = dist(rbind(op1,op2), method = "minkowski", p=1)/length(ini1)
    
    r = c(op1,op2,hd)
    hammingTable = rbind(hammingTable,r) #STORE
    
  }
  
  colnames(hammingTable)  = c(paste0(brules$genes[i],"_OP1_",c(1:29)),
      paste0(brules$genes[i],"_OP2_",c(1:29)),
      "Dist")
  hammingTable
  meanHamming = mean(hammingTable$Dist) 
  hammingDist = rbind(hammingDist, c(h,meanHamming)) #store
}
hammingDist
write.csv(hammingDist, "hammingDist_g1g2.csv", row.names = F)
write.csv(hammingTable, "hammingTable_g1g2.csv", row.names = F)
colnames(hammingDist) = c("ht", "ht+1")
plot(hammingDist$ht, hammingDist$`ht+1`, type = "b")
