setwd("/Users/manishkumar/Downloads/")
library(BoolNet)
bRules <- loadNetwork("net.txt")
bRules
plotSequence(bRules, startState = c("SNCA"=1, "SNCAagg" = 0,"ATP"=1,
                                    "MP"=1, "Ca"=1, "PTP" = 0, "ROS" = 0),
             includeAttractorStates = "all")




Normal states: (&quot;PINK1&quot; = 1, &quot;SNCAagg&quot; = 0, &quot;NCLX&quot; = 1, &quot;Ca&quot; = 0, &quot;NOX&quot; = 0, &quot;ROS&quot; = 0, &quot;PTP&quot; = 0,
                &quot;MP&quot; = 1, &quot;ATP&quot; = 1)


bRules
state <- generateState(bRules, 
                       c("SNCA"=1, "SNCAagg" = 0,"ATP"=1,
                         "MP"=1, "Ca"=0, "PTP" = 0, "ROS" = 0)
)
att1=getAttractors(bRules, startStates=list(state))
att1
att1[1]
plotAttractors(att1)
path = getPathToAttractor(bRules, state, includeAttractorStates = "all")
plotAttractors(path)

plotNetworkWiring(bRules)

plotSequence(bRules, startState = c("SNCA"=1, "SNCAagg" = 0,"ATP"=1,
                                    "MP"=1, "Ca"=1, "PTP" = 0, "ROS" = 0),
             includeAttractorStates = "all")

plotSequence(sequence = path, mode = "table")


generateTimeSeries(bRules, numSeries = 1, numMeasurements = 10)
