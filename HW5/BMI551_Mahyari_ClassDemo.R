##############################
#This part will be demonstrated in class
#These have been coded and adapted from various sources
#BMI551W16 Eisa Mahyari
#########

#Markov Chain Class Demo

#A machine works in 3 settings (states): low, medium, high.
#Every hour, it either stays in the same state or switches to another; this transition is according to a 
#probability transition matrix:

TM = matrix(c(1/3,1/3,1/3,
              0.1,.6,.3,
              .6,.1,.3), ncol=3, byrow=TRUE);



#What is the 2-step transition matrix?

TM_2 <- TM %*% TM
TM_2

#If the system is in low state at 1:00AM, what is the probability it will be in the same mode at 5:00AM?

#thus we want the 4-stem
TM_4 <- TM %*% TM %*% TM %*% TM
TM_4
paste("The likelihood is %",round(TM_4[1,1]*100,2), sep="")


#############
#Markov simulation Example:

#Start State probability vector

ST <- c(1/3,1/3,1/3)

# Markov chain simultation function
MCsim = function(trans, start, N){
  
  Allstates <- numeric()
  currentState <- sample(1:3, 1, prob=start)
  cat("**** Starting on state:", currentState, "\n**********\n");
  # Make twenty steps through the markov chain
  Allstates[1] <- currentState
  transN <- trans
  for (i in 1:N)
  {
    
    newState <- sample(1:ncol(trans), 1, prob=trans[currentState,])
    cat(i, " ***", currentState, "->", newState, "\n");
    currentState = newState;
    Allstates[i+1] <- currentState
    transN <- transN %*% trans
  }
  cat("**** Ending on state:", currentState, "\n******************\n");
  return(list(states=Allstates, finalTM=transN))
}


LSstates <- MCsim(trans=TM, start=ST, N=30)
TM




##########End of MarkovChain Example


#How to get DNA Sequences from GenBank 
#And use the HMM package 

##########First, the HMM with manual data

#install.packages("HMM")

library(HMM)


hmm = initHMM(c("AT-rich", "GC-rich", "Neither"),c("A","T","G","C"),
              transProbs=matrix(c(0.4, 0.1, 0.5,
                                  0.1, 0.4, 0.5,
                                  0.3, 0.3, 0.4),
                                byrow=TRUE, nrow=3),
              emissionProbs=matrix(c(0.4,0.4,0.1,0.1,
                                     0.1,0.1,0.4,0.4,
                                     0.25,0.25,0.25,0.25),nrow=3, byrow=T),
              startProbs=c(1/3,1/3,1/3))
print(hmm)


SimSeqA  <- simHMM(hmm, length=100)
SimSeqB  <- simHMM(hmm, length=100)


fithmm = initHMM(c("AT-rich", "GC-rich", "Neither"),
                 c("A","T","G","C"))


print(fithmm)

observation = c(SimSeqA$observation,SimSeqB$observation)

# Baum-Welch
bw = baumWelch(fithmm,observation)


print(bw$hmm)
print(bw$difference)

observation

#log backward probabilities calc
logBackwardProbabilities = backward(hmm,observation)

print(exp(logBackwardProbabilities))

# Calculate forward probablities
logForwardProbabilities = forward(hmm,observation)
print(exp(logForwardProbabilities))


# Calculate posterior probablities of the states
posterior = posterior(hmm,observation)
print(posterior)


# Calculate Viterbi path
viterbi = viterbi(hmm,observation)
print(viterbi)



##########


#How to get DNA Sequences from GenBank
#one way ... 

#install.packages("ape")
library(ape)

ref <- c("U15717")

Rampho <- read.GenBank(ref, species.names = T,
                       gene.names = T, as.character = T)

attr(Rampho, "species")


cbind(attr(Rampho, "species"), names(Rampho))




bw.Rampho = baumWelch(fithmm,toupper(Rampho[[1]]))

print(bw.Rampho$hmm)
print(bw.Rampho$difference)


#log backward probabilities calc
logBackwardProbabilities = backward(hmm,toupper(Rampho[[1]]))

print(exp(logBackwardProbabilities))

# Calculate forward probablities
logForwardProbabilities = forward(hmm,toupper(Rampho[[1]]))
print(exp(logForwardProbabilities))


# Calculate posterior probablities of the states
posterior = posterior(hmm,toupper(Rampho[[1]]))
print(posterior)


# Calculate Viterbi path
viterbi = viterbi(hmm,toupper(Rampho[[1]]))
print(viterbi)



###############

#How to get protein/DNA squences
#Goal setup model for hydrophobic/phillic states or other examples


#install.packages("seqinr")
library("seqinr")

#load from file
#myFasta <- read.fasta(file = "myFasta.fasta")


choosebank("swissprot")
leprae <- query("leprae", "AC=Q9CD83")
lepraeseq <- getSequence(leprae$req[[1]])
closebank()
lepraeseq # Display the contents of "lepraeseq"



##################END CLASS DEMO###########################











