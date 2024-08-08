#Michael Mascitti
#Dr. Leary 
#ABM Model

###################################################################################

#Setup Array

My.Variables = c("dice1","dice2","n-cell","n-cell-cn1","n-cell-cn2","n-cell-cn3",
                "ncc", "n-cell","n-cell-del","n-cell-ndel","T","P")

my.Cols = c("Cell.Number", "Y.Cord", "X.Cord","Type","Age","crbp", "egfr","psgs", "p53", "rb",
           "color","Roll1","Roll2","P","Treatment1", "Treatment2","Rep.Type", "Mutation Score",
           "Time.Since.LastSplit","Div.on.day", "Div.Dif")

N.cells = 100

Matrix.Setup = matrix(ncol = 21, nrow = N.cells)

Matrix.Setup

colnames(Matrix.Setup) = my.Cols

row.names(Matrix.Setup) = integer(N.cells)  

Cell.Array = as.data.frame(Matrix.Setup)

Cell.Array

###############################################################################

#Starting Values


for(i in 1:100){
  
  Cell.Array$Cell.Number[i] = i
  
  Cell.Array$Y.Cord[i] = 0
  
  Cell.Array$X.Cord[i] = i
  
  Cell.Array$Type[i] = "Normal"
  
  Cell.Array$Age[i] = 0
  
  Cell.Array$crbp[i] = 0 
  
  Cell.Array$egfr[i] = 0
  
  Cell.Array$psgs[i] = 0
  
  Cell.Array$p53[i] = 0
  
  Cell.Array$rb[i] = 0
  
  Cell.Array$color[i] = "Green"
  
  Cell.Array$Roll1[i] = sample(1:6,1)
  
  Cell.Array$Roll2[i]= sample(1:6, 1)
  
  Cell.Array$P[i] = Cell.Array$Roll1[i]+Cell.Array$Roll2[i]
  
  Cell.Array$Treatment1[i] = sample(1:2 , 1)
    
  Cell.Array$Treatment2[i] = sample(1:2, 1)
  
  if(Cell.Array$Treatment1[i] == 1){
    Cell.Array$Treatment1[i] = "True"
  }else{
    Cell.Array$Treatment1[i] = "False"
  }
  
  if(Cell.Array$Treatment2[i] == 1){
    Cell.Array$Treatment2[i] = "True"
  }else{
    Cell.Array$Treatment2[i] = "False"
  }
  
  T = 0
  
}

###############################################################################

#Day 1 

for(i in 1:100){
  
  Cell.Array$Age[i] = Cell.Array$Age[i] + 1
  
  Cell.Array$Type[i] = if(Cell.Array$P[i]>= 6){
    "CIN1"
  }else{
    "Normal"
  }
  
  Cell.Array$crbp[i] = 0 
  
  Cell.Array$egfr[i] = 0
  
  Cell.Array$psgs[i] = 0
  
  Cell.Array$p53[i] = 0
  
  Cell.Array$rb[i] = 0
  
  Cell.Array$color[i] = if(Cell.Array$Type[i] == "Normal"){
                           "Green"} else{ "Red"}
  T = 1
  
}


#############################################################################

#Adding Cells

add.cell <- function(current.array, row.number, T){ 
  
  
  this.cell <- current.array[row.number, ]
  this.X <- this.cell$X.Cord
  this.Y <- this.cell$Y.Cord
  current.array$Time.Since.LastSplit[row.number] <- T
  current.array$Div.Dif[row.number] <- if(current.array$Type[row.number] == "Normal"){
    3
  }else if(current.array$Type[row.number] == "CIN1"){
    3
  }else if(current.array$Type[row.number] == "CIN2"){
    2
  }else if(current.array$Type[row.number] == "CIN3" | current.array$Rep.Type[row.number] == "Damaged"){
    2
  }else if(current.array$Type[row.number] == "CNCR" | current.array$Rep.Type[row.number] == "Broken"){
    1
  }else{
    999
  }
  
  current.array$Div.on.day[row.number] <- current.array$Time.Since.LastSplit[row.number] + current.array$Div.Dif[row.number]
  
  #Moving Cells
  
  move.these.rows <- which((current.array$X.Cord == this.X) &
                             (current.array$Y.Cord > this.Y))
  for(i in move.these.rows){
    current.array$Y.Cord[i] <- current.array$Y.Cord[i] + 1
  }
  
  
  #New Cell, Add Columns
  
  new.cell <- this.cell
  
  new.cell$Cell.Number <- max(Cell.Array$Cell.Number) + 1
  
  new.cell$Age <- 0
  
  new.cell$Y.Cord <- new.cell$Y.Cord + 1
  
  new.cell$Roll1 = sample(1:6,1)
  
  new.cell$Roll2 = sample(1:6, 1)
  
  new.cell$P = new.cell$Roll1 + new.cell$Roll2 
  
  new.cell$Time.Since.LastSplit <- T
  
  new.cell$egfr <- this.cell$egfr
  new.cell$p53 <- this.cell$p53
  new.cell$rb <- this.cell$rb
  new.cell$crbp <- this.cell$crbp
  
  # Then add the new cell to current.array, and return:
  
  new.array <- rbind(current.array, new.cell)
  
  return(new.array)
  
}



###############################################################################

#Gene mutation

Gene.Mutations <- function(Cell.Array, row.number) {
  
  
    
  mutp53 <- sample(-2:1, 1, replace = TRUE) 
    
  mutcrbp <- sample(-2:1, 1, replace = TRUE) 
    
  mutrb <- sample(-2:1, 1, replace = TRUE) 
    
  mutegfr <- sample(-2:1, 1, replace = TRUE) 
    

  Cell.Array$p53[row.number] <- Cell.Array$p53[row.number] + mutp53
  
  Cell.Array$crbp[row.number] <- Cell.Array$crbp[row.number] + mutcrbp
  
  Cell.Array$rb[row.number] <- Cell.Array$rb[row.number] + mutrb
  
  Cell.Array$egfr[row.number] <- Cell.Array$egfr[row.number] + mutegfr

  
  return(Cell.Array)
  
}



#################################################################################

#Mutation Rate based on gene differences 

mutation.track <- function(Cell.Array, row.number) {
  
  Cell.Array$`Mutation Score`[row.number] = Cell.Array$egfr[row.number] + Cell.Array$rb[row.number] +
    Cell.Array$p53[row.number] + Cell.Array$crbp[row.number]
  
  Cell.Array$`Mutation Score`[is.na(Cell.Array$`Mutation Score`)] <- median(Cell.Array$`Mutation Score`, na.rm = TRUE) 

  # Return the modified row
  return(Cell.Array)
}

###############################################################################

#Mutation Rate based on gene differences 

Rep.function <- function(current.array, row.number){
  
 dat = current.array$`Mutation Score`[row.number]
  
  Broken <- if(dat > -30){
    0 
  }else if(dat <= -30 && dat > -70){
    1
  }else if(dat <= -70){
    2
  }else{
    777
  }
  
  a = runif(1)
  
 current.array$Rep.Type[row.number] = if(Broken == 0 && a >= 0.25){
    "Normal"
  }else if(Broken == 1 && a >= 0.50){
    "Damaged"
  }else if(Broken == 2 && a >= .25){
    "Broken"
  }else{
    "Normal"
  }

  # Return the modified row
  return(current.array)
}

###############################################################################

#Adding Gene information to each cell 

for(i in 1:nrow(Cell.Array)){
  
  if(Cell.Array$Age[i] >= 1){
    
    Cell.Array$crbp[i] = sample(1:2, 1) - sample(1:3, 1)
    Cell.Array$egfr[i] = sample(1:2, 1) - sample(1:3, 1)
    Cell.Array$psgs[i] = sample(1:2, 1) - sample(1:3, 1)
    Cell.Array$p53[i] = sample(1:2, 1) - sample(1:3, 1)
    Cell.Array$rb[i] = sample(1:2, 1) - sample(1:3, 1)
    
  }
  
  
  if(Cell.Array$Treatment1[i] == "True" & Cell.Array$Age[i] >= 6){
    
    Cell.Array$crbp[i] = sample(1:2, 1) - sample(1:3, 1)
    Cell.Array$egfr[i] = sample(1:2, 1) - sample(1:3, 1)
    Cell.Array$psgs[i] = sample(1:2, 1) - sample(1:3, 1)
    Cell.Array$p53[i] = sample(1:2, 1) - sample(1:3, 1)
    Cell.Array$rb[i] = sample(1:2, 1) - sample(1:3, 1)
    
  }
  
}



###############################################################################

#Setting variables for dice roll

Roll1= runif(nrow(Cell.Array), min=0, max = 0)    

Roll2= runif(nrow(Cell.Array), min=0, max = 0)

R1.R2= runif(nrow(Cell.Array), min=0, max = 0)

keep_rows <- integer()

T

for(i in 1:nrow(Cell.Array)){
  Cell.Array$Div.on.day[i] = T
  
}

##########################################################################

#Everyday (First run with 100 cells)


for(j in 1:50){  #8 - 9 days is the sweet spot 
  
  Current.State = runif(nrow(Cell.Array), min = 0, max = 0)
  Current.State = Cell.Array$Type                       
  
  
  
  for(i in 1:nrow(Cell.Array)){
    
    #Age and Dice Roll
    
    Cell.Array$Age[i] = Cell.Array$Age[i] + 1
    
    Roll1[i] = sample(1:6, 1)
    
    Roll2[i] = sample(1:6, 1)
    
    R1.R2[i] = Roll1[i] + Roll2[i]
    
 
      

    #Updating State
    
    if(Current.State[i] == "CIN3" && R1.R2[i] != 7){
        
        Cell.Array$Type[i] =  "CNCR"
        
    }
    
    
      
    if(Current.State[i] == "CIN2" && R1.R2[i] == 7 && Cell.Array$P[i] != 7 ){
        
        Cell.Array$Type[i] =  "CIN3"
        
      }
    
      
   if(Current.State[i] == "CIN1" && R1.R2[i] == 5  && Cell.Array$P[i] != 7 && Cell.Array$P[i] > 2){
        
        Cell.Array$Type[i] = "CIN2"
        
      }
    
   if(Current.State[i] == "Normal" && R1.R2[i] == 2 && Cell.Array$P[i] != 7 && Cell.Array$P[i] > 6){
     
        Cell.Array$Type[i] = "CIN1"
     
   }
    

    #Updating Color 
    
    Cell.Array$color[i] = if(Cell.Array$Type[i] == "Normal"){    
      "Green"
      }else if(Cell.Array$Type[i] == "CIN1"){
        "Blue"
      }else if(Cell.Array$Type[i] == "CIN2"){
        "Yellow"
      }else if(Cell.Array$Type[i] == "CIN3"){
        "Red"
      }else{
        "Black"
      }
    
    #Division 
    
    if(Cell.Array$Div.on.day[i] == T){                                  
      
      
      Cell.Array = add.cell(Cell.Array, i, T)      #Division
      
      Cell.Array = Gene.Mutations(Cell.Array, i)   #Gene Mutation
      
      
    }

    

   
    #Assigning Each cell Mutation Score 
    
    
    if(Cell.Array$Div.on.day[i] >= T){                                  
      
      
      Cell.Array = mutation.track(Cell.Array, i)
      
      
    }      
    
    
    #Non reproducing Mutation 
    
    
    a = runif(1)
    
    
    if(a >= .95){                                  
        
      
      Cell.Array = Gene.Mutations(Cell.Array, i)   #Gene Mutation
      
      
    }

      
    #Updating Reproduction Score 
    
   if(T > 0){
     
     Cell.Array = Rep.function(Cell.Array, i)
   
   }
    
    #Taking Away Old Cells 
    
    if(Cell.Array[i, "Age"] < 20){
      
      keep_rows <- c(keep_rows, i)
      
    }
 

    
  }
  
  T = 1 + j
  
  
  F.N.Cell = length(which(Cell.Array$Type== "Normal"))/nrow(Cell.Array)
  
  F.CN1.Cell = length(which(Cell.Array$Type == "CIN1"))/nrow(Cell.Array)
  
  F.CN2.Cell = length(which(Cell.Array$Type == "CIN2"))/nrow(Cell.Array)
  
  F.CN3.Cell = length(which(Cell.Array$Type == "CIN3"))/nrow(Cell.Array)
  
  F.CNCR.Cell = length(which(Cell.Array$Type == "CNCR"))/nrow(Cell.Array)
  
  
  
  
  if(j == 5 | j == 15 | j == 35 | j == 50){
  plot(Cell.Array$X.Cord,Cell.Array$Y.Cord, col = Cell.Array$color, pch = 16, ylim = c(0, 80), xlim = c(0,100), 
       ylab = "Y - Coordinate", xlab = "X - Coordinate", main = paste("Day", j))
    
    
  }
 
    
  if(j == 5 | j == 15 | j == 35 | j == 50){
  barplot(c(F.N.Cell, F.CN1.Cell, F.CN2.Cell, F.CN3.Cell, F.CNCR.Cell), names.arg = c("Normal", "CIN1", "CIN2", "CIN3","Cancer"),
          ylab = "Percent of Total Cells", ylim = c(0,1), col = "Black", main = paste("Day", j), xlab = "Cell Type")
    abline(h = 0)
  }
}

Cell.Array

tail(Cell.Array)
###############################################################################

#Information

length(which(Cell.Array$Type== "Normal"))

length(which(Cell.Array$Type == "CIN1"))

length(which(Cell.Array$Type == "CIN2"))

length(which(Cell.Array$Type == "CIN3"))

length(which(Cell.Array$Type == "CNCR"))

length(which(Cell.Array$Rep.Type== "Normal"))

length(which(Cell.Array$Rep.Type== "Broken"))

length(which(Cell.Array$Rep.Type== "Damaged"))

nrow(Cell.Array)

###########################################################################

#Hpv causes a orgin spot of cancer
#Manipulation of genes can cause different growth rates
#Disease hot map

CRBP - COntrols cell growth and matures the cell making it specialized 
EGFR - Growth factor, responds to enviroment 
P53 - Tumor supressor 
Rb - G1 to S phase, when in contact with E7 from HPV cell profeliaration is induced 
  


#Feb 23rd, to March 3rd 

  
Ratio of cell divison rate between cancer and normal cells (Not able to Complete)
 - Use article and talk about it
 - Gompertz equation? 
 - Refer to KI-67 test 

Add reproduction score and mutation rate AND Time since last split day column and when they are due to divide 
 - added columns 

Provide a mechanism for changing the genes values 


#March 3rd to March 9th 

make a mechanism with made up numbers that show our process for the ratio of cell divsion rate

Value of mutation rate based off gene values 

Please find the different reproduction rates of the different stages 

Update mutation rate, time since last split, dividing in T 


#March 23rd to March 30th 

whats the reproduction rate of cell i input type of cell and number of mutations 

fix abs and have it return mutation score #Mutation Score now assigned to each cell

write a function to tell us when cell is going tor reproduce #done 
-depend on cell type and mutations 

Email Leary this weekend with an update
Over the next couple weeks get together more often (Tuesday 28th 10 am)

#March 28th 

Fix up the time function and the everyday function so that the div on day and ts are changing #Done
Take mutation score and help it modify reproduction rate  #Done but highlighted out
change mutation score to walk and checkpoint style  #Done'


#### 

change the div dif to subtract based off reproduction values, do something with the 999 (done)
Add visuals for explanation 
Fix NA values (done)

#############################################################################

#Plot Practice 

library(ggplot2)
library(grid)

# Create the main scatterplot
p1 <- ggplot(mpg, aes(x = displ, y = hwy)) + geom_point()

# Create the preview scatterplot
p2 <- ggplot(mpg, aes(x = cyl, y = cty)) + geom_point()

# Convert the preview plot to a grob object
grob <- ggplotGrob(p2)

# Create a viewport for the preview within the main plot
vp <- viewport(x = 0.7, y = 0.75, width = 0.2, height = 10)

# Draw the main scatterplot
p1 + 
  
  # Add the preview viewport
  annotation_custom(grob = grob, xmin = 4.5, xmax = 8, ymin = 12, ymax = 26) + 
  
  # Add labels for the axes and the preview
  labs(x = "Engine Displacement (L)", y = "Highway MPG",
       title = "Main Scatterplot with Preview") + 
  annotate("text", x = 6, y = 35, label = "Preview") + 
  annotate("segment", x = 7, xend = 7.5, y = 31, yend = 23, size = 0.5, arrow = arrow(length = unit(0.3, "cm")))






































