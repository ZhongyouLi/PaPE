# Source code for PaPE the Shiny app
PaPE is a user friendly Shiny app to perfome Parallel analysis of Pseudomonas aeruginosa & Escherichia coli auto-annotated 
transcriptomic compendia. The transcriptomic compendia was built using the automated sample group detection algorithm followed by ANOVA
analysis. Please refer the corresponding publication for further details.

The shiny app has deployed on the free shiny server, Shinyapps.io, at https://iamsoshiny.shinyapps.io/pape/ 
For users who prefer to run the app locally or modify the script, please follow the instruction below:

1. Clone this repository into your local directories
2. install required packages in your R environment:

```R
install.packages("dplyr")
install.packages("DT")
install.packages("shiny")
install.packages("shiyjs")
```
3. Launch PaPE the Shiny app
```R
library(shiny)
runApp(launch.browser = TRUE)
```
