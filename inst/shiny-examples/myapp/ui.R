#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(plotly)
library(DT)
library(shinydashboard)

# example: data management
#load("x")
#players = x %>% distinct(y) %>% pull()
#made = x %>% distinct(t ) %>% pull()


# Define UI 
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Simulation"),
  
  # Sidebar  
  sidebarLayout(
    sidebarPanel(
       numericInput("dose",
                   h3("Number of pre-specified dose levels:"),
                   min = 2,
                   max = 10,
                   value = 5),
       #Acceptable (p_yes) and unacceptable (p_no) DLT rates used for establishing safety (select range between 0 and 1)
       fluidRow(
         box(width = 11, title = "Dose Limiting Toxicity (DLT) Rates",
             splitLayout(
              sliderInput("p_yes", h4("Acceptable DLT:"),
                   min = 0, max = 1, value = 0.15,
                   step = diff(0:1/20, 1) #animate=TRUE
                ),
              uiOutput("Udlt")
                      )
            )
         ),
       #next: Vector of true toxicities associated with each dose (select range and "by")
       #e.g. dose.tox <- c(0.05, 0.10, 0.15, 0.20, 0.30)       
       
       fluidRow(
         box(width = 11, title = "True Toxicities", 
             splitLayout(
               textInput("tox1", "1", value = 0.05),
               conditionalPanel(
                 condition = "input.dose >= 2",
                 textInput("tox2", "2", value = 0.08)
               ),
               conditionalPanel(
                 condition = "input.dose >= 3",
                 textInput("tox3", "3", value = 0.1)
               ),
               conditionalPanel(
                 condition = "input.dose >= 4",
                 textInput("tox4", "4", value = 0.2)
               ),
               conditionalPanel(
                 condition = "input.dose >= 5",
                 textInput("tox5", "5", value = 0.3)
               ),
               conditionalPanel(
                 condition = "input.dose >= 6",
                 textInput("tox6", "6", value = 0.4)
               ),
               conditionalPanel(
                 condition = "input.dose >= 7",
                 textInput("tox7", "7", value = 0.5)
               ),
               conditionalPanel(
                 condition = "input.dose >= 8",
                 textInput("tox8", "8", value = 0.6)
               ),
               conditionalPanel(
                 condition = "input.dose >= 9",
                 textInput("tox9", "9", value = 0.7)
               ),
               conditionalPanel(
                 condition = "input.dose == 10",
                 textInput("tox10", "10", value = 0.8)
               )
               
             )
         )
       ),
       #Likelihood-ratio (LR) threshold (2, 4, 8?  check paper)
       numericInput("K",
                    h3("Likelihood Ratio Threshold:"),
                    min = 1,
                    max = 32,
                    step = 1,
                    value = 2),
       #Cohort size used in stage 1: Usually always 3, right?
       numericInput("coh.size",
                    h3("Stage 1 cohort size:"),
                    min = 1,
                    max = 10,
                    step = 1,
                    value = 3),
       #Vector of true mean efficacies per dose (select range and "by")
       fluidRow(
         box(width = 11, title = "True Mean Efficacies", 
             splitLayout(
               textInput("eff1", "1", value = 5),
               conditionalPanel(
                 condition = "input.dose >= 2",
                 textInput("eff2", "2", value = 10)
                ),
               conditionalPanel(
                 condition = "input.dose >= 3",
                 textInput("eff3", "3", value = 15)
                ),
               conditionalPanel(
                 condition = "input.dose >= 4",
                 textInput("eff4", "4", value = 20)
               ),
               conditionalPanel(
                 condition = "input.dose >= 5",
                 textInput("eff5", "5", value = 25)
               ),
               conditionalPanel(
                 condition = "input.dose >= 6",
                 textInput("eff6", "6", value = 30)
               ),
               conditionalPanel(
                 condition = "input.dose >= 7",
                 textInput("eff7", "7", value = 35)
               ),
               conditionalPanel(
                 condition = "input.dose >= 8",
                 textInput("eff8", "8", value = 40)
               ),
               conditionalPanel(
                 condition = "input.dose >= 9",
                 textInput("eff9", "9", value = 45)
               ),
               conditionalPanel(
                 condition = "input.dose == 10",
                 textInput("eff10", "10", value = 50)
               )
               
             )
         )
       ),
       
       numericInput("v", h3("True Variance of Efficacies:"),
                   min = 0, max = 1, value = 0.01,
                   step = 0.01),
       helpText("Note: Variance is assumed uniform"),
      #Update Button
      actionButton("update", "Simulate")
    ),
    
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Stage 1: Establishing Safety and Estimating Efficacy",
      fluidRow(
        box(title = h4("Selected Efficacy and Toxicities"),
                 plotOutput("plot")),
        box(title = h4("Simulated Toxicity Profile"), 
                   DT::dataTableOutput("table1"))
               ),
      fluidRow(
        box(title = h4("Simulated Safe Doses"),
                    DT::dataTableOutput("table2")),
        box(title = h4("Estimated Efficacy of Safe Doses"),
            plotlyOutput("plot_safe"))
        )
                  ),
      tabPanel("Stage 2: Adaptive Randomization"),
      tabPanel("Repeated Simulation")

    )
  )
)))
