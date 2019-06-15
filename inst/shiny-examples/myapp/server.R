#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(plotly)
library(DT)
#library(iAdapt)
library(here)


# load data? and source functions 
source(here::here("R/beta.ab.R"))
source(here::here("R/eff.stg1.R"))
source(here::here("R/next.dose.R"))
source(here::here("R/rand.stg2.R"))
source(here::here("R/safe.dose.R"))
source(here::here("R/sim.trials.R"))
source(here::here("R/tox.profile.R"))


#################################################
# Define server logic 
shinyServer(function(input, output, session) {
   
  ## example: set range of x based on user choice
 # output$season_choice <- renderUI({
    #seasons = nba_shots %>% filter(player_name == input$player_choice) %>% 
    #  distinct(season) %>% pull()
    
  #  selectizeInput("season_choice", label = "Select season", choices = seasons,
   #                selected = seasons[1], multiple = TRUE)
 # })
  
  ## example: subset data by selected player using reactive expression
 # player_data = reactive({
  #  filter(nba_shots, player_name == input$player_choice)
 # })
  
## Unacceptable minimum ##
  
output$Udlt = renderUI({
  p_no_min = input$p_yes
  
  sliderInput("p_no", h4("Unacceptable DLT:"),
              min = p_no_min, max = 1, value = 0.40,
              step = diff(0:1/20, 1) #animate=TRUE
  )
})


  
##Creating toxicity and efficacy vectors##
m <- eventReactive(input$update, {
  as.numeric(c(input$eff1, input$eff2, input$eff3, input$eff4, input$eff5, 
            input$eff6, input$eff7, input$eff8, input$eff9, input$eff10)[1:input$dose])
  })
  
output$eff_range <- renderPrint({ 
    m()
  })

dose.tox <- eventReactive(input$update, {
  as.numeric(c(input$tox1, input$tox2, input$tox3, input$tox4, input$tox5, 
                                   input$tox6, input$tox7, input$tox8, input$tox9, input$tox10)[1:input$dose])
})

output$tox_range <- renderPrint({ 
  dose.tox()
})

dose <- eventReactive(input$update, {
  input$dose
})
  
##Creating datatables##
seed <- eventReactive(input$update,{
  rpois(1, 100000)
})

d1 <- eventReactive(input$update, {
  set.seed(seed())
  tox.profile(dose = input$dose, 
              dose.tox = dose.tox(),
              p1 = input$p_no, p2 = input$p_yes, K = input$K, coh.size = input$coh.size) %>% 
    as_tibble() %>% 
    mutate(logical = ifelse(V4 > (1/input$K), TRUE, FALSE)) %>% 
    rename("Dose Assignment" = V1, "Dose Limiting Toxicity Events" = V2, "Cohort" = V3, "Likelihood Ratio of Toxicity: 1 trial" = V4)
  }, ignoreNULL = FALSE)

d2 <- eventReactive(input$repeated, {
  set.seed(seed())

  sim.trials(numsims = 100, #change this to select 
            dose = input$dose, 
            dose.tox = dose.tox(),
            p1 = input$p_no, p2 = input$p_yes, K = input$K, coh.size = input$coh.size,
            m = m(), v = rep(input$v, input$dose), N = 30)$safe.d %>% 
    colMeans(.) %>% 
    as_tibble() %>% 
    mutate(Dose = row_number(), value = 100*value) %>% 
    select(Dose, `Number of Trials for Which Dose was Safe` = value)

}, ignoreNULL = FALSE)

eff <- eventReactive(input$update, {
  set.seed(seed())
  eff.stg1(dose = input$dose, 
                dose.tox = dose.tox(),
                p1 = input$p_no, p2 = input$p_yes, K = input$K, coh.size = input$coh.size, 
                m = m(),
                v = rep(input$v, input$dose), nbb = 100)
}, ignoreNULL = FALSE)


##Table outputs##

output$table1 <- DT::renderDataTable({
    d1 <- d1() 
    ###Use package DT to make it nice##
    DT::datatable(d1, rownames = FALSE, 
                  options = list(paging = FALSE, searching = FALSE, rownames = FALSE,
                                 columnDefs = list(list(visible=FALSE, targets = 4))))  %>% 
      formatStyle("logical",
                    target = 'row',
                    backgroundColor = styleEqual(c(1, 0), c('#98fb98', '#ffcccb'))
      )
     })

output$table2 <- DT::renderDataTable({
    d2 <- d2()
    ###Use package DT to make it nice##
    DT::datatable(d2, rownames = FALSE, 
                  options = list(paging = FALSE, searching = FALSE))
    })

##Plot output##

output$plot <- renderPlot({
  p2 <- tibble(efficacy = m(),
               toxicity = dose.tox(),
               dose = seq(from = 1, to = dose(), by = 1)
               ) %>% 
    ggplot() +
    geom_line(aes(y = efficacy, x = dose, colour = "Efficacy")) +
    geom_point(aes(y = efficacy, x = dose, colour = "Efficacy"))
  
  # adding the relative humidity data, transformed to match roughly the range of the temperature
  p2 <- p2 + geom_line(aes(y = toxicity*100, x = dose, colour = "Toxicity")) +
             geom_point(aes(y = toxicity*100, x = dose, colour = "Toxicity"))
  
  # now adding the secondary axis, following the example in the help file ?scale_y_continuous
  # and, very important, reverting the above transformation
  p2 <- p2 + scale_y_continuous(sec.axis = sec_axis(~./100, name = "Toxicity"))
  
  p2 <- p2 + scale_colour_manual(values = c("blue", "red"))
  p2 <- p2 + labs(y = "Mean Efficacy",
                x = "Dose Level",
                colour = " ")
  p2 <- p2 + theme_bw() +
               theme(legend.position = c(0.8, 0.9),
                   text = element_text(size = 15))
  
  p2
  
               
})

output$plot_safe <- renderPlotly({
p1 <- tibble(eff_safe = eff()$Y.safe, 
         dose_safe = eff()$d.safe) %>% 
    ggplot(aes(y = eff_safe, x = factor(as.character(eff()$d.safe)))) + 
    geom_boxplot() +
    labs(y = "Efficacy Outcomes", x = "Safe Doses")

plotly::ggplotly(p1)
})

})
