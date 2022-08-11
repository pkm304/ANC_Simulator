## library
library(shiny)
library(ggplot2)
library(dplyr)
library(deSolve)
library(numDeriv)
library(reshape2)
library(mvtnorm)
library(adaptMCMC)
library(matrixcalc)


dosechemUI <- function(id, base_bsa, count) {
  ns <- NS(id)
  fluidRow(style='padding:20px;',
    fluidRow(
      column(6, 
             p(strong(paste("Dose ", count,":"))))
    ),
    fluidRow(
      column(6, 
             numericInput(ns("chem_time"), "Days", value = (count-1)*30 ), min = 0, max = 180, step = 1),
      column(6, 
             # sliderInput(ns("bsa_at_dose"), "BSA at dose (m^2)",
             #             min = 0.3, max = 4,
             #             value = base_bsa, step = 0.01)
             radioButtons(ns("BSA_input_at_dose"), "BSA at dose (m^2)",
                                           c("Direct input" = "direct",
                                             "Compute from height & weight" = "height.weight")),
             uiOutput(ns("RUI_bsa_1_at_dose"))
      )
    ),
    fluidRow(
      column(6,
             sliderInput(ns("dose_per_bsa_pac"), "Paclitaxel dose per BSA (mg/m^2)",
                         min = 100, max = 250,
                         value = 175, step = 1)
      ),
      column(6,
             textOutput(ns("dose_pac"))
      )
    ),
    fluidRow(
      column(6,
             sliderInput(ns("dose_per_bsa_cis"), "Cisplatin dose per BSA (mg/m^2)",
                         min = 50, max = 150,
                         value = 75, step = 1)
      ),
      column(6,
             textOutput(ns("dose_cis"))
      )
    )
  )
  
  
} 




doseGCSFUI <- function(id, base_bsa, count) {
  ns <- NS(id)
  fluidRow(style='padding:20px;',
    fluidRow(
      column(6, 
             p(strong(paste("Dose ", count,":")))),
      column(6, 
             numericInput(ns("gcsf_time"), "Days", value = (count-1)*30+15 ), min = 0, max = 180, step = 1)
      
    )
    
  )
  
  
} 


dosechemServer <- function(input, output, session, values, base_bsa, height, weight) {
  ns <- session$ns
  
 # values_at_dose <- reactiveValues()

  
  
  output$RUI_bsa_1_at_dose <- renderUI({
    ns <- session$ns
    #if(!is.null(input$BSA_input_at_dose)){
      if(input$BSA_input_at_dose == "direct"){
        sliderInput(ns("bsa_at_dose"), "Direct input",
                    min = 0.3, max = 5,
                    value = base_bsa, step = 0.01)
      }else{
        if(is.null(height)){
          fluidRow(
            column(5, numericInput(ns("height_at_dose"), "Height (cm)", value = 170, min = 100, max = 230, step = 0.5)),
            column(5, numericInput(ns("weight_at_dose"), "Weight (kg)", value = 70, min = 10, max = 130, step = 0.1)),
            column(2, h5("BSA"),textOutput(ns("bsa_hw_at_dose")))  
          )
        }else{
          fluidRow(
            column(5, numericInput(ns("height_at_dose"), "Height (cm)", value = height, min = 100, max = 230, step = 0.5)),
            column(5, numericInput(ns("weight_at_dose"), "Weight (kg)", value = weight, min = 10, max = 130, step = 0.1)),
            column(2, h5("BSA"),textOutput(ns("bsa_hw_at_dose")))
            
          )
        }
      }
    #}

    })
  
  observe({
      if(!is.null(input$BSA_input_at_dose)){
        if(input$BSA_input_at_dose == "direct"){
          values[[paste0(ns("bsa_at_dose"))]] <-     input$bsa_at_dose
        }else{
          values[[paste0(ns("bsa_at_dose"))]] <-     signif(sqrt(input$weight_at_dose*input$height_at_dose/3600),3)
        }
      }

    })

    output$bsa_hw_at_dose <- renderText({
      paste0( values[[paste0(ns("bsa_at_dose"))]])
      #1
    })

  
  
  
  output$dose_pac <- renderText({paste0(values[[paste0(ns("bsa_at_dose"))]]*input$dose_per_bsa_pac)})
  output$dose_cis <- renderText({paste0(values[[paste0(ns("bsa_at_dose"))]]*input$dose_per_bsa_cis)})
  
  #return(values_at_dose)
}


shinyServer(function(input, output, session) {
  values <- reactiveValues(upload_state = 'none')
  dose_time <- reactiveVal(list())
  
  observeEvent(input$dvdata, {
    if(!is.null(input$dvdata)){
      values$upload_state <- 'uploaded'
      values$dv <- read.csv(as.character(input$dvdata$datapath), header=T)
      colnames(values$dv) <- c('time','y')
      
      
      # values$dv_sample <- read.csv(as.character(input$dvdata$datapath), header=T)
      # colnames(values$dv_sample) <- c('time','y')
      # #tt <- dv$time
    }
  })
  
  
  observe({
    # retreive dose specification

    vars.chemdose <- names(input)[grepl("chemdose",names(input))]
    vars.chemdose <- vars.chemdose[!grepl("input", vars.chemdose)]
    vars.chemdose_bsa <- names(values)[grepl("chemdose",names(values))]
    
    vars.gcsfdose <- names(input)[grepl("gcsfdose",names(input))]
    vals.chemdose <-sapply(vars.chemdose, function(x) input[[x]])
    vals.chemdose_bsa <-sapply(vars.chemdose_bsa, function(x) values[[x]])
    vals.gcsfdose <-sapply(vars.gcsfdose, function(x) input[[x]])
    num_chemdose <- input$add_chem_dose
    num_gcsfdose <- input$add_gcsf_dose

    chemdose <- data.frame(time = rep(NA, num_chemdose), dosep = rep(NA, num_chemdose), dosec = rep(NA, num_chemdose))
    gcsfdose <- data.frame(time = rep(NA, num_gcsfdose))

    if(num_chemdose != 0 & (length(vars.chemdose) !=0 & length(vars.chemdose_bsa) != 0 & length(vars.chemdose[!grepl("t_at",vars.chemdose)])/4 == num_chemdose)){
      for(i in 1:num_chemdose){
        if(!is.null(vals.chemdose_bsa[[paste0("chemdose_", i,"-bsa_at_dose")]]) & length(vals.chemdose_bsa[[paste0("chemdose_", i,"-bsa_at_dose")]]) != 0){
          chemdose$time[i] <- vals.chemdose[[paste0("chemdose_", i,"-chem_time")]]
          chemdose$dosep[i] <- vals.chemdose_bsa[[paste0("chemdose_", i,"-bsa_at_dose")]]*vals.chemdose[[paste0("chemdose_", i,"-dose_per_bsa_pac")]]
          chemdose$dosec[i] <- vals.chemdose_bsa[[paste0("chemdose_", i,"-bsa_at_dose")]]*vals.chemdose[[paste0("chemdose_", i,"-dose_per_bsa_cis")]]  
        }
        
      }
    }

    if(num_gcsfdose != 0 & (length(vars.gcsfdose) !=0 & length(vars.gcsfdose) == num_gcsfdose)){
      for(i in 1:num_gcsfdose){
        gcsfdose$time[i] <- vals.gcsfdose[[paste0("gcsfdose_", i,"-gcsf_time")]]

      }
    }


    if(num_chemdose != 0 & (length(vars.chemdose) !=0 & length(vars.chemdose_bsa) != 0 & length(vars.chemdose[!grepl("t_at",vars.chemdose)])/4 == num_chemdose)) {
      chemdose.melt <-chemdose %>% melt(value.name = "dose" ,id.vars = "time")
      chemdose.melt <- chemdose.melt %>% mutate(time = if_else(variable == "dosec", as.double(time + 1), as.double(time)))
    }else{
      chemdose.melt <- data.frame(time = 0, variable = "no", dose = 0)
    }
    if(num_gcsfdose != 0 & (length(vars.gcsfdose) !=0 & length(vars.gcsfdose) == num_gcsfdose)){
      totaldose <- rbind(chemdose.melt, data.frame(time =gcsfdose$time, variable = "doseg", dose =0))
    }else{
      totaldose <- chemdose.melt
    }

    totaldose <- totaldose %>% arrange(time)
    values$totaldose <-totaldose
  })

  observeEvent(input$estim_prm_individual,{

    print(input$estim_prm_individual)
    values$estimated.indiv.prms <- runmain_indiv_estim(totaldose =  values$totaldose,
                                                 SEX = input$sex,
                                                 DM = input$dm,
                                                 maxtime = input$TIME,
                                                 BASE = input$BASE,
                                                 DV = values$dv,
                                                 estim.variability = input$estim_indiv_var)
    values$estimated.indiv.prms$totaldose <- values$totaldose
    values$estimated.indiv.prms$dv <- values$dv
    values$estimated.indiv.prms$demo <- data.frame(BSA = values$bsa, ANC = input$BASE,Sex = if_else(input$sex == "1", "Female", "Male"), DM = if_else(input$dm == "1", "Yes", "No")  )


  })
  
  

  observeEvent(input$reset, {
    values$upload_state <- 'none'
    values$dv <- NULL
    values$estimated.indiv.prms <- NULL
  })
  
  
  observeEvent(input$add_chem_dose, {
    new_id <- paste("chemdose", input$add_chem_dose, sep = "_")
    if(input$add_chem_dose <= 6){
    insertUI(
      selector = "#placeholderchem",
      where = "beforeBegin",
      ui = dosechemUI(new_id, values$bsa, as.numeric(input$add_chem_dose))
    )
    
    module.out <- callModule(dosechemServer, new_id, values, values$bsa, input$height, input$weight)
    
    #print(module.out)
    #print(new_dose)
    #dose_time(c(dose_time, new_dose))
    }
  })
  
  observeEvent(input$add_gcsf_dose, {
    new_id <- paste("gcsfdose", input$add_gcsf_dose, sep = "_")
    if(input$add_gcsf_dose <= 10){
      insertUI(
        selector = "#placeholdergcsf",
        where = "beforeBegin",
        ui = doseGCSFUI(new_id, values$bsa, as.numeric(input$add_gcsf_dose))
      )
      
      #callModule(doseServer, new_id)
      
      #print(new_dose)
      #dose_time(c(dose_time, new_dose))
    }
  })
  
  
  observe({
    print(dose_time)
  })
  
  # output$dosep <- renderText({ paste0("Paclitaxel: ",input$bsa*175, "mg") })
  # output$dosec <- renderText({ paste0("Cisplatin: ", input$bsa*75, "mg") })
  output$estim_status <- renderText({
    if(!is.null(values$estimated.indiv.prms)){
      "MAP parameters are estimated based on:"
    }else{
      "No estimated MAP parameters yet..."
    }
  })
  
  output$estim_status_1 <- renderText({
    if(!is.null(values$estimated.indiv.prms)){
      "MAP parameters are estimated based on:"
    }else{
      "No estimated MAP parameters yet..."
    }
  })
  
  output$RUI_estim_prm_individual <- renderUI({
    if(!is.null(values$dv)){
      actionButton("estim_prm_individual", strong("Estimate MAP individual parameters"))
    }
  })
  
  output$RUI_estim_indiv_var <- renderUI({
    if(!is.null(values$dv)){
      radioButtons(inputId = "estim_indiv_var",
                   label="With posterior distribution",
                   choices = list("No"=0,
                                  "Yes (may take long..)"=1))
    }
  })

  
  
  output$RUI_plot_individual <- renderUI({
    if(!is.null(values$estimated.indiv.prms)){
      checkboxInput("plot_individual", strong("MAP individual trajectory"), value = TRUE)
    }
  })
  
  output$RUI_plot_indiv_var <- renderUI({
    if(!is.null(values$estimated.indiv.prms)){
      if(!is.null(values$estimated.indiv.prms$mat.resampled)){
        radioButtons(inputId = "indiv_var",
                     label="With posterior uncertainties",
                     choices = list("No"=0,
                                    "Yes (may take long..)"=1))
      }
      
    }
  })
  
  
  output$RUI_bsa <- renderUI({
    # if(!is.null(input$`chemdose_1-bsa_at_dose`)){
    #   tagList(sliderInput("bsa", "Baseline body surface area (m^2)",
    #                       min = 0.3, max = 5,
    #                       value = input$`chemdose_1-bsa_at_dose`, step = 0.01),
    #           numericInput("height", "Input height (cm)", value = 170, min = 100, max = 230, step = 0.5),
    #           numericInput("weight", "Input weight (kg)", value = 70, min = 10, max = 130, step = 0.1))
    #   
    # }else{
    tagList(
      radioButtons("BSA_input", "Baseline body surface area (m^2)",
                   c("Direct input" = "direct",
                     "Compute from height & weight" = "height.weight")),
      uiOutput("RUI_bsa_1")
      
      
    )
      
    # }
      
  })
  
  output$RUI_bsa_1 <- renderUI({
    if(input$BSA_input == "direct"){
      sliderInput("bsa", "Direct input",
                  min = 0.3, max = 5,
                  value = 1.7, step = 0.01)
    }else{
      fluidRow(
        column(4, numericInput("height", "Height (cm)", value = 170, min = 100, max = 230, step = 0.5)),
        column(4, numericInput("weight", "Weight (kg)", value = 70, min = 10, max = 130, step = 0.1)),
        column(4, h6("Computed BSA (m^2)"),textOutput("bsa_hw"))
      )
    }
    
  })
  
  observe({
    if(!is.null(input$BSA_input)){
      if(input$BSA_input == "direct"){
        values$bsa <-     input$bsa
      }else{
        values$bsa <-     signif(sqrt(input$weight*input$height/3600),3)
      }
    }
    
  })
  
  output$bsa_hw <- renderText({
    paste0( values$bsa)
  })
  
  output$RUI_estim_smpl_dose <- renderUI({
    if(!is.null(values$estimated.indiv.prms)){
      tabsetPanel( #L5
        tabPanel('Demographics', #L6 
                 tableOutput(outputId = "tab_estim_demo")
        ),
        tabPanel('Dosing', #L6 
                 tableOutput(outputId = "tab_estim_dose")
                 
        ), #L6 closed
        tabPanel('Sample', #L6
                 tableOutput(outputId = "tab_estim_smpl")
        )
      )#L6 closed    
    }
  })
  
  output$tab_estim_demo <- renderTable({
    values$estimated.indiv.prms$demo
  })
  
  output$tab_estim_dose <- renderTable({
    df <- values$estimated.indiv.prms$totaldose
    colnames(df) <- c("Time (Day)", "Drug", "Amount (mg)")
    df <- df %>%  mutate(Drug = if_else(Drug == "dosec", "Cisplatin", if_else(Drug == "dosep", "Paclitaxel",if_else(Drug == "doseg", "GCSF", "No") )))  
    df   
  })
  
  output$tab_estim_smpl <- renderTable({
    df <- values$estimated.indiv.prms$dv
    colnames(df) <- c("Time (Day)", "ANC (10^6/mL)")
    df
  })
  
  
  output$RUI_estim_smpl_dose_1 <- renderUI({
    if(!is.null(values$estimated.indiv.prms)){
      tabsetPanel( #L5
        tabPanel('Demographics', #L6 
                 tableOutput(outputId = "tab_estim_demo_1")
        ),
        tabPanel('Dosing', #L6 
                 tableOutput(outputId = "tab_estim_dose_1")
                 
        ), #L6 closed
        tabPanel('Sample', #L6
                 tableOutput(outputId = "tab_estim_smpl_1")
        )
      )#L6 closed    
    }
  })
  
  
  output$tab_estim_demo_1<- renderTable({
    values$estimated.indiv.prms$demo
  })
  
  
  output$tab_estim_dose_1 <- renderTable({
    df <- values$estimated.indiv.prms$totaldose
    colnames(df) <- c("Time (Day)", "Drug", "Amount (mg)")
    df <- df %>%  mutate(Drug = if_else(Drug == "dosec", "Cisplatin", if_else(Drug == "dosep", "Paclitaxel",if_else(Drug == "doseg", "GCSF", "No") )))  
    df   
  })
  
  output$tab_estim_smpl_1 <- renderTable({
    df <- values$estimated.indiv.prms$dv
    colnames(df) <- c("Time (Day)", "ANC (10^6/mL)")
    df
  
  })
  
  
  output$simPlot <- renderPlot({
    input$gen_plot
    #input$plot_typical
    #input$plot_individual
    #input$pop_var
    #input$indiv_var
    #input$estim_indiv_var
    
    plot.options <- NULL # typical, indiv, typical_indiv, no_plot
    
    isolate({
      sim_time <-  input$TIME
      maxy <- input$maxy 
    })
    
    
    
    
     if(isolate( input$plot_typical & is.null(input$plot_individual)) ){
       plot.options <- "typical"
     }else if(isolate( input$plot_typical & isolate(!(input$plot_individual) ))){
       plot.options <- "typical"
     }else if(isolate( !(input$plot_typical)) & isolate(is.null(input$plot_individual))){
       plot.options <- "no_plot"
     }else if(isolate( !(input$plot_typical)) & isolate(!(input$plot_individual))){
       plot.options <- "no_plot"
     }else if(isolate( !(input$plot_typical)) & isolate(input$plot_individual)){
       plot.options <- "indiv"
     }else if(isolate( input$plot_typical) & isolate(input$plot_individual)){
       plot.options <- "typical_indiv"
     }
      
      
    if(plot.options == "typical" | plot.options == "typical_indiv"){
      isolate({
        values$pred_typical <- runmain_typical(totaldose =  values$totaldose, 
                                               SEX = input$sex, 
                                               DM = input$dm, 
                                               maxtime = input$TIME,
                                               BASE = input$BASE,
                                               plot.variability = input$pop_var)
        df.typical <-  values$pred_typical
      })
    }
      
    if(plot.options == "indiv" | plot.options == "typical_indiv"){
      isolate({
        if(is.null(input$indiv_var)){
          values$pred_indiv <- runmain_indiv_pred(totaldose =  values$totaldose,
                                                  SEX = input$sex, 
                                                  DM = input$dm, 
                                                  maxtime = input$TIME,
                                                  BASE = input$BASE,
                                                  DV = dv,
                                                  list.estimated.prms = values$estimated.indiv.prms,
                                                  plot.variability = "0")
        }else{
          values$pred_indiv <- runmain_indiv_pred(totaldose =  values$totaldose,
                                                  SEX = input$sex, 
                                                  DM = input$dm, 
                                                  maxtime = input$TIME,
                                                  BASE = input$BASE,
                                                  DV = dv,
                                                  list.estimated.prms = values$estimated.indiv.prms,
                                                  plot.variability = input$estim_indiv_var)
        }
        
        df.indiv <- values$pred_indiv
      })
    }
    
    if(plot.options == "indiv" | plot.options == "typical_indiv"){
      dv <- isolate(values$dv)
    }
    
    colors.all <- c("Typical traj." = "black", 
                "POV of indiv. traj." = "deepskyblue", 
                "95% CI of indiv. traj." = "deepskyblue4", 
                "MAP Indiv. traj." = "blue", 
                "Observations" = "black",
                "IIV of typical traj." = "gray", 
                "95% CI of typical traj." = "gray37")
    
    
    
    gg <- ggplot()
    
    if(plot.options == "typical"){
      colors <- colors.all[1]
      
      if(isolate(input$pop_var) == "1"){
        df.typical.melt <- df.typical[,-c(2:8)] %>% melt(value.name = "anc" ,id.vars = "time")
        df.typical.melt.ci <- df.typical.melt  %>%  group_by(time) %>% summarise(ci.l = quantile(anc, 0.025 ), ci.u = quantile(anc,0.975), ci.med = quantile(anc,0.5))
        gg <- gg + geom_line(data = df.typical.melt, mapping =  aes(x=time,y=anc, group = variable, color = "IIV of typical traj."), alpha = 0.1) +
          geom_line(data = df.typical.melt.ci, mapping=  aes(x=time,y=ci.l, color = "95% CI of typical traj."), alpha = 0.5, size = 1)+
          geom_line(data = df.typical.melt.ci, mapping=  aes(x=time,y=ci.u, color = "95% CI of typical traj."), alpha = 0.5, size = 1)+
          geom_line(data = df.typical.melt.ci, mapping=  aes(x=time,y=ci.med, color = "95% CI of typical traj."), linetype = "dashed", alpha = 0.5, size = 1)
        colors <- colors.all[c(1,6,7)]
      }
      
      gg <- gg + geom_line(data = df.typical, mapping = aes(x=time,y=typical, color = "Typical traj."), size=1)
      
    }else if (plot.options == "typical_indiv"){
     
      
      if(isolate(input$pop_var) == "1"){
        df.typical.melt <- df.typical[,-c(2:8)] %>% melt(value.name = "anc" ,id.vars = "time")
        df.typical.melt.ci <- df.typical.melt  %>%  group_by(time) %>% summarise(ci.l = quantile(anc, 0.025 ), ci.u = quantile(anc,0.975), ci.med = quantile(anc,0.5))
        gg <- gg + geom_line(data = df.typical.melt, mapping =  aes(x=time,y=anc, group = variable, color = "IIV of typical traj."), alpha = 0.1) +
          geom_line(data = df.typical.melt.ci, mapping=  aes(x=time,y=ci.l, color = "95% CI of typical traj."), alpha = 0.5, size = 1)+
          geom_line(data = df.typical.melt.ci, mapping=  aes(x=time,y=ci.u, color = "95% CI of typical traj."), alpha = 0.5, size = 1)+
          geom_line(data = df.typical.melt.ci, mapping=  aes(x=time,y=ci.med, color = "95% CI of typical traj."), linetype = "dashed", alpha = 0.5, size = 1)
      }
      
      if(!is.null(isolate(input$indiv_var))){
        if(isolate(input$indiv_var) == "1"){
          df.indiv.melt <- df.indiv[,-c(2:8)] %>% melt(value.name = "anc" ,id.vars = "time")
          df.indiv.melt.ci <- df.indiv.melt %>%   group_by(time) %>% summarise(ci.l = quantile(anc, 0.025 ), ci.u = quantile(anc,0.975), ci.med = quantile(anc,0.5))
          
          gg <- gg + geom_line(data = df.indiv.melt , mapping =  aes(x=time,y=anc, group = variable, color = "POV of indiv. traj."), alpha = 0.1) +
            geom_line(data = df.indiv.melt.ci, mapping=  aes(x=time,y=ci.l, color = "95% CI of indiv. traj."), alpha = 0.5, size = 1)+
            geom_line(data = df.indiv.melt.ci, mapping=  aes(x=time,y=ci.u, color = "95% CI of indiv. traj."), alpha = 0.5, size = 1)+
            geom_line(data = df.indiv.melt.ci, mapping=  aes(x=time,y=ci.med, color = "95% CI of indiv. traj."), linetype = "dashed", alpha = 0.5, size = 1)
          
        }
      }
          
      gg <- gg + geom_line(data = df.typical, mapping = aes(x=time,y=typical, color = "Typical traj."), size=1) +
        geom_line(data = df.indiv, mapping=  aes(x=time,y=indiv, color = "MAP Indiv. traj."), size = 1)+
        geom_point(data = dv, mapping=  aes(x=time,y=y, color = "Observations"), size = 1)
  
      
    }else if(plot.options == "indiv"){
      if(!is.null(isolate(input$indiv_var))){
        if(isolate(input$indiv_var) == "1"){
          df.indiv.melt <- df.indiv[,-c(2:8)] %>% melt(value.name = "anc" ,id.vars = "time")
          df.indiv.melt.ci <- df.indiv.melt %>%   group_by(time) %>% summarise(ci.l = quantile(anc, 0.025 ), ci.u = quantile(anc,0.975), ci.med = quantile(anc,0.5))
          
          gg <- gg + geom_line(data = df.indiv.melt , mapping =  aes(x=time,y=anc, group = variable, color = "POV of indiv. traj."), alpha = 0.1) +
            geom_line(data = df.indiv.melt.ci, mapping=  aes(x=time,y=ci.l, color = "95% CI of indiv. traj."), alpha = 0.5, size = 1)+
            geom_line(data = df.indiv.melt.ci, mapping=  aes(x=time,y=ci.u, color = "95% CI of indiv. traj."), alpha = 0.5, size = 1)+
            geom_line(data = df.indiv.melt.ci, mapping=  aes(x=time,y=ci.med, color = "95% CI of indiv. traj."), linetype = "dashed", alpha = 0.5, size = 1)
        }
      }
      
      gg <- gg + geom_line(data = df.indiv, mapping=  aes(x=time,y=indiv, color = "MAP Indiv. traj."), size = 1)+
        geom_point(data = dv, mapping=  aes(x=time,y=y, color = "Observations"), size = 1)

    }
    
    
    colors <- NULL
    if(exists("dv")){
      colors <- append(colors, colors.all[5])
    }
    if(exists("df.typical")){
      colors <- append(colors, colors.all[1])
    }
    if(exists("df.typical.melt.ci")){
      colors <- append(colors, colors.all[c(6,7)])
    }
    if(exists("df.indiv")){
      colors <- append(colors, colors.all[4])
    }
    if(exists("df.indiv.melt.ci")){
      colors <- append(colors, colors.all[c(2,3)])
    }
            
    
    
    gg <- gg + ylab(paste("Absolute Neutrophil Count (10^6/mL)",sep="")) +
      xlab("Time after dose (Day)")+
      theme_bw()+
      theme(legend.position = c(0.75, 0.85)) +
      labs(color = "Legend") +
      scale_color_manual(values = colors)+
      geom_hline(yintercept=0.5, linetype="dashed", color = "red") +
      ylim(0,maxy) +
      xlim(0, sim_time)
      
    
    gg
    
    
    
    
   
  })
 
  
  output$dummy <- renderTable({
    if(is.null(values$dv)){
      dummytab <- data.frame(t = c(0,10,28,40),
                             y = c(4.721, 1.288, 5.380, 0.179))
      colnames(dummytab) <- c('Time (Day)','ANC (10^6/mL)')
      dummytab
    }else{
      dummytab <- values$dv
      colnames(dummytab) <- c('Time (Day)','ANC (10^6/mL)')
      dummytab
    }
    
  })
  
  
  output$dummy_text <- renderText({
    if(is.null(values$dv)){
      "Example:"
    }else{
      "Uploaded sample:"
    }
    
  })
  
  
 
  
  
  output$RUI_plot_prob_neut <- renderUI({
    if(!is.null(values$pred_typical) & is.null(values$pred_indiv)){
      if(ncol(values$pred_typical) >= 10){
        fluidRow(style='padding:20px;',
                 actionButton("neut_refresh", strong("Refresh!")),
                 
                 HTML('<p style="font-size:100%;color:blue"> Typical trajectory </p>'),
                 plotOutput("plot_neut_typical_IIV"),
                 tableOutput("table_neut_typical_IIV")
        )
      }
    }else if(!is.null(values$pred_typical) & !is.null(values$pred_indiv)){
      if(ncol(values$pred_typical) >= 10 & ncol(values$pred_indiv) < 10){
        fluidRow(style='padding:20px;',
                 actionButton("neut_refresh", strong("Refresh!")),
                 br(),
                 HTML('<p style="font-size:100%;color:blue"> Typical trajectory </p>'),
                 plotOutput("plot_neut_typical_IIV"),
                 tableOutput("table_neut_typical_IIV")
        )
      }else if(ncol(values$pred_typical) < 10 & ncol(values$pred_indiv) >= 10){
        fluidRow(style='padding:20px;',
                 actionButton("neut_refresh", strong("Refresh!")),
                 br(),
                 HTML('<p style="font-size:100%;color:blue"> Individual trajectory </p>'),
                 plotOutput("plot_neut_indiv_posterior"),
                 tableOutput("table_neut_indiv_posterior")
        )
      }else if(ncol(values$pred_typical) >= 10 & ncol(values$pred_indiv) >= 10){
        fluidRow(style='padding:20px;',
                 actionButton("neut_refresh", strong("Refresh!")),
                 br(),
                 HTML('<p style="font-size:100%;color:blue"> Typical trajectory </p>'),
                 plotOutput("plot_neut_typical_IIV"),
                 tableOutput("table_neut_typical_IIV"),
                 HTML('<p style="font-size:100%;color:blue"> Individual trajectory </p>'),
                 plotOutput("plot_neut_indiv_posterior"),
                 tableOutput("table_neut_indiv_posterior")
        )
      }
      
    }else if(is.null(values$pred_typical) & !is.null(values$pred_indiv)){
      if(ncol(values$pred_indiv) >= 10){
        fluidRow(style='padding:20px;',
                 actionButton("neut_refresh", strong("Refresh!")),
                 br(),
                 HTML('<p style="font-size:100%;color:blue"> Individual trajectory </p>'),
                 plotOutput("plot_neut_indiv_posterior"),
                 tableOutput("table_neut_indiv_posterior"))
      }
    }
    
    

  })
    
  output$plot_neut_typical_IIV <- renderPlot({
    input$neut_refresh

        isolate({
        dosetime <- values$totaldose %>% filter(variable == "dosep") %>% select(time)
        cut.nadir <- cut(values$pred_typical$time, c(dosetime$time, dosetime$time[length(dosetime$time)]+30, input$TIME),include.lowest = T)
        
        values$pred_typical %>% mutate(cut.nadir = cut.nadir) -> temp.pred
        temp.pred <- temp.pred %>% select(-c(1:7)) 
        temp.grade4 <- temp.pred %>% group_by(cut.nadir) %>% summarise(across(1:(ncol(temp.pred)-1), function(x) any(x < 0.5)))
        temp.grade4 <- cbind(temp.grade4[,1], prob = rowSums(temp.grade4[,-1])/ncol(temp.grade4[,-1]))
        
        temp.grade3 <- temp.pred %>% group_by(cut.nadir) %>% summarise(across(1:(ncol(temp.pred)-1), function(x) any(x < 1 & x >= 0.5)))
        temp.grade3 <- cbind(temp.grade3[,1], prob = rowSums(temp.grade3[,-1])/ncol(temp.grade3[,-1]))
        
        
        df <- data.frame(time = temp.grade3$cut.nadir, Grade_3 = temp.grade3$prob, Grade_4 = temp.grade4$prob)
        df.melt <- df %>% melt(value.name = "Probability", id.vars = "time")
        ggplot(df.melt) + geom_bar(aes(x = time, y = Probability, fill = variable ), stat="identity", position = position_dodge()) + ggtitle("Probability of neutropenia based on prior")+
          theme_bw()+
          theme(legend.position = c(0.75, 0.85)) 
        })
  
  })
  
  output$table_neut_typical_IIV <- renderTable({
    input$neut_refresh
    
    dosetime <- values$totaldose %>% filter(variable == "dosep") %>% select(time)
    cut.nadir <- cut(values$pred_typical$time, c(dosetime$time, dosetime$time[length(dosetime$time)]+30, input$TIME),include.lowest = T)
    
    values$pred_typical %>% mutate(cut.nadir = cut.nadir) -> temp.pred
    temp.pred <- temp.pred %>% select(-c(1:7)) 
    temp.grade4 <- temp.pred %>% group_by(cut.nadir) %>% summarise(across(1:(ncol(temp.pred)-1), function(x) any(x < 0.5)))
    temp.grade4 <- cbind(temp.grade4[,1], prob = rowSums(temp.grade4[,-1])/ncol(temp.grade4[,-1]))
    
    temp.grade3 <- temp.pred %>% group_by(cut.nadir) %>% summarise(across(1:(ncol(temp.pred)-1), function(x) any(x < 1 & x >= 0.5)))
    temp.grade3 <- cbind(temp.grade3[,1], prob = rowSums(temp.grade3[,-1])/ncol(temp.grade3[,-1]))
    
    
    df <- data.frame(time = temp.grade3$cut.nadir, Grade_3 = signif(temp.grade3$prob,3), Grade_4 = signif(temp.grade4$prob,3))
    colnames(df) <- c("Time interval (Day)", "Grade 3", "Grade 4")
    time_int <- c(dosetime$time, dosetime$time[length(dosetime$time)]+30, input$TIME)
    time_int <- paste0(time_int[1:(length(time_int)-1)], " ~ ", time_int[2:length(time_int)])
    df[,1] <- time_int
    df
  })
    
  
      
  
  
  output$plot_neut_indiv_posterior <- renderPlot({
    input$neut_refresh
    
    isolate({
          dosetime <- values$totaldose %>% filter(variable == "dosep") %>% select(time)
          cut.nadir <- cut(values$pred_indiv$time, c(dosetime$time, dosetime$time[length(dosetime$time)]+30, input$TIME),include.lowest = T)
          
          values$pred_indiv %>% mutate(cut.nadir = cut.nadir) -> temp.pred
          temp.pred <- temp.pred %>% select(-c(1:7)) 
          temp.grade4 <- temp.pred %>% group_by(cut.nadir) %>% summarise(across(1:(ncol(temp.pred)-1), function(x) any(x < 0.5)))
          temp.grade4 <- cbind(temp.grade4[,1], prob = rowSums(temp.grade4[,-1])/ncol(temp.grade4[,-1]))
          
          temp.grade3 <- temp.pred %>% group_by(cut.nadir) %>% summarise(across(1:(ncol(temp.pred)-1), function(x) any(x < 1 & x >= 0.5)))
          temp.grade3 <- cbind(temp.grade3[,1], prob = rowSums(temp.grade3[,-1])/ncol(temp.grade3[,-1]))
          
          
          df <- data.frame(time = temp.grade3$cut.nadir, Grade_3 = temp.grade3$prob, Grade_4 = temp.grade4$prob)
          df.melt <- df %>% melt(value.name = "Probability", id.vars = "time")
          ggplot(df.melt) + geom_bar(aes(x = time, y = Probability, fill = variable ), stat="identity", position = position_dodge()) + ggtitle("Probability of neutropenia based on posterior")+
            theme_bw()+
            theme(legend.position = c(0.75, 0.85))
    })
  })
  
  
  output$table_neut_indiv_posterior <- renderTable({
    input$neut_refresh
    
    dosetime <- values$totaldose %>% filter(variable == "dosep") %>% select(time)
    cut.nadir <- cut(values$pred_indiv$time, c(dosetime$time, dosetime$time[length(dosetime$time)]+30, input$TIME),include.lowest = T)
    
    values$pred_indiv %>% mutate(cut.nadir = cut.nadir) -> temp.pred
    temp.pred <- temp.pred %>% select(-c(1:7)) 
    temp.grade4 <- temp.pred %>% group_by(cut.nadir) %>% summarise(across(1:(ncol(temp.pred)-1), function(x) any(x < 0.5)))
    temp.grade4 <- cbind(temp.grade4[,1], prob = rowSums(temp.grade4[,-1])/ncol(temp.grade4[,-1]))
    
    temp.grade3 <- temp.pred %>% group_by(cut.nadir) %>% summarise(across(1:(ncol(temp.pred)-1), function(x) any(x < 1 & x >= 0.5)))
    temp.grade3 <- cbind(temp.grade3[,1], prob = rowSums(temp.grade3[,-1])/ncol(temp.grade3[,-1]))
    
    
    df <- data.frame(time = temp.grade3$cut.nadir, Grade_3 = signif(temp.grade3$prob,3), Grade_4 = signif(temp.grade4$prob,3))
    colnames(df) <- c("Time interval (Day)", "Grade 3", "Grade 4")
    time_int <- c(dosetime$time, dosetime$time[length(dosetime$time)]+30, input$TIME)
    time_int <- paste0(time_int[1:(length(time_int)-1)], " ~ ", time_int[2:length(time_int)])
    df[,1] <- time_int
    df
  })
  
  
  output$pred_spec_demo <- renderTable({
    if(is.null(values$bsa)){
      df<- data.frame(BSA = 1.7, ANC = input$BASE,Sex = if_else(input$sex == "1", "Female", "Male"), DM = if_else(input$dm == "1", "Yes", "No")  )
      
    }else(
      df<- data.frame(BSA = values$bsa, ANC = input$BASE,Sex = if_else(input$sex == "1", "Female", "Male"), DM = if_else(input$dm == "1", "Yes", "No")  )
      
    )
    colnames(df)[2] <- "ANC"
    df
    
    
  })
  
  output$pred_spec_dosing <- renderTable({
    
    if(!is.null(values$totaldose)){
      df <- values$totaldose  
      colnames(df) <- c("Time (Day)", "Drug", "Amount (mg)")
      df <- df %>%  mutate(Drug = if_else(Drug == "dosec", "Cisplatin", if_else(Drug == "dosep", "Paclitaxel",if_else(Drug == "doseg", "GCSF", "No") )))  
      df  
    }
    
    
    
  })
  
  output$estim_spec_demo <- renderTable({
    
    if(is.null(values$bsa)){
      df<- data.frame(BSA = 1.7, ANC = input$BASE,Sex = if_else(input$sex == "1", "Female", "Male"), DM = if_else(input$dm == "1", "Yes", "No")  )
      
    }else(
      df<- data.frame(BSA = values$bsa, ANC = input$BASE,Sex = if_else(input$sex == "1", "Female", "Male"), DM = if_else(input$dm == "1", "Yes", "No")  )
      
    )
    colnames(df)[2] <- "ANC"
    df
    
    
  })
  
  output$estim_spec_dosing <- renderTable({
    
    if(!is.null(values$totaldose)){
      df <- values$totaldose  
      colnames(df) <- c("Time (Day)", "Drug", "Amount (mg)")
      df <- df %>%  mutate(Drug = if_else(Drug == "dosec", "Cisplatin", if_else(Drug == "dosep", "Paclitaxel",if_else(Drug == "doseg", "GCSF", "No") )))  
      df  
    }
    
    
    
  })
  
  
  output$dosing_summary <- renderTable({
    
    if(!is.null(values$totaldose)){
      df <- values$totaldose  
      colnames(df) <- c("Time (Day)", "Drug", "Amount (mg)")
      df <- df %>%  mutate(Drug = if_else(Drug == "dosec", "Cisplatin", if_else(Drug == "dosep", "Paclitaxel",if_else(Drug == "doseg", "GCSF", "No") )))  
      df  
    }
    
    
    
  })
  
})



