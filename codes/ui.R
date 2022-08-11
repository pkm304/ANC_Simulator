library(shiny)

shinyUI( #L1
  fluidPage( #L2
    title = "ANC simulator in NSCLC Chemotherapy",
    fluidRow(style='padding:20px;', #L3
      column(12, #L4
            
             img(src="banner.png", style = "max-width: 1000px; width: 100%; height: auto"),
             HTML('<b style="padding:5px;font-size:100%;color:dark blue"> <a href=“https://github.com/pkm304/ANC_Simulator/"> https://github.com/pkm304/ANC_Simulator </a> </b>'),
            
      ) #L4 closed
    ), #L3 closed
    
    fluidRow(
      
      column(8,
             HTML('<br>'),
             HTML('<b style="font-size:120%;color:blue"> Plot </b><br><br>'),
             plotOutput('simPlot'),
             
             br(),
             HTML('<b style="font-size:100%;color:blue"> Abbreviations </b>'),
             HTML('<p style="font-size:100%;color:black"> indiv., individual; traj., trajectory; IIV, interindividual variabilities; POV, posterior variabilities; MAP, maximum a posteriori; ANC, absolute neutrophil count; BSA, body surface area </p>'),
      ),
      column(4, 
             HTML('<br>'),
             HTML('<b style="font-size:120%;color:blue"> Graphical setting </b><br><br>'),
             
             numericInput(inputId = "TIME",
                          label="Simulation time (Days)",
                          min=30, max=2000, value=120),
             numericInput(inputId = "maxy",
                          label="Y-axis range (ANC: 10^6 cells/mL)",
                          min=0.1, max=50, value=10))
    ),
    
    fluidRow(style='padding:20px;',
      br(),
      tabsetPanel(
        tabPanel( HTML('<b style="font-size:120%;color:blue"> Predictions </b>'),style='padding:20px;',
                  fluidRow(
                    br(),
                    column(8,
                           HTML('<b style="font-size:100%;color:blue"> Predict and plot  trajectories </b>'),
                           HTML('<p style="font-size:100%;color:black"> Click the button below to predict and plot trajectories based on the specified baseline demographics and dosing and/or individual parameter estimations.</p>'),
                           
                           fluidRow(
                             
                             
                             column(4,
                                    
                                    actionButton("gen_plot",  HTML('<b> Execute! </b>'))
                                    
                             )
                           ),
                           
                           fluidRow(
                             column(4,
                                    checkboxInput("plot_typical",label = HTML('<b> Typical trajectory </b>'),value = TRUE),
                                    radioButtons(inputId = "pop_var",
                                                 label="With interindividual variabilities",
                                                 choices = list("NO"=0,
                                                                "Yes (may take long..)"=1)),
                                    
                             ),
                             
                             
                             column(4,
                                    uiOutput("RUI_plot_individual"),
                                    uiOutput("RUI_plot_indiv_var")
                                    
                                    
                             )
                           ),
                           br(),
                           HTML('<b style="font-size:100%;color:blue"> Probability of neutropenia  </b>'),
                           HTML('<p style="font-size:100%;color:black"> Probabilites of developing Grade 3 and 4 neutropenia are automatically estimated and depicted below once variabilities in typical and/or individual trajectories are predicted.</p>'),
                           
                           
                           uiOutput("RUI_plot_prob_neut")
                             
                             
                  
                           
                          
                             
                             
                    ),
                    column(4,
                           HTML('<b style="font-size:100%;color:blue"> Specified baseline demographics </b>'),
                           tableOutput("pred_spec_demo"),
                           br(),
                           HTML('<b style="font-size:100%;color:blue"> Specified dosing </b>'),
                           tableOutput("pred_spec_dosing"),
                           br(),
                           HTML('<b style="font-size:100%;color:blue"> Individualized estimation status </b><br><br>'),
                           textOutput("estim_status_1"),
                           uiOutput("RUI_estim_smpl_dose_1"),
                           br(),
                           HTML('<p style="font-size:100%;color:black"> Predictions will be made based on specifications above. You can adjust them in "Patient demographics", "Dosing", and "Individual estimations" tabs.</p>'),
                           
                           
                    )
                    
                    
                    
                    
                    
                  ),
                  
                  
                  hr()
        ),#tabPanel close
        tabPanel(HTML('<b style="font-size:120%;color:blue"> Patient demographics </b>'),style='padding:20px;',
                 fluidRow(
                 
                   br(),
                   column(4,#L7
                          uiOutput("RUI_bsa")
                   ),
                          
                   column(4, #L7
                          sliderInput("BASE", "Baseline ANC (10^6 / mL)",
                                      min = 0.1, max = 20,
                                      value = 5.324, step = 0.01)
                   )
                 ),
                 fluidRow(
                   column(4,#L7
                          radioButtons(inputId = "sex",
                                       label="Sex",
                                       choices = list("Male"=0,
                                                      "Female"=1),
                                       selected = 0)), #L7 closed
                   column(4, #L7
                          radioButtons(inputId = "dm",
                                       label = "DM",
                                       choices = list("No" = 0,
                                                      "Yes" = 1),
                                       selected = 0)
                   )
                   
                 ),
                 hr()
        ),#tabPanel close
        tabPanel(
          HTML('<b style="font-size:120%;color:blue"> Dosing </b>'),style='padding:20px;',
          fluidRow( 
            br(),
            column(5,
                  
                   HTML('<b style="font-size:100%;color:blue"> Pacltaxel/Cisplatin </b>'),
                   div(id="placeholderchem"),
                   actionButton("add_chem_dose", strong("Add Pacltaxel/Cisplatin Doses"))
                   
            ),
            column(3,
                   HTML('<b style="font-size:100%;color:blue"> GCSF </b>'),
                   div(id="placeholdergcsf"),
                   actionButton("add_gcsf_dose", strong("Add GCSF Doses"))
            ),
            column(4,
                   HTML('<b style="font-size:100%;color:blue"> Dosing summary </b>'),
                   tableOutput("dosing_summary")
            )
            
          ),
          hr()
          
         
          
          
        ),
        tabPanel(
          HTML('<b style="font-size:120%;color:blue"> Individual estimations </b>'),style='padding:20px;',
          fluidRow(
            br(),
            column(4,
                   HTML('<b style="font-size:100%;color:blue"> Individual sample data file upload </b><br><br>'),
                   
                   HTML('<p style="color:black"> The data file should be comma-separated (i.e. csv files) and have the following format arranged longitudinally :</p>'),
                   textOutput('dummy_text'),
                   tableOutput('dummy'),
                   br(),
                   
                   fileInput('dvdata', 'Choose CSV File', multiple = FALSE,
                             accept = c('text/csv','text/comma-separated-values,text/plain','.csv')),
                   HTML('<br>'),
                   HTML('<b style="font-size:100%;color:blue">  Click the button below to reset data.</b>'),
                   HTML('<br>'),
                   actionButton("reset", "Reset")
                   ),
            column(4,
                   HTML('<b style="font-size:100%;color:blue"> Specified baseline demographics </b>'),
                   tableOutput("estim_spec_demo"),
                   br(),
                   HTML('<b style="font-size:100%;color:blue"> Specified dosing </b>'),
                   tableOutput("estim_spec_dosing"),
                   br(),
                   HTML('<b style="font-size:100%;color:blue"> Estimate Individual parameters </b><br><br>'),
                   HTML('<p style="color:black"> Click the button below once appears. Individualized parameters (and their posterior distribution) will be estimated based on the specified baseline demographics and dosing and the uploaded sample data. </p>'),
                   
                   uiOutput("RUI_estim_prm_individual"),
                   uiOutput("RUI_estim_indiv_var")
                   
                     
                   ),
            column(4,
                   HTML('<b style="font-size:100%;color:blue"> Estimation status </b><br><br>'),
                   textOutput("estim_status"),
                   uiOutput("RUI_estim_smpl_dose"),
            )
                   
            
          
           
          ),
          hr()
          
        )#,
        # tabPanel( HTML('<b style="font-size:120%;color:blue"> Tutorial </b>'),style='padding:20px;',
        #          hr()
        # ),
        # tabPanel(HTML('<b style="font-size:120%;color:blue"> About </b>'),style='padding:20px;',
        #          hr()
        # )
        
      )# tabsetPanel close
    )
    
  )
  
  # tabPanel('Dosing simulation', #L6             
  #          column(4, #L7
  #                 HTML('<br><br>'),
  #                 HTML('<b style="font-size:120%;color:blue"> Dosing </b><br><br>'),  
  #          
  #                 numericInput(inputId = "DOSE",
  #                       label="Dosing Amount [IU/dL]",
  #                       min=2, max=10000, value=1000),
  #          
  #                 numericInput(inputId = "II",
  #                       label="Dosing interval [day]",
  #                       min=0.5, max=30, value=3),
  #                 
  #                 numericInput(inputId = "ADDL",
  #                              label="Additional doses",
  #                              min=0, max=9999, value=0)
  #                 ) #L7 closedTYPE
  #          ), #L6 closed
  # tabPanel('Graphics', #L6 
  #          column(4, #L7
  #                 HTML('<br><br>'),
  #                 HTML('<b style="font-size:120%;color:blue"> Graphical Display </b><br><br>'), 
  #          
  #                 numericInput(inputId = "TIME",
  #                       label="simulation time",
  #                       min=24, max=1200, value=120),
  #                 numericInput(inputId = "maxy",
  #                       label="Y-axis range",
  #                       min=10, max=300, value=100)
  #                 ) #L7 closed
  #        ) #L6 closed
  #) #L5 closed
  #) #L4 closed
  #) #L3 closed
) #L2 closed
#L1 closed


