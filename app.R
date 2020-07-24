library(shinythemes); library(shinyBS)
library(shinycssloaders)
library(tidyverse)

source("immune_annotator_functions.R")

immune_key <- make_immune_key()
subset_key <- make_immune_subset_key()
name_order_key <- make_immune_name_order()


ui <- navbarPage("Immune Subset Annotator", theme = shinytheme("flatly"), collapsible = T,
                 header = 
                   tags$head(
                     includeHTML("google-analytics.js"),
                     tags$style(HTML("
                        #test {
                          padding: 100px;
                        }
                        .navbar {
                          margin: 0px;
                        }
                        .footer {
                            position: relative;
                            left: 0;
                            bottom: 0;
                            width: 100%;
                            background-color: #d7dfea;
                            # color: white;
                            text-align: center;
                        }
                        "))
                   ),
                 
                 tabPanel("Annotator", id="test", 
                          sidebarLayout(
                            sidebarPanel(width = 2,
                                         textAreaInput(inputId="caption", 
                                                       label=HTML("<h3>Put text here</h5><p>1 cell per line</p>"), 
                                                       width="100%", height="400px", value="", placeholder = "1000 character limit"),
                                         actionButton("go", label = "Submit")
                            ),
                            mainPanel(
                              h3("Annotation output"),
                              withSpinner(DT::dataTableOutput("mytable")),
                              downloadButton('downloadData', 'Download annotated data')
                            )
                          )
                 )
)

server <- function(input, output) {
  data.browse <- eventReactive(input$go,{
    values <- substr(input$caption, 1,1000)

    df <- tibble(value = unlist(str_split(values, "\n"))) 
    
    df %>% filter(grepl(paste(immune_key$pattern, collapse = "|"), value))  %>% 
      mutate(lineage = map_chr(value, function(x) pattern_match(x, immune_key))) %>% 
      select(value, lineage) %>% 
      left_join(x=df,y=., by = "value") %>% 
      mutate(lineage = str_remove(lineage, "_monocyte|monocyte_")) %>% 
      mutate(subset = 
               map2_chr(value, lineage, 
                        function(x, y) pattern_match_3(x, y, subset_key, name_order_key))) 
    
  })
  output$mytable = DT::renderDataTable({data.browse()}, rownames = F)
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("immune_annotations_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(data.browse(), file,row.names = FALSE)
    }
  )

}

shinyApp(ui, server)
