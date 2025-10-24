library(shiny)
library(shinydashboard)
library(readxl)
library(writexl)
library(DT)
library(dplyr)
library(stringr)

# ---------- Theoretical Mass Calculation Function ----------
add_fragment_mass <- function(df) {
  # Base monoisotopic masses (RNA)
  base_mass <- c(
    A = 329.05252,
    C = 305.04129,
    G = 345.04743,
    U = 306.02530
  )
  
  # Get the fragment sequence column name
  frag_col <- grep("fragment.*sequence", names(df), ignore.case = TRUE, value = TRUE)
  if (length(frag_col) == 0) stop("No column matching 'Fragment sequence' found.")
  
  seqs <- as.character(df[[frag_col]])
  seqs <- trimws(seqs)
  
  compute_mass <- function(seq) {
    if (is.na(seq) || nchar(seq) == 0) return(NA_real_)
    seq <- toupper(seq)
    bases <- unlist(strsplit(seq, ""))
    
    if (!all(bases %in% names(base_mass))) {
      cat("Non-ACGU base found in:", seq, "\n")
      return(NA_real_)
    }
    
    sum(base_mass[bases])
  }
  
  # Compute masses
  df$mass <- sapply(seqs, compute_mass)
  
  # # Add chemical adjustments
  # df$mass <- df$mass + 97.9769 
  # if (nrow(df) > 1) {
  #   df$mass[-nrow(df)] <- df$mass[-nrow(df)] - 18.01
  # }
  # 
  return(df)
}

# ---------- Cleave Simulation Fuction ----------
cleavage <- function(seq_ref, enzyme = "MazF") {
  stopifnot(is.character(seq_ref), length(seq_ref) == 1)
  n <- nchar(seq_ref)
  
  if (n == 0) {
    return(data.frame(
      Cleavage_site = integer(0),
      Cleave_position = integer(0),
      Fragment_size = integer(0),
      Fragment_sequence = character(0)
    ))
  }
  
  # Define recognition motif and where to cut relative to motif
  enzyme_rules <- list(
    MazF = list(motif = "ACA", cut_before = TRUE, offset = 0),
    RNaseT1 = list(motif = "G", cut_before = FALSE, offset = 0),
    RNaseA = list(motif = c("U", "C"), cut_before = FALSE, offset = 0)
  )
  
  if (!enzyme %in% names(enzyme_rules)) {
    stop(paste("Unsupported enzyme:", enzyme))
  }
  
  rule <- enzyme_rules[[enzyme]]
  
  # Find cleavage motif positions
  if (length(rule$motif) == 1) {
    pos <- gregexpr(rule$motif, seq_ref, fixed = TRUE)[[1]]
  } else {
    # multiple motifs
    pos_list <- lapply(rule$motif, function(m) gregexpr(m, seq_ref, fixed = TRUE)[[1]])
    pos <- sort(unlist(pos_list))
  }
  
  # No cleavage motif found → single full-length fragment
  if (identical(pos, -1L)) {
    return(data.frame(
      Cleavage_site = 0L,
      Cleave_position = 0L,
      Fragment_size = n,
      Fragment_sequence = substr(seq_ref, 1, n),
      stringsAsFactors = FALSE
    ))
  }
  
  # Determine cleavage boundary
  cleave_positions <- if (rule$cut_before) {
    pos + rule$offset
  } else {
    pos + nchar(if (length(rule$motif)==1) rule$motif else substr(seq_ref, pos, pos)) + rule$offset
  }
  
  # Build fragments
  starts <- c(1L, cleave_positions)
  ends <- c(cleave_positions - 1L, n)
  keep <- (ends - starts + 1L) > 0L
  starts <- starts[keep]
  ends <- ends[keep]
  cleave_pos_for_frag <- c(0L, cleave_positions)[keep]
  
  df <- data.frame(
    Cleavage_site = seq_along(starts) - 1L,
    Cleave_position = replace(cleave_pos_for_frag, 1, cleave_pos_for_frag[2]),
    Fragment_size = ends - starts + 1L,
    Fragment_sequence = substring(seq_ref, starts, ends),
    stringsAsFactors = FALSE
  )
  
  df <- add_fragment_mass(df)
  names(df)[names(df) == "mass"] <- "Theoretic_mass"
  return(df)
}

# ---------- UI Design ----------
header <- dashboardHeader(
  title = span("RNA Cleavage Simulation Tool", style = "font-weight:600; letter-spacing:0.3px;")
)

sidebar <- dashboardSidebar(
  width = 230,
  sidebarMenu(
    id = "tabs",
    menuItem("Description", tabName = "desc", icon = icon("info-circle")),
    menuItem("Input", tabName = "input", icon = icon("sliders")),
    menuItem("Results", tabName = "results", icon = icon("table"))
  )
)

body <- dashboardBody(
  tags$head(tags$style(HTML("
    body, .content-wrapper, .right-side, .main-sidebar { font-size: 16px !important; line-height: 1.6 !important; }
    .box-title { font-size: 18px !important; font-weight: 600; }
    .content-wrapper, .right-side { background:#f7f9fc; }
    .box { border-radius:12px; overflow: hidden; }
    .skin-blue .main-header .logo { background:#3c8dbc !important; }
    .skin-blue .main-header .navbar { background:#3c8dbc !important; }
    .skin-blue .sidebar-menu>li.active>a { border-left-color:#3c8dbc; }
    .sidebar-mini.sidebar-collapse .main-sidebar { overflow: visible !important; }
  "))),
  
  tabItems(
    # ----- Description -----
    tabItem(tabName = "desc",
            fluidRow(
              style = "display: flex; justify-content: center; padding-top: 20px;",
              box(title = "Overview", width = 8, solidHeader = TRUE, status = "primary",
                  style = "font-size:16px; line-height:1.6;",
                  p(HTML("<b>Brief Description.</b> This application simulates RNA cleavage by selected endoribonucleases and computes the theoretical monoisotopic masses of the resulting fragments. 
                    Given a user-supplied RNA sequence (A, C, G, U only), the algorithm identifies enzyme-specific recognition sites, 
                    generates contiguous fragments, and reports fragment size, sequence, and theoretical mass.")),
                  tags$ul(
                    tags$li("Supported enzymes: MazF (cuts before 'ACA'), RNase T1 (after G), RNase A (after U or C)."),
                    tags$li(HTML("<b>Cleavage_site</b> indicates the index of the fragment, corresponding to the fragment generated after the n-th cleavage site.")),
                    tags$li(HTML("<b>Cleave_position</b> specifies the nucleotide position in the RNA sequence where the cleavage event for that fragment occurs.")),
                    tags$li("Navigate to ", tags$b("Input"), " to enter a sequence and enzyme; view outputs in ", tags$b("Results"), ".")
                  ),
                  p(HTML(
                    "This tool is built with Shiny, and the source code is available on ",
                    "<a href='https://github.com/Shangsi-Lin/mRNA_sequencing' target='_blank'><b>GitHub</b></a>.",
                    "The Shiny application specifically can be found under <code>Source/cleavage_app</code> in the repository."
                  ))
              ),
            )
    ),
    
    # ----- Input -----
    tabItem(tabName = "input",
            fluidRow(
              box(title = "Sequence & Enzyme", width = 8, solidHeader = TRUE, status = "primary",
                  textAreaInput("seq_input", "Reference Sequence (5')",
                                placeholder = "Paste an RNA sequence (A/C/G/U only)", height = "160px"),
                  selectInput("enzyme", "Enzyme", choices = c("MazF", "RNaseT1", "RNaseA"), selected = "MazF", selectize = FALSE),
                  actionButton("run_btn", "Run Simulation", class = "btn btn-primary")
              ),
              box(title = "Notes", width = 4, solidHeader = TRUE, status = "info",
                  p("• Large inputs may take longer to process in the browser."),
                  p("• Only canonical RNA bases (A/C/G/U) are accepted.")
              )
            )
    ),
    
    # ----- Output -----
    tabItem(tabName = "results",
            fluidRow(
              box(title = "Cleavage Results", width = 12, solidHeader = TRUE, status = "primary",
                  div(style = "margin-bottom:10px;",
                      downloadButton("download_table", "Download .xlsx")
                  ),
                  DTOutput("result_table")
              )
            )
    )
  )
)

ui <- dashboardPage(header, sidebar, body, skin = "blue")

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  # Start on Description
  observeEvent(TRUE, { updateTabItems(session, "tabs", "desc") }, once = TRUE)
  
  results <- eventReactive(input$run_btn, {
    req(input$seq_input)
    seq <- gsub("\\s+", "", toupper(input$seq_input))
    validate(need(nchar(seq) > 0, "Please paste a non-empty RNA sequence."))
    cleavage(seq, enzyme = input$enzyme)
  })
  
  output$result_table <- renderDT({
    req(results())
    datatable(results(),
              options = list(pageLength = 10, scrollX = TRUE),
              rownames = FALSE)
  })
  
  output$download_table <- downloadHandler(
    filename = function() paste0("cleavage_results_", input$enzyme, ".xlsx"),
    content  = function(file) write_xlsx(results(), file)
  )
}

shinyApp(ui = ui, server = server)