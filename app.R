library(shiny)
library(shinyMobile)
library(shinyjs)
library(future)
library(promises)
library(dplyr)
library(tidyr)
library(purrr)
library(Hmisc) # for capitalize
library(ggplot2)
library(readr)
library(highcharter)
library(zoo)
library(httr)
library(jsonlite)
library(leaflet)
library(leaflet.providers)
library(leaflet.extras)

plan(multisession)

# Modelos de Carvalho & Del Ponte (2026), Plant Pathology 75:e70173

# Helper function for date handling
as_date_safe <- function(x) {
  if (is.null(x)) return(Sys.Date())
  d <- tryCatch(as.Date(x), error = function(e) NA)
  if (all(is.na(d))) return(Sys.Date())
  d
}

# --- DATA FETCHING FUNCTIONS (Open-Meteo) ---

# Fetch Archive Data (Historical)
fetch_archive_data <- function(lat, lon, start_date, end_date) {
  url <- "https://archive-api.open-meteo.com/v1/archive"
  s_str <- as.character(start_date)
  e_str <- as.character(end_date)
  params <- list(
    latitude = lat, longitude = lon,
    start_date = s_str, end_date = e_str,
    hourly = "temperature_2m,relative_humidity_2m,dew_point_2m,precipitation,rain",
    timezone = "auto"
  )
  resp <- tryCatch({ GET(url, query = params, timeout(30)) }, error = function(e) NULL)
  if (is.null(resp) || status_code(resp) != 200) return(NULL)
  dat <- tryCatch(fromJSON(content(resp, "text", encoding = "UTF-8")), error = function(e) NULL)
  if (is.null(dat) || is.null(dat$hourly)) return(NULL)
  df <- data.frame(
    time = as.POSIXct(dat$hourly$time, format="%Y-%m-%dT%H:%M", tz="UTC"),
    T2M = dat$hourly$temperature_2m,
    RH2M = dat$hourly$relative_humidity_2m,
    DEW = dat$hourly$dew_point_2m,
    PREC = dat$hourly$precipitation,
    SOURCE = "archive"
  )
  return(df)
}

# Fetch Forecast Data (Future)
fetch_forecast_data <- function(lat, lon) {
  url <- "https://api.open-meteo.com/v1/forecast"
  params <- list(
    latitude = lat, longitude = lon,
    hourly = "temperature_2m,relative_humidity_2m,dew_point_2m,precipitation,rain",
    past_days = 92, forecast_days = 16,
    timezone = "auto"
  )
  resp <- tryCatch({ GET(url, query = params, timeout(30)) }, error = function(e) NULL)
  if (is.null(resp) || status_code(resp) != 200) return(NULL)
  dat <- tryCatch(fromJSON(content(resp, "text", encoding = "UTF-8")), error = function(e) NULL)
  if (is.null(dat) || is.null(dat$hourly)) return(NULL)
  df <- data.frame(
    time = as.POSIXct(dat$hourly$time, format="%Y-%m-%dT%H:%M", tz="UTC"),
    T2M = dat$hourly$temperature_2m,
    RH2M = dat$hourly$relative_humidity_2m,
    DEW = dat$hourly$dew_point_2m,
    PREC = dat$hourly$precipitation,
    SOURCE = "forecast"
  )
  return(df)
}

# Unified function to get stitched weather series
fetch_all_weather <- function(lat, lon, start_req, end_req) {
  df_fore <- fetch_forecast_data(lat, lon)
  archive_end <- Sys.Date() - 5
  df_arch <- NULL
  if (start_req < archive_end) { df_arch <- fetch_archive_data(lat, lon, start_req, archive_end) }
  df_all <- rbind(df_arch, df_fore) 
  if (is.null(df_all) || nrow(df_all) == 0) return(data.frame())
  df_all <- df_all |> 
    distinct(time, .keep_all = TRUE) |> 
    arrange(time) |>
    filter(as.Date(time) >= (start_req - 2) & as.Date(time) <= (end_req + 2))
  df_daily <- df_all |>
    mutate(DATE = as.Date(time)) |>
    group_by(DATE) |>
    summarise(
      T2M_MIN = min(T2M, na.rm = TRUE),
      T2M_MAX = max(T2M, na.rm = TRUE),
      T2M_MEAN = mean(T2M, na.rm = TRUE),
      RH2M = mean(RH2M, na.rm = TRUE),
      T2M_DEW = mean(DEW, na.rm = TRUE),
      PRECTOTCORR = sum(PREC, na.rm = TRUE),
      N_HOURS = sum(!is.na(T2M)),
      SOURCE = if (any(SOURCE == "forecast")) "forecast" else "archive",
      .groups = "drop"
    ) |>
    # descarta dias incompletos (ex.: borda do fuso horário) para não enviesar min/soma
    filter(N_HOURS >= 20) |>
    mutate(
      YEAR = as.numeric(format(DATE, "%Y")),
      MO   = as.numeric(format(DATE, "%m")),
      DY   = as.numeric(format(DATE, "%d"))
    )
  return(df_daily)
}

# --- RISK MODEL: Carvalho & Del Ponte (2026), coeficientes de fhb_models.rds (repo giberela2) ---
# Implementado sem o pacote rms para não depender da versão do rms no servidor.
# rcs() de Harrell (rms/Hmisc), normalização padrão: divide por (t_k - t_1)^2
rcs_basis <- function(x, knots) {
  k <- length(knots); t1 <- knots[1]; tk <- knots[k]; tk1 <- knots[k - 1]
  pos3 <- function(z) pmax(z, 0)^3
  sapply(seq_len(k - 2), function(j) {
    tj <- knots[j]
    (pos3(x - tj) - pos3(x - tk1) * (tk - tj) / (tk - tk1) + pos3(x - tk) * (tk1 - tj) / (tk - tk1)) / (tk - t1)^2
  }, simplify = "matrix")
}
fhb_predict <- function(tmin, rh, dew, prec2) {
  # nós (knots) dos splines, iguais aos de fhb_models.rds
  KNOTS_RH  <- c(64.64567, 74.33071, 80.35433, 87.40157)
  KNOTS_DEW <- c(10.22943, 14.02429, 17.11000)
  b_rh  <- matrix(rcs_basis(rh, KNOTS_RH), ncol = 2)
  b_dew <- matrix(rcs_basis(dew, KNOTS_DEW), ncol = 1)
  lp1 <- -2.87818740 + 0.53883988 * tmin - 0.07932981 * rh + 0.34703191 * b_rh[, 1] - 1.21812658 * b_rh[, 2]
  lp2 <- 5.740245274 - 0.117068699 * rh + 0.482242121 * b_rh[, 1] - 2.020216240 * b_rh[, 2] -
         0.006828326 * dew + 0.532288414 * b_dew[, 1]
  lp3 <- -8.45717860 + 0.57942204 * tmin + 0.01700918 * prec2
  p1 <- plogis(lp1); p2 <- plogis(lp2); p3 <- plogis(lp3)
  stk <- plogis(-2.8931398 + 0.8888407 * p1 + 2.7411160 * p2 + 2.3117934 * p3)
  data.frame(p1 = p1, p2 = p2, p3 = p3, ensemble_mean = (p1 + p2 + p3) / 3, ensemble_stacked = stk)
}

# --- Preditores e janela de risco ---
# Dia 0 = início da antese (data de florescimento simulada). Preditores:
#   tmin  = média da temperatura mínima diária, dias 2 a 10
#   rh    = média da UR média diária, dias 5 a 10
#   dew   = média do ponto de orvalho médio diário, dias 4 a 10
#   prec2 = chuva acumulada (mm), dias 6 a 10
# Só calcula o risco quando todos os dias 2 a 10 existem; caso contrário devolve NA
# (por exemplo, quando a janela passa do fim da previsão de 16 dias).
calc_predictors <- function(weather_df, d) {
  d <- as.Date(d)
  w <- weather_df[weather_df$DATE >= d + 2 & weather_df$DATE <= d + 10, ]
  if (!all((d + 2:10) %in% w$DATE)) return(NULL)
  sel <- function(a, b) w$DATE >= d + a & w$DATE <= d + b
  data.frame(
    tmin  = mean(w$T2M_MIN[sel(2, 10)]),
    rh    = mean(w$RH2M[sel(5, 10)]),
    dew   = mean(w$T2M_DEW[sel(4, 10)]),
    prec2 = sum(w$PRECTOTCORR[sel(6, 10)]),
    uses_forecast = any(w$SOURCE == "forecast")
  )
}

run_window <- function(weather_df, start_date, end_date) {
  dates <- seq(as.Date(start_date), as.Date(end_date), by = "day")
  rows <- lapply(dates, function(d) {
    x <- calc_predictors(weather_df, d)
    if (is.null(x) || anyNA(x[, c("tmin", "rh", "dew", "prec2")])) {
      return(data.frame(flowering_date = d, tmin = NA, rh = NA, dew = NA, prec2 = NA,
                        uses_forecast = NA, p1 = NA, p2 = NA, p3 = NA,
                        ensemble_mean = NA, ensemble_stacked = NA))
    }
    data.frame(flowering_date = d, x, fhb_predict(x$tmin, x$rh, x$dew, x$prec2))
  })
  dplyr::arrange(do.call(rbind, rows), flowering_date)
}


# ======================= UI =======================

ui <- f7Page(
  title = "Monitor Giberela",
  options = list(theme = "ios", dark = FALSE, color = "#4682B4"),
  
  tags$head(tags$style(HTML('
    body { margin: 0; padding: 0; overflow: hidden; background-color: #f4f7f9; }
    #splash_screen { margin: 0; padding: 0; }
    
    .block-title {
      color: #2c5270 !important;
      font-weight: 600 !important;
      font-size: 16px !important;
      margin-top: 20px !important;
      margin-bottom: 10px !important;
      text-transform: none !important;
    }

    /* Navigation Bar Styling */
    .navbar, .navbar-inner, .toolbar, .toolbar-inner {
      background-color: #2c5270 !important;
      color: #ffffff !important;
    }
    
    .navbar .title, .navbar .left, .navbar .right, .navbar a {
      color: #ffffff !important;
    }
    
    .toolbar a, .tab-link, .tab-link-active {
      color: rgba(255, 255, 255, 0.7) !important;
    }
    
    .tab-link-active {
      color: #ffffff !important;
    }
    
    .tab-link-active .icon, .tab-link-active .tabbar-label {
      color: #ffffff !important;
    }

    .toast {
      background-color: rgba(255, 255, 255, 0.95) !important;
      color: #333 !important;
      box-shadow: 0 4px 12px rgba(0,0,0,0.15) !important;
      border: 1px solid #eee !important;
      border-radius: 8px !important;
    }
    .toast .toast-text {
      color: #333 !important;
      font-weight: 500;
      font-family: -apple-system, system-ui, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    }
  '))),
  
  tags$div(
    id = "splash_screen",
    style = "
      position: fixed; top: 0; left: 0; width: 100%; height: 100%;
      background: url('splash_wheat_blue_sky.png') center center / cover no-repeat;
      z-index: 99999;
      display: flex; flex-direction: column; align-items: center; justify-content: center;
    ",
    tags$div(
      style = "
        position: absolute; top: 0; left: 0; width: 100%; height: 100%;
        display: flex; flex-direction: column; align-items: center; justify-content: center;
        backdrop-filter: blur(1px);
        background: rgba(44, 82, 112, 0.4);
      ",
      tags$div(
        style = "
           text-align: center; 
           color: #ffffff;
           max-width: 90%;
        ",
        tags$h3("Monitor climático: giberela do trigo", 
          style = "font-family: 'Roboto', sans-serif; font-size: 26px; font-weight: 700; color: #ffffff; text-shadow: 2px 2px 4px rgba(0,0,0,0.8); margin-bottom: 5px;"
        ),
        tags$p("Universidade Federal de Viçosa", 
          style = "font-size: 14px; font-weight: 500; text-shadow: 1px 1px 2px rgba(0,0,0,0.8); margin: 0; line-height: 1.2;"
        ),
        tags$p("Rede Fitossanidade Tropical", 
          style = "font-size: 14px; font-weight: 500; text-shadow: 1px 1px 2px rgba(0,0,0,0.8); margin: 0; line-height: 1.2; margin-bottom: 20px;"
        ),
        tags$style(HTML("
          #enter_app {
            background-color: rgba(255,255,255,0.2) !important;
            border: 1px solid #fff !important;
            color: #fff !important;
            text-shadow: 1px 1px 2px rgba(0,0,0,0.5);
            font-weight: 600;
            width: 200px;
            margin: 0 auto;
            backdrop-filter: blur(5px);
          }
        ")),
        f7Button(
          inputId = "enter_app", 
          label = "Entrar no sistema", 
          fill = TRUE, 
          color = "gray", 
          size = "large"
        )
      )
    )
  ),
  
  f7TabLayout(
    useShinyjs(),
    navbar = f7Navbar(title = "", hairline = FALSE, shadow = FALSE, leftPanel = FALSE),
    f7Tabs(
      id = "tabs", animated = TRUE, swipeable = FALSE,
      f7Tab(
        title = "Mapa", tabName = "tab_mapa", icon = f7Icon("map"), active = TRUE,
        tags$style(HTML("
          #map { height: calc(100vh - 100px) !important; width: 100%; border:0; }
          #tab_mapa .page-content { padding: 0 !important; overflow: hidden; }
          .center-pin {
            position: absolute; top: 50%; left: 50%; transform: translate(-50%, -100%);
            z-index: 2000; pointer-events: none; display: flex; flex-direction: column; align-items: center;
          }
          .pin-drop {
            width: 32px; height: 32px; background: #e74c3c; border-radius: 50% 50% 50% 0; transform: rotate(-45deg);
            border: 2px solid white; box-shadow: 0 2px 4px rgba(0,0,0,0.3);
          }
          .pin-drop::after {
            content: ''; width: 10px; height: 10px; background: white; position: absolute;
            border-radius: 50%; top: 50%; left: 50%; transform: translate(-50%, -50%);
          }
          .floating-date-btn { position: absolute; top: 70px; right: 15px; z-index: 9999; display: flex; flex-direction: column; gap: 10px; }
          #run_btn, #fab_config {
            background-color: rgba(70, 130, 180, 0.9) !important; color: white !important; border: 1px solid rgba(70, 130, 180, 1) !important;
            box-shadow: 0 4px 6px rgba(0,0,0,0.3) !important; text-transform: none !important; font-weight: 600 !important;
          }
        ")),
        tags$div(class = "center-pin", tags$div(class = "pin-drop")),
        leafletOutput("map"),
        tags$div(class = "floating-date-btn",
          f7Button(inputId = "run_btn", label = "Simular", icon = f7Icon("play_fill"), fill = TRUE, color = "green", rounded = TRUE, size = "large", shadow = TRUE),
          f7Button(inputId = "fab_config", label = "Data", icon = f7Icon("calendar"), fill = TRUE, color = "pink", rounded = TRUE, size = "large", shadow = TRUE),
          div(style = "display: none;",
            f7DatePicker(inputId = "ref_date", label = "Data", value = Sys.Date(), openIn = "auto", toolbar = TRUE)
          )
        )
      ),

      f7Tab(title = "Gráficos", tabName = "tab_graficos", icon = f7Icon("chart_bar_fill"),
        f7BlockTitle("Risco Estimado"),
        highchartOutput("risk_hc", height = "300px"),
        f7BlockTitle("Temperatura (°C)"),
        highchartOutput("meteo_temp", height = "250px"),
        f7BlockTitle("Umidade e Chuva"),
        highchartOutput("meteo_rh_prec", height = "250px")
      ),
      f7Tab(title = "Sobre", tabName = "tab_sobre", icon = f7Icon("info_circle"),
        f7Block(
          f7Card(
            tags$div(
              style = "font-size: 15px; line-height: 1.6; color: #333; text-align: justify;",
              tags$p("O aplicativo implementa o modelo preditivo de risco proposto por Carvalho & Del Ponte (2026), utilizando o ensemble empilhado (stacked): um meta-modelo logístico que combina as probabilidades dos três modelos logísticos individuais (temperatura mínima, umidade relativa, ponto de orvalho e chuva nos dias seguintes à antese). Valores acima de 0,5 indicam risco de epidemia (severidade > 10%)."),
              tags$p("Desenvolvido pelo Laboratório de Epidemiologia da Universidade Federal de Viçosa (UFV) em parceria com a Rede Fitossanidade Tropical (RFT), o app integra dados meteorológicos históricos e de previsão derivados da reanálise ERA5, disponibilizados via API da Open-Meteo, para estimar o risco de ocorrência da doença de forma contínua no tempo e no espaço."),
              tags$p("A ferramenta foi concebida para apoiar pesquisadores, técnicos e produtores na interpretação do risco climático associado à doença, permitindo consultas rápidas, simulações para datas passadas e projeções baseadas na previsão do tempo."),
              tags$p("Contato: delponte@ufv.br", style = "margin-top: 20px; font-weight: 500;")
            )
          )
        )
      )
    )
  )
)


# ======================= SERVER =======================

server <- function(input, output, session) {
  
  # Splash Screen Removal
  observeEvent(input$enter_app, {
    shinyjs::hide("splash_screen")
  })
  
  # Reactive values exactly as app_exemplo.R
  coords <- reactiveValues(lat = -29.2620, lon = -52.4060, zoom = 6)
  risk_result <- reactiveVal(NULL)
  weather_result <- reactiveVal(NULL)
  
  output$map <- renderLeaflet({
    leaflet(options = leafletOptions(zoomControl = FALSE, attributionControl = FALSE)) |> 
      addProviderTiles(providers$Stadia.Outdoors) %>%
      setView(lng = isolate(coords$lon), lat = isolate(coords$lat), zoom = isolate(coords$zoom)) |>
      addSearchOSM(options = searchOptions(autoCollapse = TRUE, minLength = 2, hideMarkerOnCollapse = TRUE))
  })
  
  # Background monitoring - PROTECTED BY TAB CHECK
  observe({
    # TRAP: Only update if the map is the active tab. 
    # This prevents hidden map resets from corrupting the coordinates.
    req(input$tabs == "tab_mapa")
    
    center <- input$map_center
    zoom <- input$map_zoom
    if (!is.null(center)) {
      coords$lat <- center$lat
      coords$lon <- center$lng
    }
    if (!is.null(zoom)) {
      coords$zoom <- zoom
    }
  })
  
  observeEvent(input$fab_config, { 
    shinyjs::runjs("if(typeof app !== 'undefined') { var cal = app.calendar.get('#ref_date'); if(cal) cal.open(); else $('#ref_date').click(); } else { $('#ref_date').click(); }")
  })
  
  observeEvent(input$run_btn, {
    # SNAPSHOT MAP STATE (FIRST STEP)
    center <- input$map_center
    zoom <- input$map_zoom
    if (!is.null(center)) {
      coords$lat <- center$lat
      coords$lon <- center$lng
    }
    if (!is.null(zoom)) {
      coords$zoom <- zoom
    }

    safe_date <- input$ref_date
    date_is_valid <- FALSE
    if (!is.null(safe_date) && length(safe_date) > 0) {
      if (!any(is.na(safe_date)) && as.character(safe_date) != "") date_is_valid <- TRUE
    }
    if (!date_is_valid) { f7Toast("Erro: Data inválida.", position = "bottom", color = "red"); return() }
    
    # 1. Tab change first (Quietly, no toast as per user request)
    updateF7Tabs(session, id = "tabs", selected = "tab_graficos")
    
    # 2. Processamento
    lat_sel <- coords$lat
    lon_sel <- coords$lon
    d_val   <- tryCatch({ as.Date(input$ref_date) }, error = function(e) Sys.Date())
    vals    <- list(start = d_val - 10, end = d_val + 10)
    
    future({
       library(dplyr); library(tidyr); library(purrr); library(Hmisc); library(zoo); library(httr); library(jsonlite); library(tibble)
       df <- fetch_all_weather(lat_sel, lon_sel, vals$start, vals$end + 10)
       risk <- if (nrow(df) > 0) tryCatch(run_window(df, vals$start, vals$end), error=function(e) NULL) else NULL
       list(w = df, r = risk)
    }, seed = TRUE) %...>% 
      (function(res) {
          weather_result(res$w)
          risk_result(res$r)
          # REMOVED "Cálculo concluído" TOAST TO PREVENT UI REDRAWS
      }) %...!% 
      (function(e) {
          f7Toast("Erro no processamento.", type = "error", position = "bottom")
          print(e)
      })
    NULL
  })
  
  output$risk_hc <- renderHighchart({
    rs <- risk_result()
    shiny::validate(need(!is.null(rs) && nrow(rs) > 0, "Aguardando cálculo..."))
    dates <- as_date_safe(rs$flowering_date)
    rs$x <- datetime_to_timestamp(as.POSIXct(dates, tz = "UTC"))
    rs$y <- round(rs$ensemble_stacked, 3)  # ensemble empilhado (stacked)
    df_obs  <- rs |> dplyr::filter(!is.na(y) & !uses_forecast) |> dplyr::select(x, y)
    df_prev <- rs |> dplyr::filter(!is.na(y) &  uses_forecast) |> dplyr::select(x, y)
    hc <- highchart() |>
      hc_chart(zoomType = "x", backgroundColor = "transparent") |>
      hc_title(text = "") |>
      hc_xAxis(type = "datetime") |>
      hc_yAxis(min = 0, max = 1, gridLineColor = "#e0e0e0",
               plotLines = list(list(value = 0.5, color = "red", width = 2, dashStyle = "shortdash", label = list(text = "Limiar", style = list(color="red", fontWeight="bold"))))) |>
      hc_plotOptions(column = list(pointPadding = 0.05, groupPadding = 0.1, borderWidth = 0, grouping = FALSE)) |>
      hc_tooltip(valueDecimals = 2)
    if (nrow(df_obs) > 0)
      hc <- hc |> hc_add_series(data = list_parse2(df_obs), type = "column", name = "Risco (clima observado)", color = "rgba(70, 130, 180, 0.8)")
    if (nrow(df_prev) > 0)
      hc <- hc |> hc_add_series(data = list_parse2(df_prev), type = "column", name = "Risco (inclui previsão)", color = "rgba(230, 126, 34, 0.7)")
    hc |> hc_add_theme(hc_theme_flat())
  })
  
  output$meteo_temp <- renderHighchart({
    df0 <- weather_result()
    shiny::validate(need(!is.null(df0) && nrow(df0) > 0, "Sem dados."))
    df0 <- df0 |> dplyr::mutate(DATE = as_date_safe(DATE))
    highchart() |> hc_chart(backgroundColor = "transparent") |> hc_xAxis(type = "datetime") |> 
      hc_add_series(data = list_parse2(tibble(x = datetime_to_timestamp(as.POSIXct(df0$DATE, tz="UTC")), y = df0$T2M_MIN)), type = "line", name = "Min") |> 
      hc_add_series(data = list_parse2(tibble(x = datetime_to_timestamp(as.POSIXct(df0$DATE, tz="UTC")), y = df0$T2M_MAX)), type = "line", name = "Max") |> 
      hc_add_theme(hc_theme_flat())
  })
  
  output$meteo_rh_prec <- renderHighchart({
    df0 <- weather_result()
    shiny::validate(need(!is.null(df0) && nrow(df0) > 0, "Sem dados."))
    df0 <- df0 |> dplyr::mutate(DATE = as_date_safe(DATE))
    highchart() |> hc_chart(backgroundColor = "transparent") |> hc_xAxis(type = "datetime") |> 
      hc_yAxis_multiples(list(title = list(text = "UR (%)"), height = "50%"), list(title = list(text = "Chuva (mm)"), top = "55%", height = "45%")) |> 
      hc_add_series(data = list_parse2(tibble(x = datetime_to_timestamp(as.POSIXct(df0$DATE, tz="UTC")), y = df0$RH2M)), type = "line", name = "UR", yAxis = 0) |> 
      hc_add_series(data = list_parse2(tibble(x = datetime_to_timestamp(as.POSIXct(df0$DATE, tz="UTC")), y = df0$PRECTOTCORR)), type = "column", name = "Chuva", yAxis = 1) |> 
      hc_add_theme(hc_theme_flat())
  })
  
  # ROBUST RESTORATION WITH DELAY (Crucial for mobile tabs)
  observeEvent(input$tabs, {
    if (input$tabs == "tab_mapa") {
       # 1. Force recalculation of container size
       shinyjs::runjs("setTimeout(function(){ var map = HTMLWidgets.find('#map'); if(map) map.getMap().invalidateSize(); }, 200);")
       
       # 2. Restore view with a safety delay
       shinyjs::delay(350, {
         leafletProxy("map") %>%
           setView(lng = coords$lon, lat = coords$lat, zoom = coords$zoom)
       })
    }
  }, ignoreInit = TRUE)
}

shinyApp(ui, server)