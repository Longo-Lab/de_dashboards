# Run in base_dir
# ml R/4.2.2

library(rsconnect)
library(pkgconfig)
library(stringr)

typeouts <- list.files(pattern = 'dashboard_files\\.rdata', recursive = T) %>%
  str_extract('[^/]+')

rsconnect::setAccountInfo(
  name = 'longo-stanford', 
  token = 'XXXXXXXXXX', 
  secret = 'XXXXXXXXXX'
)

deployApp(
  appFiles = c(
    'app.R',
    'www/style.css',
    'www/script.js',
    list.files(typeouts, pattern = 'dashboard_files\\.rdata|Modules_Up-Down\\.full\\.logfdr\\.png|TREAT-AD\\.correlations\\.png', full.names = T)
  ),
  appName = 'APP_NAME',
  server = 'shinyapps.io'
)
