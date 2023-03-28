# Run in base_dir
# ml R/4.0

library(rsconnect)
library(pkgconfig)

typeouts <- c()

rsconnect::setAccountInfo(
  name = 'longo-stanford', 
  token = 'XXXXXXXXXX', 
  secret = 'XXXXXXXXXX'
)

deployApp(
  appFiles = c(
    'app.R',
    'www/style.css',
    list.files(typeouts, pattern = '\\.rdata|Modules_Up-Down\\.full\\.logfdr\\.png|TREAT-AD\\.correlations\\.png', full.names = T)
  ),
  appName = 'APP_NAME',
  server = 'shinyapps.io'
)
