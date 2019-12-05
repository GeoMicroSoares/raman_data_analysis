library(tidyverse)
library(ggplot2)

# Create a tibble to store each spectra
spec <- tibble()
# and a data.frame to store the fit parameters
fitpar <- data.frame()

# # Find all the files in "Data"
# norm_spls <- list.files(path="Data_means", pattern = "PBH481_10mM")

norm_spls = c("PBS1_05mM_r1_means.txt", "PBS1_1mM_r1_means.txt", "PBS1_1mM_r2_means.txt", 
              "PBS1_5mM_r1_means.txt", "PBS1_10mM_r1_means.txt","PBH481_05mM_r1_means.txt", 
              "PBH481_1mM_r1_means.txt", "PBH481_1mM_r2_means.txt", "PBH481_1mM_r3_means.txt",
              "PBH481_5mM_r1_means.txt", "PBH481_10mM_r1_means.txt")

for(i in seq_along(norm_spls)){#i <- 1
  d <- read_tsv(file.path("Data_means",norm_spls[i]), 
                col_names = c("w", "Int") )
  d = d%>%filter(w>205 & w<290)
  bg = baseline(d$w, d$Int)$bg
  xmax  <- d$w[which.max(d$Int)] 
  Guess <- c(y0 = 0.01, # constant background
             x = c(230,235,250), # positions
             FWHM = c(10, 10, 10), # full width at half maximum
             A = rep(max(d$Int)*2, 3))
  fit <- nls(data =d, 
             Int-bg ~ y0 + 
               A1*Lor(w,x1,FWHM1) + 
               A2*Lor(w,x2,FWHM2) + 
               A3*Lor(w,x3,FWHM3), 
             start=as.list(Guess), 
             lower=as.list(Guess*0),
             algorithm = "port")
  p  <- coef(fit)
  y0 <- p['y0']
  y1 <- y0 + p['A1']*Lor(d$w, x0=p['x1'], FWHM=p['FWHM1'])
  y2 <- y0 + p['A2']*Lor(d$w, x0=p['x2'], FWHM=p['FWHM2'])
  y3 <- y0 + p['A3']*Lor(d$w, x0=p['x3'], FWHM=p['FWHM3'])
  ytot <- y0+y1+y2
  spec <- rbind(spec, tibble(w      = d$w,
                             Int    = d$Int, 
                             Int_n  = d$Int/max(d$Int), 
                             y1     = y1,
                             y1_n   = y1/max(d$Int),
                             y2     = y2,
                             y2_n   = y2/max(d$Int),
                             y3     = y3,
                             y3_n   = y3/max(d$Int),
                             ytot   = ytot,
                             ytot_n = ytot/max(d$Int),
                             name   = gsub(".txt","",norm_spls[i]),
                             P      = round(Pruby(coef(fit)['x1']),2)
  ))
  fitpar <- rbind(fitpar, data.frame(param = names(coef(fit)),
                                     fit   = coef(fit),
                                     name  = gsub(".txt","",norm_spls[i]),
                                     P     = round(Pruby(coef(fit)['x1']),2)
  ))
}

#PBS1_05mM_r2
PBS1_05mM_r2_d <- read_tsv("Data_means/PBS1_05mM_r2_means.txt", 
              col_names = c("w", "Int"))
PBS1_05mM_r2_d = PBS1_05mM_r2_d%>%filter(w>205 & w<290)
PBS1_05mM_r2_bg = baseline(PBS1_05mM_r2_d$w, PBS1_05mM_r2_d$Int)$bg
Guess <- c(y0 = 0.01, # constant background
           x = c(230,235,250), # positions
           FWHM = c(9, 10, 10), # full width at half maximum
           A = rep(max(PBS1_05mM_r2_d$Int)*2, 3))
PBS1_05mM_r2_fit <- nls(data =PBS1_05mM_r2_d, 
           Int-PBS1_05mM_r2_bg ~ y0 + 
             A1*Lor(w,x1,FWHM1) + 
             A2*Lor(w,x2,FWHM2) + 
             A3*Lor(w,x3,FWHM3), 
           start=as.list(Guess), 
           lower=as.list(Guess*0),
           algorithm = "port")
PBS1_05mM_r2_p  <- coef(PBS1_05mM_r2_fit)
PBS1_05mM_r2_y0 <- PBS1_05mM_r2_p['y0']
PBS1_05mM_r2_y1 <- PBS1_05mM_r2_y0 + 
                   PBS1_05mM_r2_p['A1']*Lor(PBS1_05mM_r2_d$w, 
                   x0=PBS1_05mM_r2_p['x1'], 
                   FWHM=PBS1_05mM_r2_p['FWHM1'])
PBS1_05mM_r2_y2 <- PBS1_05mM_r2_y0 + 
                   PBS1_05mM_r2_p['A2']*Lor(PBS1_05mM_r2_d$w, 
                   x0=PBS1_05mM_r2_p['x2'], 
                   FWHM=PBS1_05mM_r2_p['FWHM2'])
PBS1_05mM_r2_y3 <- PBS1_05mM_r2_y0 + 
                   PBS1_05mM_r2_p['A3']*Lor(PBS1_05mM_r2_d$w, 
                   x0=PBS1_05mM_r2_p['x3'], 
                   FWHM=PBS1_05mM_r2_p['FWHM3'])
PBS1_05mM_r2_ytot <- PBS1_05mM_r2_y0+PBS1_05mM_r2_y1+PBS1_05mM_r2_y2
spec <- rbind(spec, tibble(w      = PBS1_05mM_r2_d$w,
                           Int    = PBS1_05mM_r2_d$Int, 
                           Int_n  = PBS1_05mM_r2_d$Int/max(PBS1_05mM_r2_d$Int), 
                           y1     = PBS1_05mM_r2_y1,
                           y1_n   = PBS1_05mM_r2_y1/max(PBS1_05mM_r2_d$Int),
                           y2     = PBS1_05mM_r2_y2,
                           y2_n   = PBS1_05mM_r2_y2/max(PBS1_05mM_r2_d$Int),
                           y3     = PBS1_05mM_r2_y3,
                           y3_n   = PBS1_05mM_r2_y3/max(PBS1_05mM_r2_d$Int),
                           ytot   = PBS1_05mM_r2_ytot,
                           ytot_n = PBS1_05mM_r2_ytot/max(PBS1_05mM_r2_d$Int),
                           name   = gsub(".txt","","PBS1_05mM_r2_means.txt"),
                           P      = round(Pruby(coef(PBS1_05mM_r2_fit)['x1']),2)
))
fitpar <- rbind(fitpar, data.frame(param = names(coef(PBS1_05mM_r2_fit)),
                                   fit   = coef(PBS1_05mM_r2_fit),
                                   name  = gsub(".txt","","PBS1_05mM_r2_means.txt"),
                                   P     = round(Pruby(coef(PBS1_05mM_r2_fit)['x1']),2)
))

#PBS1_05mM_r3
PBS1_05mM_r3_d <- read_tsv("Data_means/PBS1_05mM_r3_means.txt", 
                           col_names = c("w", "Int"))
PBS1_05mM_r3_d = PBS1_05mM_r3_d%>%filter(w>205 & w<290)
PBS1_05mM_r3_bg = baseline(PBS1_05mM_r3_d$w, PBS1_05mM_r3_d$Int)$bg
Guess <- c(y0 = 0.01, # constant background
           x = c(230,235,250), # positions
           FWHM = c(9, 10, 10), # full width at half maximum
           A = rep(max(PBS1_05mM_r3_d$Int)*2, 3))
PBS1_05mM_r3_fit <- nls(data =PBS1_05mM_r3_d, 
                        Int-PBS1_05mM_r3_bg ~ y0 + 
                          A1*Lor(w,x1,FWHM1) + 
                          A2*Lor(w,x2,FWHM2) + 
                          A3*Lor(w,x3,FWHM3), 
                        start=as.list(Guess), 
                        lower=as.list(Guess*0),
                        algorithm = "port")
PBS1_05mM_r3_p  <- coef(PBS1_05mM_r3_fit)
PBS1_05mM_r3_y0 <- PBS1_05mM_r3_p['y0']
PBS1_05mM_r3_y1 <- PBS1_05mM_r3_y0 + 
  PBS1_05mM_r3_p['A1']*Lor(PBS1_05mM_r3_d$w, 
                           x0=PBS1_05mM_r3_p['x1'], 
                           FWHM=PBS1_05mM_r3_p['FWHM1'])
PBS1_05mM_r3_y2 <- PBS1_05mM_r3_y0 + 
  PBS1_05mM_r3_p['A2']*Lor(PBS1_05mM_r3_d$w, 
                           x0=PBS1_05mM_r3_p['x2'], 
                           FWHM=PBS1_05mM_r3_p['FWHM2'])
PBS1_05mM_r3_y3 <- PBS1_05mM_r3_y0 + 
  PBS1_05mM_r3_p['A3']*Lor(PBS1_05mM_r3_d$w, 
                           x0=PBS1_05mM_r3_p['x3'], 
                           FWHM=PBS1_05mM_r3_p['FWHM3'])
PBS1_05mM_r3_ytot <- PBS1_05mM_r3_y0+PBS1_05mM_r3_y1+PBS1_05mM_r3_y2
spec <- rbind(spec, tibble(w      = PBS1_05mM_r3_d$w,
                           Int    = PBS1_05mM_r3_d$Int, 
                           Int_n  = PBS1_05mM_r3_d$Int/max(PBS1_05mM_r3_d$Int), 
                           y1     = PBS1_05mM_r3_y1,
                           y1_n   = PBS1_05mM_r3_y1/max(PBS1_05mM_r3_d$Int),
                           y2     = PBS1_05mM_r3_y2,
                           y2_n   = PBS1_05mM_r3_y2/max(PBS1_05mM_r3_d$Int),
                           y3     = PBS1_05mM_r3_y3,
                           y3_n   = PBS1_05mM_r3_y3/max(PBS1_05mM_r3_d$Int),
                           ytot   = PBS1_05mM_r3_ytot,
                           ytot_n = PBS1_05mM_r3_ytot/max(PBS1_05mM_r3_d$Int),
                           name   = gsub(".txt","","PBS1_05mM_r3_means.txt"),
                           P      = round(Pruby(coef(PBS1_05mM_r3_fit)['x1']),2)
))
fitpar <- rbind(fitpar, data.frame(param = names(coef(PBS1_05mM_r3_fit)),
                                   fit   = coef(PBS1_05mM_r3_fit),
                                   name  = gsub(".txt","","PBS1_05mM_r3_means.txt"),
                                   P     = round(Pruby(coef(PBS1_05mM_r3_fit)['x1']),2)
))

#PBS1_1mM_r3
PBS1_1mM_r3_d <- read_tsv("Data_means/PBS1_1mM_r3_means.txt", 
                           col_names = c("w", "Int"))
PBS1_1mM_r3_d = PBS1_1mM_r3_d%>%filter(w>205 & w<290)
PBS1_1mM_r3_bg = baseline(PBS1_1mM_r3_d$w, PBS1_1mM_r3_d$Int)$bg
Guess <- c(y0 = 0.01, # constant background
           x = c(230,235,250), # positions
           FWHM = c(8, 10, 10), # full width at half maximum
           A = rep(max(PBS1_1mM_r3_d$Int)*2, 3))
PBS1_1mM_r3_fit <- nls(data =PBS1_1mM_r3_d, 
                        Int-PBS1_1mM_r3_bg ~ y0 + 
                          A1*Lor(w,x1,FWHM1) + 
                          A2*Lor(w,x2,FWHM2) + 
                          A3*Lor(w,x3,FWHM3), 
                        start=as.list(Guess), 
                        lower=as.list(Guess*0),
                        algorithm = "port")
PBS1_1mM_r3_p  <- coef(PBS1_1mM_r3_fit)
PBS1_1mM_r3_y0 <- PBS1_1mM_r3_p['y0']
PBS1_1mM_r3_y1 <- PBS1_1mM_r3_y0 + 
  PBS1_1mM_r3_p['A1']*Lor(PBS1_1mM_r3_d$w, 
                           x0=PBS1_1mM_r3_p['x1'], 
                           FWHM=PBS1_1mM_r3_p['FWHM1'])
PBS1_1mM_r3_y2 <- PBS1_1mM_r3_y0 + 
  PBS1_1mM_r3_p['A2']*Lor(PBS1_1mM_r3_d$w, 
                           x0=PBS1_1mM_r3_p['x2'], 
                           FWHM=PBS1_1mM_r3_p['FWHM2'])
PBS1_1mM_r3_y3 <- PBS1_1mM_r3_y0 + 
  PBS1_1mM_r3_p['A3']*Lor(PBS1_1mM_r3_d$w, 
                           x0=PBS1_1mM_r3_p['x3'], 
                           FWHM=PBS1_1mM_r3_p['FWHM3'])
PBS1_1mM_r3_ytot <- PBS1_1mM_r3_y0+PBS1_1mM_r3_y1+PBS1_1mM_r3_y2
spec <- rbind(spec, tibble(w      = PBS1_1mM_r3_d$w,
                           Int    = PBS1_1mM_r3_d$Int, 
                           Int_n  = PBS1_1mM_r3_d$Int/max(PBS1_1mM_r3_d$Int), 
                           y1     = PBS1_1mM_r3_y1,
                           y1_n   = PBS1_1mM_r3_y1/max(PBS1_1mM_r3_d$Int),
                           y2     = PBS1_1mM_r3_y2,
                           y2_n   = PBS1_1mM_r3_y2/max(PBS1_1mM_r3_d$Int),
                           y3     = PBS1_1mM_r3_y3,
                           y3_n   = PBS1_1mM_r3_y3/max(PBS1_1mM_r3_d$Int),
                           ytot   = PBS1_1mM_r3_ytot,
                           ytot_n = PBS1_1mM_r3_ytot/max(PBS1_1mM_r3_d$Int),
                           name   = gsub(".txt","","PBS1_1mM_r3_means.txt"),
                           P      = round(Pruby(coef(PBS1_1mM_r3_fit)['x1']),2)
))
fitpar <- rbind(fitpar, data.frame(param = names(coef(PBS1_1mM_r3_fit)),
                                   fit   = coef(PBS1_1mM_r3_fit),
                                   name  = gsub(".txt","","PBS1_1mM_r3_means.txt"),
                                   P     = round(Pruby(coef(PBS1_1mM_r3_fit)['x1']),2)
))

#PBS1_5mM_r2
PBS1_5mM_r2_d <- read_tsv("Data_means/PBS1_5mM_r2_means.txt", 
                           col_names = c("w", "Int"))
PBS1_5mM_r2_d = PBS1_5mM_r2_d%>%filter(w>205 & w<290)
PBS1_5mM_r2_bg = baseline(PBS1_5mM_r2_d$w, PBS1_5mM_r2_d$Int)$bg
Guess <- c(y0 = 0.01, # constant background
           x = c(230,235,250), # positions
           FWHM = c(10, 10, 9), # full width at half maximum
           A = rep(max(PBS1_5mM_r2_d$Int)*2, 3))
PBS1_5mM_r2_fit <- nls(data =PBS1_5mM_r2_d, 
                        Int-PBS1_5mM_r2_bg ~ y0 + 
                          A1*Lor(w,x1,FWHM1) + 
                          A2*Lor(w,x2,FWHM2) + 
                          A3*Lor(w,x3,FWHM3), 
                        start=as.list(Guess), 
                        lower=as.list(Guess*0),
                        algorithm = "port")
PBS1_5mM_r2_p  <- coef(PBS1_5mM_r2_fit)
PBS1_5mM_r2_y0 <- PBS1_5mM_r2_p['y0']
PBS1_5mM_r2_y1 <- PBS1_5mM_r2_y0 + 
  PBS1_5mM_r2_p['A1']*Lor(PBS1_5mM_r2_d$w, 
                           x0=PBS1_5mM_r2_p['x1'], 
                           FWHM=PBS1_5mM_r2_p['FWHM1'])
PBS1_5mM_r2_y2 <- PBS1_5mM_r2_y0 + 
  PBS1_5mM_r2_p['A2']*Lor(PBS1_5mM_r2_d$w, 
                           x0=PBS1_5mM_r2_p['x2'], 
                           FWHM=PBS1_5mM_r2_p['FWHM2'])
PBS1_5mM_r2_y3 <- PBS1_5mM_r2_y0 + 
  PBS1_5mM_r2_p['A3']*Lor(PBS1_5mM_r2_d$w, 
                           x0=PBS1_5mM_r2_p['x3'], 
                           FWHM=PBS1_5mM_r2_p['FWHM3'])
PBS1_5mM_r2_ytot <- PBS1_5mM_r2_y0+PBS1_5mM_r2_y1+PBS1_5mM_r2_y2
spec <- rbind(spec, tibble(w      = PBS1_5mM_r2_d$w,
                           Int    = PBS1_5mM_r2_d$Int, 
                           Int_n  = PBS1_5mM_r2_d$Int/max(PBS1_5mM_r2_d$Int), 
                           y1     = PBS1_5mM_r2_y1,
                           y1_n   = PBS1_5mM_r2_y1/max(PBS1_5mM_r2_d$Int),
                           y2     = PBS1_5mM_r2_y2,
                           y2_n   = PBS1_5mM_r2_y2/max(PBS1_5mM_r2_d$Int),
                           y3     = PBS1_5mM_r2_y3,
                           y3_n   = PBS1_5mM_r2_y3/max(PBS1_5mM_r2_d$Int),
                           ytot   = PBS1_5mM_r2_ytot,
                           ytot_n = PBS1_5mM_r2_ytot/max(PBS1_5mM_r2_d$Int),
                           name   = gsub(".txt","","PBS1_5mM_r2_means.txt"),
                           P      = round(Pruby(coef(PBS1_5mM_r2_fit)['x1']),2)
))
fitpar <- rbind(fitpar, data.frame(param = names(coef(PBS1_5mM_r2_fit)),
                                   fit   = coef(PBS1_5mM_r2_fit),
                                   name  = gsub(".txt","","PBS1_5mM_r2_means.txt"),
                                   P     = round(Pruby(coef(PBS1_5mM_r2_fit)['x1']),2)
))

#PBS1_5mM_r3
PBS1_5mM_r3_d <- read_tsv("Data_means/PBS1_5mM_r3_means.txt", 
                           col_names = c("w", "Int"))
PBS1_5mM_r3_d = PBS1_5mM_r3_d%>%filter(w>205 & w<290)
PBS1_5mM_r3_bg = baseline(PBS1_5mM_r3_d$w, PBS1_5mM_r3_d$Int)$bg
Guess <- c(y0 = 0.01, # constant background
           x = c(230,235,250), # positions
           FWHM = c(13, 10, 11), # full width at half maximum
           A = rep(max(PBS1_5mM_r3_d$Int)*2, 3))
PBS1_5mM_r3_fit <- nls(data =PBS1_5mM_r3_d, 
                        Int-PBS1_5mM_r3_bg ~ y0 + 
                          A1*Lor(w,x1,FWHM1) + 
                          A2*Lor(w,x2,FWHM2) + 
                          A3*Lor(w,x3,FWHM3), 
                        start=as.list(Guess), 
                        lower=as.list(Guess*0),
                        algorithm = "port")
PBS1_5mM_r3_p  <- coef(PBS1_5mM_r3_fit)
PBS1_5mM_r3_y0 <- PBS1_5mM_r3_p['y0']
PBS1_5mM_r3_y1 <- PBS1_5mM_r3_y0 + 
  PBS1_5mM_r3_p['A1']*Lor(PBS1_5mM_r3_d$w, 
                           x0=PBS1_5mM_r3_p['x1'], 
                           FWHM=PBS1_5mM_r3_p['FWHM1'])
PBS1_5mM_r3_y2 <- PBS1_5mM_r3_y0 + 
  PBS1_5mM_r3_p['A2']*Lor(PBS1_5mM_r3_d$w, 
                           x0=PBS1_5mM_r3_p['x2'], 
                           FWHM=PBS1_5mM_r3_p['FWHM2'])
PBS1_5mM_r3_y3 <- PBS1_5mM_r3_y0 + 
  PBS1_5mM_r3_p['A3']*Lor(PBS1_5mM_r3_d$w, 
                           x0=PBS1_5mM_r3_p['x3'], 
                           FWHM=PBS1_5mM_r3_p['FWHM3'])
PBS1_5mM_r3_ytot <- PBS1_5mM_r3_y0+PBS1_5mM_r3_y1+PBS1_5mM_r3_y2
spec <- rbind(spec, tibble(w      = PBS1_5mM_r3_d$w,
                           Int    = PBS1_5mM_r3_d$Int, 
                           Int_n  = PBS1_5mM_r3_d$Int/max(PBS1_5mM_r3_d$Int), 
                           y1     = PBS1_5mM_r3_y1,
                           y1_n   = PBS1_5mM_r3_y1/max(PBS1_5mM_r3_d$Int),
                           y2     = PBS1_5mM_r3_y2,
                           y2_n   = PBS1_5mM_r3_y2/max(PBS1_5mM_r3_d$Int),
                           y3     = PBS1_5mM_r3_y3,
                           y3_n   = PBS1_5mM_r3_y3/max(PBS1_5mM_r3_d$Int),
                           ytot   = PBS1_5mM_r3_ytot,
                           ytot_n = PBS1_5mM_r3_ytot/max(PBS1_5mM_r3_d$Int),
                           name   = gsub(".txt","","PBS1_5mM_r3_means.txt"),
                           P      = round(Pruby(coef(PBS1_5mM_r3_fit)['x1']),2)
))
fitpar <- rbind(fitpar, data.frame(param = names(coef(PBS1_5mM_r3_fit)),
                                   fit   = coef(PBS1_5mM_r3_fit),
                                   name  = gsub(".txt","","PBS1_5mM_r3_means.txt"),
                                   P     = round(Pruby(coef(PBS1_5mM_r3_fit)['x1']),2)
))

#PBS1_10mM_r2
PBS1_10mM_r2_d <- read_tsv("Data_means/PBS1_10mM_r2_means.txt", 
                           col_names = c("w", "Int"))
PBS1_10mM_r2_d = PBS1_10mM_r2_d%>%filter(w>205 & w<290)
PBS1_10mM_r2_bg = baseline(PBS1_10mM_r2_d$w, PBS1_10mM_r2_d$Int)$bg
Guess <- c(y0 = 0.01, # constant background
           x = c(230,235,250), # positions
           FWHM = c(16, 10, 13), # full width at half maximum
           A = rep(max(PBS1_10mM_r2_d$Int)*2, 3))
PBS1_10mM_r2_fit <- nls(data =PBS1_10mM_r2_d, 
                        Int-PBS1_10mM_r2_bg ~ y0 + 
                          A1*Lor(w,x1,FWHM1) + 
                          A2*Lor(w,x2,FWHM2) + 
                          A3*Lor(w,x3,FWHM3), 
                        start=as.list(Guess), 
                        lower=as.list(Guess*0),
                        algorithm = "port")
PBS1_10mM_r2_p  <- coef(PBS1_10mM_r2_fit)
PBS1_10mM_r2_y0 <- PBS1_10mM_r2_p['y0']
PBS1_10mM_r2_y1 <- PBS1_10mM_r2_y0 + 
  PBS1_10mM_r2_p['A1']*Lor(PBS1_10mM_r2_d$w, 
                           x0=PBS1_10mM_r2_p['x1'], 
                           FWHM=PBS1_10mM_r2_p['FWHM1'])
PBS1_10mM_r2_y2 <- PBS1_10mM_r2_y0 + 
  PBS1_10mM_r2_p['A2']*Lor(PBS1_10mM_r2_d$w, 
                           x0=PBS1_10mM_r2_p['x2'], 
                           FWHM=PBS1_10mM_r2_p['FWHM2'])
PBS1_10mM_r2_y3 <- PBS1_10mM_r2_y0 + 
  PBS1_10mM_r2_p['A3']*Lor(PBS1_10mM_r2_d$w, 
                           x0=PBS1_10mM_r2_p['x3'], 
                           FWHM=PBS1_10mM_r2_p['FWHM3'])
PBS1_10mM_r2_ytot <- PBS1_10mM_r2_y0+PBS1_10mM_r2_y1+PBS1_10mM_r2_y2
spec <- rbind(spec, tibble(w      = PBS1_10mM_r2_d$w,
                           Int    = PBS1_10mM_r2_d$Int, 
                           Int_n  = PBS1_10mM_r2_d$Int/max(PBS1_10mM_r2_d$Int), 
                           y1     = PBS1_10mM_r2_y1,
                           y1_n   = PBS1_10mM_r2_y1/max(PBS1_10mM_r2_d$Int),
                           y2     = PBS1_10mM_r2_y2,
                           y2_n   = PBS1_10mM_r2_y2/max(PBS1_10mM_r2_d$Int),
                           y3     = PBS1_10mM_r2_y3,
                           y3_n   = PBS1_10mM_r2_y3/max(PBS1_10mM_r2_d$Int),
                           ytot   = PBS1_10mM_r2_ytot,
                           ytot_n = PBS1_10mM_r2_ytot/max(PBS1_10mM_r2_d$Int),
                           name   = gsub(".txt","","PBS1_10mM_r2_means.txt"),
                           P      = round(Pruby(coef(PBS1_10mM_r2_fit)['x1']),2)
))
fitpar <- rbind(fitpar, data.frame(param = names(coef(PBS1_10mM_r2_fit)),
                                   fit   = coef(PBS1_10mM_r2_fit),
                                   name  = gsub(".txt","","PBS1_10mM_r2_means.txt"),
                                   P     = round(Pruby(coef(PBS1_10mM_r2_fit)['x1']),2)
))

#PBS1_10mM_r3
PBS1_10mM_r3_d <- read_tsv("Data_means/PBS1_10mM_r3_means.txt", 
                           col_names = c("w", "Int"))
PBS1_10mM_r3_d = PBS1_10mM_r3_d%>%filter(w>205 & w<290)
PBS1_10mM_r3_bg = baseline(PBS1_10mM_r3_d$w, PBS1_10mM_r3_d$Int)$bg
Guess <- c(y0 = 0.01, # constant background
           x = c(230,235,250), # positions
           FWHM = c(16, 10, 13), # full width at half maximum
           A = rep(max(PBS1_10mM_r3_d$Int)*2, 3))
PBS1_10mM_r3_fit <- nls(data =PBS1_10mM_r3_d, 
                        Int-PBS1_10mM_r3_bg ~ y0 + 
                          A1*Lor(w,x1,FWHM1) + 
                          A2*Lor(w,x2,FWHM2) + 
                          A3*Lor(w,x3,FWHM3), 
                        start=as.list(Guess), 
                        lower=as.list(Guess*0),
                        algorithm = "port")
PBS1_10mM_r3_p  <- coef(PBS1_10mM_r3_fit)
PBS1_10mM_r3_y0 <- PBS1_10mM_r3_p['y0']
PBS1_10mM_r3_y1 <- PBS1_10mM_r3_y0 + 
  PBS1_10mM_r3_p['A1']*Lor(PBS1_10mM_r3_d$w, 
                           x0=PBS1_10mM_r3_p['x1'], 
                           FWHM=PBS1_10mM_r3_p['FWHM1'])
PBS1_10mM_r3_y2 <- PBS1_10mM_r3_y0 + 
  PBS1_10mM_r3_p['A2']*Lor(PBS1_10mM_r3_d$w, 
                           x0=PBS1_10mM_r3_p['x2'], 
                           FWHM=PBS1_10mM_r3_p['FWHM2'])
PBS1_10mM_r3_y3 <- PBS1_10mM_r3_y0 + 
  PBS1_10mM_r3_p['A3']*Lor(PBS1_10mM_r3_d$w, 
                           x0=PBS1_10mM_r3_p['x3'], 
                           FWHM=PBS1_10mM_r3_p['FWHM3'])
PBS1_10mM_r3_ytot <- PBS1_10mM_r3_y0+PBS1_10mM_r3_y1+PBS1_10mM_r3_y2
spec <- rbind(spec, tibble(w      = PBS1_10mM_r3_d$w,
                           Int    = PBS1_10mM_r3_d$Int, 
                           Int_n  = PBS1_10mM_r3_d$Int/max(PBS1_10mM_r3_d$Int), 
                           y1     = PBS1_10mM_r3_y1,
                           y1_n   = PBS1_10mM_r3_y1/max(PBS1_10mM_r3_d$Int),
                           y2     = PBS1_10mM_r3_y2,
                           y2_n   = PBS1_10mM_r3_y2/max(PBS1_10mM_r3_d$Int),
                           y3     = PBS1_10mM_r3_y3,
                           y3_n   = PBS1_10mM_r3_y3/max(PBS1_10mM_r3_d$Int),
                           ytot   = PBS1_10mM_r3_ytot,
                           ytot_n = PBS1_10mM_r3_ytot/max(PBS1_10mM_r3_d$Int),
                           name   = gsub(".txt","","PBS1_10mM_r3_means.txt"),
                           P      = round(Pruby(coef(PBS1_10mM_r3_fit)['x1']),2)
))
fitpar <- rbind(fitpar, data.frame(param = names(coef(PBS1_10mM_r3_fit)),
                                   fit   = coef(PBS1_10mM_r3_fit),
                                   name  = gsub(".txt","","PBS1_10mM_r3_means.txt"),
                                   P     = round(Pruby(coef(PBS1_10mM_r3_fit)['x1']),2)
))

#PBH481_05mM_r2
PBH481_05mM_r2_d <- read_tsv("Data_means/PBH481_05mM_r2_means.txt", 
                           col_names = c("w", "Int"))
PBH481_05mM_r2_d = PBH481_05mM_r2_d%>%filter(w>205 & w<290)
PBH481_05mM_r2_bg = baseline(PBH481_05mM_r2_d$w, PBH481_05mM_r2_d$Int)$bg
Guess <- c(y0 = 0.01, # constant background
           x = c(230,235,250), # positions
           FWHM = c(12, 10, 10), # full width at half maximum
           A = rep(max(PBH481_05mM_r2_d$Int)*2, 3))
PBH481_05mM_r2_fit <- nls(data =PBH481_05mM_r2_d, 
                        Int-PBH481_05mM_r2_bg ~ y0 + 
                          A1*Lor(w,x1,FWHM1) + 
                          A2*Lor(w,x2,FWHM2) + 
                          A3*Lor(w,x3,FWHM3), 
                        start=as.list(Guess), 
                        lower=as.list(Guess*0),
                        algorithm = "port")
PBH481_05mM_r2_p  <- coef(PBH481_05mM_r2_fit)
PBH481_05mM_r2_y0 <- PBH481_05mM_r2_p['y0']
PBH481_05mM_r2_y1 <- PBH481_05mM_r2_y0 + 
  PBH481_05mM_r2_p['A1']*Lor(PBH481_05mM_r2_d$w, 
                           x0=PBH481_05mM_r2_p['x1'], 
                           FWHM=PBH481_05mM_r2_p['FWHM1'])
PBH481_05mM_r2_y2 <- PBH481_05mM_r2_y0 + 
  PBH481_05mM_r2_p['A2']*Lor(PBH481_05mM_r2_d$w, 
                           x0=PBH481_05mM_r2_p['x2'], 
                           FWHM=PBH481_05mM_r2_p['FWHM2'])
PBH481_05mM_r2_y3 <- PBH481_05mM_r2_y0 + 
  PBH481_05mM_r2_p['A3']*Lor(PBH481_05mM_r2_d$w, 
                           x0=PBH481_05mM_r2_p['x3'], 
                           FWHM=PBH481_05mM_r2_p['FWHM3'])
PBH481_05mM_r2_ytot <- PBH481_05mM_r2_y0+PBH481_05mM_r2_y1+PBH481_05mM_r2_y2
spec <- rbind(spec, tibble(w      = PBH481_05mM_r2_d$w,
                           Int    = PBH481_05mM_r2_d$Int, 
                           Int_n  = PBH481_05mM_r2_d$Int/max(PBH481_05mM_r2_d$Int), 
                           y1     = PBH481_05mM_r2_y1,
                           y1_n   = PBH481_05mM_r2_y1/max(PBH481_05mM_r2_d$Int),
                           y2     = PBH481_05mM_r2_y2,
                           y2_n   = PBH481_05mM_r2_y2/max(PBH481_05mM_r2_d$Int),
                           y3     = PBH481_05mM_r2_y3,
                           y3_n   = PBH481_05mM_r2_y3/max(PBH481_05mM_r2_d$Int),
                           ytot   = PBH481_05mM_r2_ytot,
                           ytot_n = PBH481_05mM_r2_ytot/max(PBH481_05mM_r2_d$Int),
                           name   = gsub(".txt","","PBH481_05mM_r2_means.txt"),
                           P      = round(Pruby(coef(PBH481_05mM_r2_fit)['x1']),2)
))
fitpar <- rbind(fitpar, data.frame(param = names(coef(PBH481_05mM_r2_fit)),
                                   fit   = coef(PBH481_05mM_r2_fit),
                                   name  = gsub(".txt","","PBH481_05mM_r2_means.txt"),
                                   P     = round(Pruby(coef(PBH481_05mM_r2_fit)['x1']),2)
))

#PBH481_05mM_r3
PBH481_05mM_r3_d <- read_tsv("Data_means/PBH481_05mM_r3_means.txt", 
                             col_names = c("w", "Int"))
PBH481_05mM_r3_d = PBH481_05mM_r3_d%>%filter(w>205 & w<290)
PBH481_05mM_r3_bg = baseline(PBH481_05mM_r3_d$w, PBH481_05mM_r3_d$Int)$bg
Guess <- c(y0 = 0.01, # constant background
           x = c(230,235,250), # positions
           FWHM = c(16, 10, 13), # full width at half maximum
           A = rep(max(PBH481_05mM_r3_d$Int)*2, 3))
PBH481_05mM_r3_fit <- nls(data =PBH481_05mM_r3_d, 
                          Int-PBH481_05mM_r3_bg ~ y0 + 
                            A1*Lor(w,x1,FWHM1) + 
                            A2*Lor(w,x2,FWHM2) + 
                            A3*Lor(w,x3,FWHM3), 
                          start=as.list(Guess), 
                          lower=as.list(Guess*0),
                          algorithm = "port")
PBH481_05mM_r3_p  <- coef(PBH481_05mM_r3_fit)
PBH481_05mM_r3_y0 <- PBH481_05mM_r3_p['y0']
PBH481_05mM_r3_y1 <- PBH481_05mM_r3_y0 + 
  PBH481_05mM_r3_p['A1']*Lor(PBH481_05mM_r3_d$w, 
                             x0=PBH481_05mM_r3_p['x1'], 
                             FWHM=PBH481_05mM_r3_p['FWHM1'])
PBH481_05mM_r3_y2 <- PBH481_05mM_r3_y0 + 
  PBH481_05mM_r3_p['A2']*Lor(PBH481_05mM_r3_d$w, 
                             x0=PBH481_05mM_r3_p['x2'], 
                             FWHM=PBH481_05mM_r3_p['FWHM2'])
PBH481_05mM_r3_y3 <- PBH481_05mM_r3_y0 + 
  PBH481_05mM_r3_p['A3']*Lor(PBH481_05mM_r3_d$w, 
                             x0=PBH481_05mM_r3_p['x3'], 
                             FWHM=PBH481_05mM_r3_p['FWHM3'])
PBH481_05mM_r3_ytot <- PBH481_05mM_r3_y0+PBH481_05mM_r3_y1+PBH481_05mM_r3_y2
spec <- rbind(spec, tibble(w      = PBH481_05mM_r3_d$w,
                           Int    = PBH481_05mM_r3_d$Int, 
                           Int_n  = PBH481_05mM_r3_d$Int/max(PBH481_05mM_r3_d$Int), 
                           y1     = PBH481_05mM_r3_y1,
                           y1_n   = PBH481_05mM_r3_y1/max(PBH481_05mM_r3_d$Int),
                           y2     = PBH481_05mM_r3_y2,
                           y2_n   = PBH481_05mM_r3_y2/max(PBH481_05mM_r3_d$Int),
                           y3     = PBH481_05mM_r3_y3,
                           y3_n   = PBH481_05mM_r3_y3/max(PBH481_05mM_r3_d$Int),
                           ytot   = PBH481_05mM_r3_ytot,
                           ytot_n = PBH481_05mM_r3_ytot/max(PBH481_05mM_r3_d$Int),
                           name   = gsub(".txt","","PBH481_05mM_r3_means.txt"),
                           P      = round(Pruby(coef(PBH481_05mM_r3_fit)['x1']),2)
))
fitpar <- rbind(fitpar, data.frame(param = names(coef(PBH481_05mM_r3_fit)),
                                   fit   = coef(PBH481_05mM_r3_fit),
                                   name  = gsub(".txt","","PBH481_05mM_r3_means.txt"),
                                   P     = round(Pruby(coef(PBH481_05mM_r3_fit)['x1']),2)
))

#PBH481_5mM_r2
PBH481_5mM_r2_d <- read_tsv("Data_means/PBH481_5mM_r2_means.txt", 
                             col_names = c("w", "Int"))
PBH481_5mM_r2_d = PBH481_5mM_r2_d%>%filter(w>205 & w<290)
PBH481_5mM_r2_bg = baseline(PBH481_5mM_r2_d$w, PBH481_5mM_r2_d$Int)$bg
Guess <- c(y0 = 0.01, # constant background
           x = c(230,235,250), # positions
           FWHM = c(16, 10, 13), # full width at half maximum
           A = rep(max(PBH481_5mM_r2_d$Int)*2, 3))
PBH481_5mM_r2_fit <- nls(data =PBH481_5mM_r2_d, 
                          Int-PBH481_5mM_r2_bg ~ y0 + 
                            A1*Lor(w,x1,FWHM1) + 
                            A2*Lor(w,x2,FWHM2) + 
                            A3*Lor(w,x3,FWHM3), 
                          start=as.list(Guess), 
                          lower=as.list(Guess*0),
                          algorithm = "port")
PBH481_5mM_r2_p  <- coef(PBH481_5mM_r2_fit)
PBH481_5mM_r2_y0 <- PBH481_5mM_r2_p['y0']
PBH481_5mM_r2_y1 <- PBH481_5mM_r2_y0 + 
  PBH481_5mM_r2_p['A1']*Lor(PBH481_5mM_r2_d$w, 
                             x0=PBH481_5mM_r2_p['x1'], 
                             FWHM=PBH481_5mM_r2_p['FWHM1'])
PBH481_5mM_r2_y2 <- PBH481_5mM_r2_y0 + 
  PBH481_5mM_r2_p['A2']*Lor(PBH481_5mM_r2_d$w, 
                             x0=PBH481_5mM_r2_p['x2'], 
                             FWHM=PBH481_5mM_r2_p['FWHM2'])
PBH481_5mM_r2_y3 <- PBH481_5mM_r2_y0 + 
  PBH481_5mM_r2_p['A3']*Lor(PBH481_5mM_r2_d$w, 
                             x0=PBH481_5mM_r2_p['x3'], 
                             FWHM=PBH481_5mM_r2_p['FWHM3'])
PBH481_5mM_r2_ytot <- PBH481_5mM_r2_y0+PBH481_5mM_r2_y1+PBH481_5mM_r2_y2
spec <- rbind(spec, tibble(w      = PBH481_5mM_r2_d$w,
                           Int    = PBH481_5mM_r2_d$Int, 
                           Int_n  = PBH481_5mM_r2_d$Int/max(PBH481_5mM_r2_d$Int), 
                           y1     = PBH481_5mM_r2_y1,
                           y1_n   = PBH481_5mM_r2_y1/max(PBH481_5mM_r2_d$Int),
                           y2     = PBH481_5mM_r2_y2,
                           y2_n   = PBH481_5mM_r2_y2/max(PBH481_5mM_r2_d$Int),
                           y3     = PBH481_5mM_r2_y3,
                           y3_n   = PBH481_5mM_r2_y3/max(PBH481_5mM_r2_d$Int),
                           ytot   = PBH481_5mM_r2_ytot,
                           ytot_n = PBH481_5mM_r2_ytot/max(PBH481_5mM_r2_d$Int),
                           name   = gsub(".txt","","PBH481_5mM_r2_means.txt"),
                           P      = round(Pruby(coef(PBH481_5mM_r2_fit)['x1']),2)
))
fitpar <- rbind(fitpar, data.frame(param = names(coef(PBH481_5mM_r2_fit)),
                                   fit   = coef(PBH481_5mM_r2_fit),
                                   name  = gsub(".txt","","PBH481_5mM_r2_means.txt"),
                                   P     = round(Pruby(coef(PBH481_5mM_r2_fit)['x1']),2)
))

#PBH481_5mM_r3
PBH481_5mM_r3_d <- read_tsv("Data_means/PBH481_5mM_r3_means.txt", 
                             col_names = c("w", "Int"))
PBH481_5mM_r3_d = PBH481_5mM_r3_d%>%filter(w>205 & w<290)
PBH481_5mM_r3_bg = baseline(PBH481_5mM_r3_d$w, PBH481_5mM_r3_d$Int)$bg
Guess <- c(y0 = 0.01, # constant background
           x = c(230,235,250), # positions
           FWHM = c(16, 10, 13), # full width at half maximum
           A = rep(max(PBH481_5mM_r3_d$Int)*2, 3))
PBH481_5mM_r3_fit <- nls(data =PBH481_5mM_r3_d, 
                          Int-PBH481_5mM_r3_bg ~ y0 + 
                            A1*Lor(w,x1,FWHM1) + 
                            A2*Lor(w,x2,FWHM2) + 
                            A3*Lor(w,x3,FWHM3), 
                          start=as.list(Guess), 
                          lower=as.list(Guess*0),
                          algorithm = "port")
PBH481_5mM_r3_p  <- coef(PBH481_5mM_r3_fit)
PBH481_5mM_r3_y0 <- PBH481_5mM_r3_p['y0']
PBH481_5mM_r3_y1 <- PBH481_5mM_r3_y0 + 
  PBH481_5mM_r3_p['A1']*Lor(PBH481_5mM_r3_d$w, 
                             x0=PBH481_5mM_r3_p['x1'], 
                             FWHM=PBH481_5mM_r3_p['FWHM1'])
PBH481_5mM_r3_y2 <- PBH481_5mM_r3_y0 + 
  PBH481_5mM_r3_p['A2']*Lor(PBH481_5mM_r3_d$w, 
                             x0=PBH481_5mM_r3_p['x2'], 
                             FWHM=PBH481_5mM_r3_p['FWHM2'])
PBH481_5mM_r3_y3 <- PBH481_5mM_r3_y0 + 
  PBH481_5mM_r3_p['A3']*Lor(PBH481_5mM_r3_d$w, 
                             x0=PBH481_5mM_r3_p['x3'], 
                             FWHM=PBH481_5mM_r3_p['FWHM3'])
PBH481_5mM_r3_ytot <- PBH481_5mM_r3_y0+PBH481_5mM_r3_y1+PBH481_5mM_r3_y2
spec <- rbind(spec, tibble(w      = PBH481_5mM_r3_d$w,
                           Int    = PBH481_5mM_r3_d$Int, 
                           Int_n  = PBH481_5mM_r3_d$Int/max(PBH481_5mM_r3_d$Int), 
                           y1     = PBH481_5mM_r3_y1,
                           y1_n   = PBH481_5mM_r3_y1/max(PBH481_5mM_r3_d$Int),
                           y2     = PBH481_5mM_r3_y2,
                           y2_n   = PBH481_5mM_r3_y2/max(PBH481_5mM_r3_d$Int),
                           y3     = PBH481_5mM_r3_y3,
                           y3_n   = PBH481_5mM_r3_y3/max(PBH481_5mM_r3_d$Int),
                           ytot   = PBH481_5mM_r3_ytot,
                           ytot_n = PBH481_5mM_r3_ytot/max(PBH481_5mM_r3_d$Int),
                           name   = gsub(".txt","","PBH481_5mM_r3_means.txt"),
                           P      = round(Pruby(coef(PBH481_5mM_r3_fit)['x1']),2)
))
fitpar <- rbind(fitpar, data.frame(param = names(coef(PBH481_5mM_r3_fit)),
                                   fit   = coef(PBH481_5mM_r3_fit),
                                   name  = gsub(".txt","","PBH481_5mM_r3_means.txt"),
                                   P     = round(Pruby(coef(PBH481_5mM_r3_fit)['x1']),2)
))

#PBH481_10mM_r2
PBH481_10mM_r2_d <- read_tsv("Data_means/PBH481_10mM_r2_means.txt", 
                            col_names = c("w", "Int"))
PBH481_10mM_r2_d = PBH481_10mM_r2_d%>%filter(w>205 & w<290)
PBH481_10mM_r2_bg = baseline(PBH481_10mM_r2_d$w, PBH481_10mM_r2_d$Int)$bg
Guess <- c(y0 = 0.01, # constant background
           x = c(230,235,250), # positions
           FWHM = c(13, 10, 14), # full width at half maximum
           A = rep(max(PBH481_10mM_r2_d$Int)*2, 3))
PBH481_10mM_r2_fit <- nls(data =PBH481_10mM_r2_d, 
                         Int-PBH481_10mM_r2_bg ~ y0 + 
                           A1*Lor(w,x1,FWHM1) + 
                           A2*Lor(w,x2,FWHM2) + 
                           A3*Lor(w,x3,FWHM3), 
                         start=as.list(Guess), 
                         lower=as.list(Guess*0),
                         algorithm = "port")
PBH481_10mM_r2_p  <- coef(PBH481_10mM_r2_fit)
PBH481_10mM_r2_y0 <- PBH481_10mM_r2_p['y0']
PBH481_10mM_r2_y1 <- PBH481_10mM_r2_y0 + 
  PBH481_10mM_r2_p['A1']*Lor(PBH481_10mM_r2_d$w, 
                            x0=PBH481_10mM_r2_p['x1'], 
                            FWHM=PBH481_10mM_r2_p['FWHM1'])
PBH481_10mM_r2_y2 <- PBH481_10mM_r2_y0 + 
  PBH481_10mM_r2_p['A2']*Lor(PBH481_10mM_r2_d$w, 
                            x0=PBH481_10mM_r2_p['x2'], 
                            FWHM=PBH481_10mM_r2_p['FWHM2'])
PBH481_10mM_r2_y3 <- PBH481_10mM_r2_y0 + 
  PBH481_10mM_r2_p['A3']*Lor(PBH481_10mM_r2_d$w, 
                            x0=PBH481_10mM_r2_p['x3'], 
                            FWHM=PBH481_10mM_r2_p['FWHM3'])
PBH481_10mM_r2_ytot <- PBH481_10mM_r2_y0+PBH481_10mM_r2_y1+PBH481_10mM_r2_y2
spec <- rbind(spec, tibble(w      = PBH481_10mM_r2_d$w,
                           Int    = PBH481_10mM_r2_d$Int, 
                           Int_n  = PBH481_10mM_r2_d$Int/max(PBH481_10mM_r2_d$Int), 
                           y1     = PBH481_10mM_r2_y1,
                           y1_n   = PBH481_10mM_r2_y1/max(PBH481_10mM_r2_d$Int),
                           y2     = PBH481_10mM_r2_y2,
                           y2_n   = PBH481_10mM_r2_y2/max(PBH481_10mM_r2_d$Int),
                           y3     = PBH481_10mM_r2_y3,
                           y3_n   = PBH481_10mM_r2_y3/max(PBH481_10mM_r2_d$Int),
                           ytot   = PBH481_10mM_r2_ytot,
                           ytot_n = PBH481_10mM_r2_ytot/max(PBH481_10mM_r2_d$Int),
                           name   = gsub(".txt","","PBH481_10mM_r2_means.txt"),
                           P      = round(Pruby(coef(PBH481_10mM_r2_fit)['x1']),2)
))
fitpar <- rbind(fitpar, data.frame(param = names(coef(PBH481_10mM_r2_fit)),
                                   fit   = coef(PBH481_10mM_r2_fit),
                                   name  = gsub(".txt","","PBH481_10mM_r2_means.txt"),
                                   P     = round(Pruby(coef(PBH481_10mM_r2_fit)['x1']),2)
))

#PBH481_10mM_r3
PBH481_10mM_r3_d <- read_tsv("Data_means/PBH481_10mM_r3_means.txt", 
                             col_names = c("w", "Int"))
PBH481_10mM_r3_d = PBH481_10mM_r3_d%>%filter(w>205 & w<290)
PBH481_10mM_r3_bg = baseline(PBH481_10mM_r3_d$w, PBH481_10mM_r3_d$Int)$bg
Guess <- c(y0 = 0.01, # constant background
           x = c(230,235,250), # positions
           FWHM = c(13, 10, 13), # full width at half maximum
           A = rep(max(PBH481_10mM_r3_d$Int)*2, 3))
PBH481_10mM_r3_fit <- nls(data =PBH481_10mM_r3_d, 
                          Int-PBH481_10mM_r3_bg ~ y0 + 
                            A1*Lor(w,x1,FWHM1) + 
                            A2*Lor(w,x2,FWHM2) + 
                            A3*Lor(w,x3,FWHM3), 
                          start=as.list(Guess), 
                          lower=as.list(Guess*0),
                          algorithm = "port")
PBH481_10mM_r3_p  <- coef(PBH481_10mM_r3_fit)
PBH481_10mM_r3_y0 <- PBH481_10mM_r3_p['y0']
PBH481_10mM_r3_y1 <- PBH481_10mM_r3_y0 + 
  PBH481_10mM_r3_p['A1']*Lor(PBH481_10mM_r3_d$w, 
                             x0=PBH481_10mM_r3_p['x1'], 
                             FWHM=PBH481_10mM_r3_p['FWHM1'])
PBH481_10mM_r3_y2 <- PBH481_10mM_r3_y0 + 
  PBH481_10mM_r3_p['A2']*Lor(PBH481_10mM_r3_d$w, 
                             x0=PBH481_10mM_r3_p['x2'], 
                             FWHM=PBH481_10mM_r3_p['FWHM2'])
PBH481_10mM_r3_y3 <- PBH481_10mM_r3_y0 + 
  PBH481_10mM_r3_p['A3']*Lor(PBH481_10mM_r3_d$w, 
                             x0=PBH481_10mM_r3_p['x3'], 
                             FWHM=PBH481_10mM_r3_p['FWHM3'])
PBH481_10mM_r3_ytot <- PBH481_10mM_r3_y0+PBH481_10mM_r3_y1+PBH481_10mM_r3_y2
spec <- rbind(spec, tibble(w      = PBH481_10mM_r3_d$w,
                           Int    = PBH481_10mM_r3_d$Int, 
                           Int_n  = PBH481_10mM_r3_d$Int/max(PBH481_10mM_r3_d$Int), 
                           y1     = PBH481_10mM_r3_y1,
                           y1_n   = PBH481_10mM_r3_y1/max(PBH481_10mM_r3_d$Int),
                           y2     = PBH481_10mM_r3_y2,
                           y2_n   = PBH481_10mM_r3_y2/max(PBH481_10mM_r3_d$Int),
                           y3     = PBH481_10mM_r3_y3,
                           y3_n   = PBH481_10mM_r3_y3/max(PBH481_10mM_r3_d$Int),
                           ytot   = PBH481_10mM_r3_ytot,
                           ytot_n = PBH481_10mM_r3_ytot/max(PBH481_10mM_r3_d$Int),
                           name   = gsub(".txt","","PBH481_10mM_r3_means.txt"),
                           P      = round(Pruby(coef(PBH481_10mM_r3_fit)['x1']),2)
))
fitpar <- rbind(fitpar, data.frame(param = names(coef(PBH481_10mM_r3_fit)),
                                   fit   = coef(PBH481_10mM_r3_fit),
                                   name  = gsub(".txt","","PBH481_10mM_r3_means.txt"),
                                   P     = round(Pruby(coef(PBH481_10mM_r3_fit)['x1']),2)
))