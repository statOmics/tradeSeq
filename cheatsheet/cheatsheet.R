library(tidyverse)
library(knitr)
library(kableExtra)
theme_set(new = theme_classic())
theme_update(axis.text = element_text(size = 12),
             axis.line = element_line(size = 2),
             axis.title = element_text(size = 15))
if (!file.exists("figure")) dir.create(file.path("figures"))

labels <- labs(x = "pseudotime", y = "count (log + 1 scale)")
im1 <- data.frame(time1 = 1:100, fit1 = 1 - exp(-(1:100)/10),
                  time2 = seq(1,80,length.out = 100),
                  fit2 = 1 - exp(-(1:100)/10) + 1/40 *
                    cos(seq(1, 80, length.out = 100)/10))
ggsave(filename = "figures/im1.pdf",
       ggplot(im1) +
         geom_line(aes(x = time1, y = fit1), col = "darkgreen", size = 6) + 
         geom_line(aes(x = time2, y = fit2), col = "brown", size = 6) +
         lims(y = c(0, 1.2)) + labels)


im2 <- data.frame(time1 = 1:100, fit1 = 1 - exp(-(1:100)/10),
                  time2 = seq(1,80,length.out = 100),
                  fit2 = 1 - exp(-1/10) + 1/40 *
                    cos(seq(1,80,length.out = 100)/10))

ggsave(filename = "figures/im2.pdf",
       ggplot(im2) +
         geom_line(aes(x = time1, y = fit1), col = "darkgreen", size = 6) + 
         geom_line(aes(x = time2, y = fit2), col = "brown", size = 6) +
         lims(y = c(0, 1.2)) )

im3 <- data.frame(time1 = 1:100,
                  fit1 = 0.5 + 1/2 * cos(2*pi*1:100/100),
                  time2 = seq(1,80,length.out = 100),
                  fit2 = 0.5  + 1/2 *
                    cos(2 * pi * seq(1, 80, length.out = 100)/80))

ggsave(filename = "figures/im3.pdf",
       ggplot(im3) +
         geom_line(aes(x = time1, y = fit1), col = "darkgreen", size = 6) + 
         geom_line(aes(x = time2, y = fit2), col = "brown", size = 6) +
         lims(y = c(0,1.2)) + labels)

im4 <- data.frame(time1 = 1:100, fit1 = 1 - exp(-(1:100)/20),
                  time2 = seq(1,80,length.out = 100),
                  fit2 = 1/2 - 1/2 * exp(-seq(1,80,length.out = 100)/20))
ggsave(filename = "figures/im4.pdf",
       ggplot(im4) +
         geom_line(aes(x = time1, y = fit1), col = "darkgreen", size = 6) + 
         geom_line(aes(x = time2, y = fit2), col = "brown", size = 6) +
         lims(y = c(0,1.2)) + labels )

im5 <- data.frame(time1 = 1:100,
                  fit1 = .5 + .5 * cos(pi*0:99/100),
                  time2 = seq(1,80,length.out = 100),
                  fit2 = exp(-seq(0,80,length.out = 100)/10))  
ggsave(filename = "figures/im5.pdf",
       ggplot(im5) +
         geom_line(aes(x = time1, y = fit1), col = "darkgreen", size = 6) + 
         geom_line(aes(x = time2, y = fit2), col = "brown", size = 6) +
         lims(y = c(0,1.2)) + labels )

im6 <- data.frame(time1 = 1:100,
                  fit1 = .5 + .5 * cos(pi*0:99/100),
                  time2 = c(1:40, seq(41,80,length.out = 60)),
                  fit2 = c(.5 + .5 * cos(pi*0:39/100) - .01,
                           (.5 + .5 * cos(pi * 0.39)) *
                             exp(-seq(1,39,length.out = 60)/10)))  
ggsave(filename = "figures/im6.pdf",
       ggplot(im6) +
         geom_line(aes(x = time1, y = fit1), col = "darkgreen", size = 6) + 
         geom_line(aes(x = time2, y = fit2), col = "brown", size = 6) +
         lims(y = c(0,1.2)) + labels )

