library(ggseqlogo)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
print(args)


plot_seq_logo <- function(d, bkg = NA, max = NA, title = "", left_side =13, right_side = 5) {
    if(length(bkg)<2){
        p <- ggseqlogo(d, seq_type = 'dna')
        p <- p + theme_classic()
    } else {
        
        bg_cor <- tibble(seq = d) %>%
            mutate(seq = str_split(seq, "", n = left_side + right_side, simplify = TRUE)) %>%
        as.matrix() %>% as_tibble() %>% 
        pivot_longer(everything(), values_to = "base") %>% 
        filter(base != "N") %>%
        mutate(pos = as.numeric(str_remove(name, "seq."))) %>%
        group_by(pos, base) %>% summarize(n=n(), .groups = "drop") %>% 
        group_by(pos) %>% 
        mutate(freq = n/sum(n), 
               total_inf = 2--sum(freq*log2(freq)),
               entropy = freq*total_inf) %>%
        left_join(bkg, by = "base") %>% 
        mutate(rel_ent = freq* log2(freq/prior)) %>%
        mutate(height = rel_ent) %>% select(pos, base, height) %>%
        pivot_wider(names_from = pos, values_from = height) %>%
        column_to_rownames("base") %>%
        as.matrix()
        p <- ggseqlogo(bg_cor, method = "custom", seq_type = 'dna')
        p <- p + theme_classic() + labs(y = "Relative Entropy")
    }
    p <- p + scale_y_continuous(limits = c(0, max)) + labs(title = title)
    if(left_side == 13 & right_side == 5){
    p <- p + annotate(geom = "rect", xmin = -8.5 + 13, xmax = -0.5 + 13, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "red") + 
    annotate(geom = "rect", xmin= -10.5 + 13, xmax = -8.5 +13, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "blue") +
    annotate(geom = "rect", xmin = -0.5 + 13, xmax = 1.5 + 13, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "green") +
    annotate(geom = "rect", xmin = 1.5 + 13, xmax = 5.5 + 13, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "purple")
    }
    p$layers <- rev(p$layers)
    p$scales$scales[[1]] <- scale_x_continuous(breaks = seq(1, left_side+ right_side, by = 1), 
                                           labels= c(seq(-left_side,-1,by=1), seq(1, right_side, by = 1)))
    list(p)
}

d <- read_tsv(args[1], col_names = c("seq"))
prior <- read_delim(args[3], col_names = c("base", "prior"), comment = "#", delim = " ")
if(length(args) > 3){
    left_side <- as.numeric(args[4]) + 1
    right_side <- as.numeric(args[5])
}else{
    left_side <- 13
    right_side <- 5
}

# Manually calculate the PWM
#bg_cor <- d %>% mutate(seq = str_split(seq, "", n = 18, simplify = TRUE)) %>% 
#    as.matrix() %>% as_tibble() %>% pivot_longer(everything()) %>% filter(value != "N") %>%
#    mutate(name = as.numeric(str_remove(name, "seq."))) %>%
#    group_by(name, value) %>% summarize(n=n()) %>% 
#    group_by(name) %>% mutate(freq = n/sum(n)) %>% 
#    left_join(prior, by = "value") %>% 
#    group_by(name) %>% 
#    mutate(total_inf = sum(freq*log2(freq/prior))) %>% ungroup() %>% 
#    mutate(rel_ent = pmax(freq* log2(freq/prior), 0)) %>%
#    mutate(height = rel_ent) %>% select(name, value, height) %>% 
#    pivot_wider(names_from = name, values_from = height) %>%
#    column_to_rownames("value") %>%
#    as.matrix()
#
#p <- ggseqlogo(d, seq_type = 'dna')
#p <- p + theme_classic() 
#p$scales$scales[[1]] <- scale_x_continuous(breaks = seq(1, 18, by = 1), 
#                                           labels= c(seq(-13,-1,by=1), seq(1, 5, by = 1)))

p1 <- plot_seq_logo(d %>% pull(seq), left_side = left_side, right_side = right_side)
ggsave(str_c(args[2], ".pdf"),  p1[[1]], width = min(20, round((left_side + right_side)/3.6)), height = 2)

#p <- ggseqlogo(bg_cor, method = "custom", seq_type = 'dna')
#p <- p + theme_classic() + labs(y = "Relative Entropy")
#p$scales$scales[[1]] <- scale_x_continuous(breaks = seq(1, 18, by = 1), 
#                                           labels= c(seq(-13,-1,by=1), seq(1, 5, by = 1)))





p2 <- plot_seq_logo(d %>% pull(seq), bkg = prior, left_side = left_side, right_side = right_side)
ggsave(str_c(args[2], "_bg_corrected.pdf"), p2[[1]], width=min(round((left_side + right_side)/3.6), 20), height = 2)
