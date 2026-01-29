By Jingxuan Wang
library(dplyr)
library(ggplot2)
library(GGally)
library(cogxwalkr)
library(table1)
library(haven)
library(stringr)
library(stringi)   
library(purrr)
library(tidyr)

round_2 = function(x) format(round(x, 2), nsmall = 2)

extract_cw_coef_CI = function(path,cogs){
  main = readRDS(path)
  cog.temp = map(path, ~ cogs[stri_detect_fixed(.x, cogs)]) %>% unlist()
  
  str.var = str_remove_all(path, fixed(cog.temp %>% paste0(collapse = '_'))) %>% str_extract_all(., "(?<=_)[^_]+(?=_)") %>% unlist()
  if(length(str.var) == 0) str.var = 'main'
  
  # names in main in inconsistent
  diffs = main[[ if ("cw.diffs" %in% names(main)) "cw.diffs" else "diffs" ]]
  boot  = main[[ if ("cw.boot"  %in% names(main)) "cw.boot"  else "boot"  ]]
  
  cw.coef = lm(paste0(cog.temp[2],'~',cog.temp[1],'-1'),data = diffs)$coef
  cw.coef.se = sd(boot$dist)
  cw.coef.CI = paste0(round_2(cw.coef),' (',round_2(cw.coef-1.96*cw.coef.se),', ',round_2(cw.coef+1.96*cw.coef.se),')')
  
  res = data.frame(pair=paste0(cog.temp[1],'->',cog.temp[2]),str = str.var,coef = cw.coef.CI)
  return(res)
}

generate_path_str = function(cog1,cog2){
  return(c(paste0('boot_',cog1,'_',cog2,'.RDS'), paste0('boot_',c(paste0('age',1:4),paste0('sex',1:2),paste0('race',1:3),paste0('educ',1:4),paste0('marital',1:4),paste0('work',1:2),paste0('income',1:2),paste0('apoe',1:2),paste0('rural',1:2)),'_',cog1,'_',cog2,'.RDS')))
}




setwd('/Users/jwang30/Postdoc/Research/Other/Crosswalk/PsyMCA2')

cogs = c("r1mmse_score", "cog_20", "cog_27", "cog_30", "cog_35", "cog_40")
pairs = t(combn(cogs, 2))  # each row is (A, B)
pair_names  =  apply(pairs, 1, function(x) paste(x[1], x[2], sep = "_"))

################################################################################
# Table 3
cw.result = data.frame()
for(i in seq_len(nrow(pairs))){
  
  cog1 = pairs[i,1]
  cog2 = pairs[i,2]
  
  path.list = generate_path_str(cog1,cog2)
  
  temp = map(path.list, ~ extract_cw_coef_CI(.x, cogs)) %>% do.call(rbind,.)
  cw.result = cw.result %>% bind_rows(temp)
}

write.csv(cw.result,'cw_result.csv',row.names = T)






################################################################################
# Figure 1


df = read_dta('HRS-HCAPcrosswalkdata.dta')
df = df[complete.cases(df[, c("r1mmse_score", "cog_20", "cog_27", "cog_30", "cog_35", "cog_40")]),]


get_pair_df <- function(c1, c2) {
  nm     <- paste(c1, c2, sep = "_")
  nm_rev <- paste(c2, c1, sep = "_")
  p1     <- sprintf("boot_%s.RDS", nm)
  p2     <- sprintf("boot_%s.RDS", nm_rev)
  
  if (file.exists(p1)) {
    res  <- readRDS(p1)
    pair <- c(c1, c2)
  } else if (file.exists(p2)) {
    res  <- readRDS(p2)
    pair <- c(c2, c1)
  } else {
    stop("Neither file found: ", p1, " or ", p2)
  }
  
  diffs_name <- if ("cw.diffs" %in% names(res)) "cw.diffs" else "diffs"
  cand <- res[[diffs_name]]
  
  # ensure data.frame and drop 'iteration' if present
  if (!is.data.frame(cand)) cand <- tibble::as_tibble(cand)
  if ("iteration" %in% names(cand)) cand$iteration <- NULL
  
  # require the expected columns and order them according to 'pair'
  missing <- setdiff(pair, names(cand))
  if (length(missing)) stop("Missing expected columns in diffs: ", paste(missing, collapse = ", "))
  
  cand[, ..pair, drop = FALSE]
}

get_pair_se <- function(c1, c2) {
  nm     <- paste(c1, c2, sep = "_")
  nm_rev <- paste(c2, c1, sep = "_")
  p1     <- sprintf("boot_%s.RDS", nm)
  p2     <- sprintf("boot_%s.RDS", nm_rev)
  
  if (file.exists(p1)) {
    res  <- readRDS(p1)
    pair <- c(c1, c2)
  } else if (file.exists(p2)) {
    res  <- readRDS(p2)
    pair <- c(c2, c1)
  } else {
    stop("Neither file found: ", p1, " or ", p2)
  }
  
  diffs_name <- if ("cw.boot" %in% names(res)) "cw.boot" else "boot"
  sd(res[[diffs_name]]$dist)
}

var.names.slides = cogs %>% rev
var.label.slides = c("MMSE","TICS-20","TICS-27","TICS-30","TICS-35","TICS-40") %>% rev

df6 = df %>% select(all_of(var.names.slides))


# LOWER: scatter + no-intercept fit line
lower_fun <- function(data, mapping, ...) {
  xvar <- rlang::as_name(mapping$x)
  yvar <- rlang::as_name(mapping$y)
  d <- get_pair_df(xvar, yvar)
  
  ggplot(d, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_point(alpha = 0.3, size = 0.5) +
    # geom_smooth(method = "lm", formula = y ~ x + 0, se = FALSE) +
    theme_minimal(base_size = 8) +
    theme(
      axis.title = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(2, 2, 2, 2)
    )
}

# UPPER: print slope from lm(y ~ x + 0)
upper_fun <- function(data, mapping, ...) {
  xvar <- rlang::as_name(mapping$x)
  yvar <- rlang::as_name(mapping$y)
  d <- get_pair_df(xvar, yvar)
  se = get_pair_se(xvar, yvar)
  
  slope <- 0
  if (nrow(d)) {
    slope <- coef(lm(as.formula(paste0(yvar, " ~ ", xvar, " + 0")), data = d))[1]
  }
  
  xlabel = var.label.slides[which(xvar == var.names.slides)]
  ylabel = var.label.slides[which(yvar == var.names.slides)]
  
  lbl = paste0(xlabel,' -> ',ylabel,'\n',round_2(slope),'\n(',round_2(slope-1.96*se),', ',round_2(slope+1.96*se),')')

  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = lbl, size = 4) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void()  
}

diag_fun <- function(data, mapping, ...) {
  # lb_ub = quantile(data[,1:6] %>% unlist %>% as.numeric,c(.1,.9))
  ggplot(data = data, mapping = mapping) +
    geom_density(aes(y=..scaled..),adjust=3) + 
    # xlim(lb_ub) +
    theme(text = element_text(size=25),
          axis.title = element_text(size = 20),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size = 20)) 
}





p.slides = GGally::ggpairs(
  df6,
  columns = seq_along(var.names.slides),   
  upper = list(continuous = wrap(upper_fun)),
  lower = list(continuous = wrap(lower_fun)),
  diag  = list(continuous = wrap(diag_fun)),
  columnLabels = var.label.slides,
  switch = 'y'
) +
  theme(
    strip.text = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
    axis.text.y = element_text(size = 10)
  )



gt <- GGally::ggmatrix_gtable(p.slides)


# find top and left strips
lay <- gt$layout
top_ids  <- which(grepl("strip-t", lay$name))   # column labels (top)
left_ids <- which(grepl("strip-l", lay$name))   # row labels (left)

# order them in natural left-to-right / top-to-bottom order
top_ids  <- top_ids[order(lay$l[top_ids])]
left_ids <- left_ids[order(lay$t[left_ids])]

# helper: recursively set fill on any rect grobs inside a grob
set_strip_colors <- function(grob, fill, text_color) {
  if (inherits(grob, "rect")) {
    grob$gp$fill <- fill
  }
  if (inherits(grob, "text")) {
    grob$gp$col <- text_color
  }
  if (!is.null(grob$children)) {
    for (i in seq_along(grob$children)) {
      grob$children[[i]] <- set_strip_colors(grob$children[[i]], fill, text_color)
    }
  }
  if (!is.null(grob$grobs)) {
    for (i in seq_along(grob$grobs)) {
      grob$grobs[[i]] <- set_strip_colors(grob$grobs[[i]], fill, text_color)
    }
  }
  grob
}

# color logic: background + matching text color
col_for_index <- function(i) {
  if (i <= 5) list(fill = "#FBE6D4", text = "black")
  else if (i > 5) list(fill = "#CCE3E1", text = "black")
  else list(fill = "grey90", text = "black")
}

# recolor top strips (columns)
for (i in seq_along(top_ids)) {
  idx <- top_ids[i]
  colors <- col_for_index(i)
  gt$grobs[[idx]] <- set_strip_colors(gt$grobs[[idx]], colors$fill, colors$text)
}

# recolor left strips (rows)
# left_ids <- left_ids[order(lay$b[left_ids])]  # sort by bottom position, not top

# Then recolor row (left) strips correctly
for (i in seq_along(left_ids)) {
  idx <- left_ids[i]
  colors <- col_for_index(i)  # same color mapping: 1–5 pink, last 2 blue
  gt$grobs[[idx]] <- set_strip_colors(gt$grobs[[idx]], colors$fill, colors$text)
}


pdf('/Users/jwang30/Postdoc/Research/Other/Crosswalk/PsyMCA2/p2f1.pdf',width = 10,height = 10)
grid::grid.newpage()
grid::grid.draw(gt)
dev.off()






















