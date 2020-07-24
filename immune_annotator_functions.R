make_immune_key <-  function() {
  # Primary patterns
  tcell <- c("t-cell", "tcell", " t cell", "^t cell", "t-lymphocyte", "t lymphocyte",
             "^cd4$", " th[12]", "^th[12]", "cd4+ ", "cd4t", "cd8t", "helper", "cd4\\+",
             "cd8[ _psa\\+]", "cd8$", "cd3\\+", "cd3 p", "treg", "regulatory[ _]t",
             "regulatory cell", "teff", "thymocyte", "mait", "cd4-[pt]",
             "[( ^]th0[ )]", "^th0$", "temra", "tfh") %>% paste(.,collapse = "|")
  bcell <- c("b-cell", "bcell", " b cell", "^b cell", "b lympho", "b memory", 
             "plasmablast", "plasma cell", "cd19\\+", "cd19 ", 
             "cd19-b") %>% paste(.,collapse = "|")
  ilc <- c("innate lymph", "^ilc", " ilc") %>% paste(.,collapse = "|")
  nk <- c("nk cell", "nkcell", "natural killer", "cd56[db\\+]") %>% paste(.,collapse = "|")
  dc <- c("dendritic.cell", "-cdc", " cdc", "plasmacytoid", "^pdc", 
          "[ _(]pdc", "myeloid dc", "+ dc", "predc", "langerhans cells", 
          "modc", "cd11c\\+", "dentritic") %>% paste(.,collapse = "|")
  macs <- c("cd11b[ \\+]", "macrop") %>% paste(.,collapse = "|")
  mono <- c("monocyte") %>% paste(.,collapse = "|")
  pmn <- c("^pmn", "neutroph") %>% paste(.,collapse = "|")
  baso <- c("basophils", "^basophil$") %>% paste(.,collapse = "|")
  eos <- c("eosinophil ", "eosinophils", "eosinophil$") %>% paste(.,collapse = "|")
  mast <- c("mast cell")
  
  #Key assembly
  immune_key <- as_tibble(matrix(c(
    "t-lymphocyte", tcell,
    "b-lymphocyte", bcell,
    "innate lymphoid cell", ilc,
    "natural killer cell", nk,
    "dendritic cell", dc,
    "macrophage", macs,
    "monocyte", mono,
    "neutrophil", pmn,
    "basophil", baso,
    "eosinophil", eos,
    "mast cell", mast
  ), ncol = 2, byrow = T, dimnames = list(NULL ,c("id", "pattern"))))
  
  immune_key
}


pattern_match <- function(input, key){
  key %>%
    mutate(match = map2_chr(pattern, id, 
                            function(pattern, id) ifelse(grepl(pattern, input), id, NA))) %>%
    filter(!is.na(match)) %>%
    pull(match) %>%
    sort() %>% 
    paste(., collapse = "_") %>%
    ifelse(.=="", NA, .)
}

make_immune_subset_key <- function(){
  #T-cell subtypes
  t_treg <- c("treg", "regulatory[ _]t", "t regu") %>% paste(.,collapse = "|")
  t_cd4 <- c("cd4[ \\+t_]", "cd4$", "cd4\\-[pt]", "help") %>% paste(.,collapse = "|")
  t_cd8 <- c("cd8[^n\\-]", "cd8$", "cd8\\-[pt]", "cytotoxic") %>% paste(.,collapse = "|")
  t_tfh <- c("follicular", "tfh") %>% paste(.,collapse = "|")
  th0 <- c("th0", "0 cell") %>% paste(.,collapse = "|")
  th1 <- c("th1$", "th1[a-z _\\-]", "1 cell") %>% paste(.,collapse = "|")
  th2 <- c("th2", "2 cell") %>% paste(.,collapse = "|")
  th17 <- c("th17", "17 cell")  %>% paste(.,collapse = "|")
  t_eff <- c("effector", "tem$", "tem ", "teff") %>% paste(.,collapse = "|")
  t_central <- c("central", "tcm") %>% paste(.,collapse = "|")
  t_mem <- c("memory",  "tem$", "tem ", "tcm", "cd45ro[ \\+]") %>% paste(.,collapse = "|")
  t_naive <- c("na[iï]ve", "cd45ra[ \\+]") %>% paste(.,collapse = "|")
  t_mait <- c("mait", "invariant") %>% paste(.,collapse = "|")
  t_inkt <- c("inkt")
  t_gd <- c("gamma[ \\-]*delta($| )", "gd t", "γδ") %>% paste(.,collapse = "|")
  t_temra <- c("temra")
  t_prog <- c("prec","prog", "cd34[ \\+]") %>% paste(.,collapse = "|")
  
  tcell_key <- as_tibble(matrix(c(
    "cd4+ regulatory", t_treg,
    "cd4+", t_cd4,
    "cd8+", t_cd8,
    "follicular", t_tfh,
    "th0 cd4+", th0,
    "th1 cd4+", th1,
    "th2 cd4+", th2,
    "th17 cd4+", th17,
    "effector", t_eff,
    "central", t_central,
    "memory", t_mem,
    "naive", t_naive,
    "mait", t_mait,
    "inkt", t_inkt,
    "gamma-delta", t_gd,
    "temra memory", t_temra,
    "progenitor", t_prog
  ), ncol = 2, byrow = T, dimnames = list(NULL ,c("id", "pattern")))) %>% 
    mutate(id = paste0(id, " t-lymphocyte")) %>% 
    mutate(group = "t-lymphocyte") 
  
  #B-cell subtypes
  b_mem <- c("memory",  "cd27\\+") %>% paste(.,collapse = "|")
  b_plasma <- c("plasma ")
  b_plasmablast <- c("plasmablast")
  b_naive <- c("na[iï]ve", "nave") %>% paste(.,collapse = "|")
  b_germ <- c("germ")
  b_fol <- c("follic")
  b_marg <- c("marginal")
  b_igg <- c("igg[ \\+]")
  b_igd <- c("igd[ \\+]")
  b_igm <- c("igm[h \\+]")
  b_iga <- c("iga[ \\+]")
  b_prog <- c("prec","prog", "cd34\\+", "cd34 (^n|p)") %>% paste(.,collapse = "|")
  
  bcell_key <- as_tibble(matrix(c(
    "naive", b_naive,
    "memory", b_mem,
    "plasma cell", b_plasma,
    "plasmablast", b_plasmablast,
    "germinal center", b_germ,
    "follicular", b_fol,
    "marginal zone", b_marg,
    "igg+", b_igg,
    "igd+", b_igd,
    "igm+", b_igm,
    "iga+", b_iga,
    "progenitor", b_prog
  ), ncol = 2, byrow = T, dimnames = list(NULL ,c("id", "pattern")))) %>% 
    mutate(id = paste0(id, " b-lymphocyte"))%>% 
    mutate(group = "b-lymphocyte") 
  
  #ILC subtypes
  ilc_1 <- c("group 1", "ilc1", " 1$") %>% paste(.,collapse = "|")
  ilc_2 <- c("group 2", "ilc2") %>% paste(.,collapse = "|")
  ilc_3 <- c("group 2", "ilc3", " 3$") %>% paste(.,collapse = "|")
  
  ilc_key <- as_tibble(matrix(c(
    "group 1", ilc_1,
    "group 2", ilc_2,
    "group 3", ilc_3
  ), ncol = 2, byrow = T, dimnames = list(NULL ,c("id", "pattern")))) %>% 
    mutate(id = paste0(id, " innate lymphoid cell"))%>% 
    mutate(group = "innate lymphoid cell") 
  
  #DC subtypes
  dc_pdc <- c("plasmacytoid", "pdc") %>% paste(.,collapse = "|")
  dc_prog <- c("pre[dc-]","prog") %>% paste(.,collapse = "|")
  dc_mono <- c("monocyte", "modc") %>% paste(.,collapse = "|")
  dc_cdc <- c("conv", "cdc", "classic") %>% paste(.,collapse = "|")
  dc_cd141 <- c("cd141")
  dc_cd1c <- c("cd1c")
  dc_bdca1 <- c("bdca1")
  dc_lang <- c("langerhans cells", "lang", "epiderm") %>% paste(.,collapse = "|")
  dc_my <- c("myeloid")
  
  dc_key <- as_tibble(matrix(c(
    "plasmacytoid", dc_pdc,
    "progenitor", dc_prog,
    "monocyte-derived", dc_mono,
    "langerhans", dc_lang,
    "conventional", dc_cdc,
    "cd141+ conventional", dc_cd141,
    "cd1c+ conventional", dc_cd1c,
    "bdca1+ conventional", dc_bdca1,
    "myeloid", dc_my
  ), ncol = 2, byrow = T, dimnames = list(NULL ,c("id", "pattern")))) %>% 
    mutate(id = paste0(id, " dendritic cell")) %>% 
    mutate(group = "dendritic cell") 
  
  #NK subtypes
  nk_56b <- c("cd56high", "cd56\\+", "cd56bright") %>% paste(.,collapse = "|")
  nk_56d <- c("cd56dim", "cd56\\-") %>% paste(.,collapse = "|")
  nk_16p <- c("cd16\\+") %>% paste(.,collapse = "|")
  nk_16n <- c("cd16\\-") %>% paste(.,collapse = "|")
  nk_inv <- c("invariant")
  
  nk_key <- as_tibble(matrix(c(
    "cd56bright", nk_56b,
    "cd56dim", nk_56d,
    "cd16+", nk_16p,
    "cd16-", nk_16n,
    "invarient", nk_inv
  ), ncol = 2, byrow = T, dimnames = list(NULL ,c("id", "pattern")))) %>% 
    mutate(id = paste0(id, " natural killer cell"))%>% 
    mutate(group = "natural killer cell") 
  
  #monocyte subtypes
  mono_classic <- c("^class", " class") %>% paste(.,collapse = "|")
  mono_nonclassic <- c("nonclass", "non class", "non\\-class") %>% paste(.,collapse = "|")
  mono_cd14p <- c("cd14\\+", "cd14[ p]") %>% paste(.,collapse = "|")
  mono_cd14d <- c("cd14\\-[c ]", "cd14dim") %>% paste(.,collapse = "|")
  mono_cd16p <- c("cd16\\+", "cd16\\-p", "cd16p") %>% paste(.,collapse = "|")
  mono_cd16n <- c("cd16\\-[ n]", "cd16n") %>% paste(.,collapse = "|")
  
  mono_key <- as_tibble(matrix(c(
    "classical", mono_classic,
    "non-classical", mono_nonclassic,
    "cd14+", mono_cd14p,
    "cd14dim", mono_cd14d,
    "cd16+", mono_cd16p,
    "cd16-", mono_cd16n
  ), ncol = 2, byrow = T, dimnames = list(NULL ,c("id", "pattern")))) %>% 
    mutate(id = paste0(id, " monocyte"))%>% 
    mutate(group = "monocyte")
  
  #macrophage subtypes
  mac_m1 <- c("m1", "classic") %>% paste(.,collapse = "|")
  mac_m2 <- c("m2", "alternative") %>% paste(.,collapse = "|")
  mac_m0 <- c("m0") %>% paste(.,collapse = "|")
  mac_prog <- c("prec","prog") %>% paste(.,collapse = "|")
  
  mac_key <- as_tibble(matrix(c(
    "m1", mac_m1,
    "m2", mac_m2,
    "m0", mac_m0,
    "progenitor", mac_prog
  ), ncol = 2, byrow = T, dimnames = list(NULL ,c("id", "pattern")))) %>% 
    mutate(id = paste0(id, " macrophage"))%>% 
    mutate(group = "macrophage") 
  
  #Compile
  subset_key <- rbind(tcell_key, bcell_key, ilc_key, dc_key, mono_key, nk_key, mac_key)
  subset_key
}

make_immune_name_order <- function(){
  t_order <- list(
    c("progenitor",  
      "temra", "naive", "effector", "central", "memory", 
      "th0", "th1", "th2", "th17", "follicular", "regulatory", 
      "cd4+", "cd8+", 
      "mait", "inkt", "gamma-delta",
      "t-lymphocyte"))
  
  b_order <- list(
    c("progenitor",  
      "igd+", "igm+", "igg+", "iga+",
      "naive", "memory", 
      "follicular", "marginal", "zone", "germinal", "center",
      "plasma", "cell", "plasmablast", 
      "b-lymphocyte"))
  
  ilc_order <- list(c("group", "1", "2", "3", "innate", "lymphoid", "cell"))
  
  dc_order <- list(
    c("progenitor",  
      "cd1c+", "cd141+", "bdca1+", 
      "plasmacytoid", "monocyte-derived", "langerhans", "conventional", "myeloid",
      "dendritic", "cell"))
  
  nk_order <- list(
    c("cd56bright", "cd56dim", "cd16+", "cd16-",
      "invariant",
      "natural", "killer", "cell"))
  
  mono_order <- list(
    c("cd14dim", "cd14+", "cd16+", "cd16-",
      "classical", "non-classical",
      "monocyte"))
  
  mac_order <- list(
    c("progenitor",  
      "m0", "m1", "m2",
      "macrophage"))
  
  order_key <- tibble(
    "t-lymphocyte" = t_order,
    "b-lymphocyte" = b_order,
    "innate lymphoid cell" = ilc_order,
    "dendritic cell" = dc_order,
    "natural killer cell" = nk_order,
    "monocyte" = mono_order,
    "macrophage" = mac_order,
    index = 1) %>% 
    pivot_longer(-index, names_to = "id", values_to = "order") %>% 
    select(-index)
  order_key
}



pattern_match_3 <- function(input_value, group_value, pattern_key, order_key){
  
  order <- order_key %>% 
    filter(id == group_value) %>% 
    pull(order) %>% unlist()
  
  pattern_match <- pattern_key %>%
    mutate(match = map2_chr(pattern, id, 
                            function(pattern, id) ifelse(grepl(pattern, input_value), id, NA))) %>%
    filter(!is.na(match) & group == group_value) %>%
    pull(match) 
  
  pattern_match <- unique(
    unlist(str_split(pattern_match, " "))
  ) 
  
  sort(factor(pattern_match, levels = order, ordered=TRUE)) %>%
    paste(., collapse = " ")
}  