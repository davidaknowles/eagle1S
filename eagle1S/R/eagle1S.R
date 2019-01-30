
logit=function(p) {log(p/(1-p))}
inv_logit=function(g) {1/(1+exp(-g))}

phased_types=c("0|0","0|1","1|0","1|1")
phased_hets=c("0|1","1|0")

get_lead_helper=function(results_table) {
  results_table %>%
  mutate(p=pchisq(lrt, df = df, lower.tail = F) ) %>%
  group_by(exonic_snp_pos) %>%
  mutate( bonf_p = p * length(p) ) %>%
  top_n(1, -bonf_p) %>%
  ungroup() %>%
  mutate(q=p.adjust(pmin(1, bonf_p), method="BH")) %>%
  select(reg_snp_pos, exonic_snp_pos, p, bonf_p, q, conc, starts_with("beta_"))
}

#' Get lead eQTLs
#'
#' @param results_table from `eagle1S`
#' @return Tibble with per exonic SNP top/lead eSNPs
#' @import dplyr
#' @export
get_lead_eQTL = function(results_table) {
  results_table %>% mutate(lrt = 2.0*(l_geno - l0)) %>% get_lead_helper()
}

#' Get lead gxeQTLs
#'
#' @param results_table from `eagle1S`
#' @return Tibble with per exonic SNP top/lead gxeSNPs
#' @import dplyr
#' @export
get_lead_gxeQTL = function(results_table) {
  results_table %>% mutate(lrt = 2.0*(l_interact - l_geno)) %>% get_lead_helper()
}

#' Map eQTLs and gxeQTLs using allele specific expression
#'
#' @param ase A table with SNP position ("POS"), alternative allele ("ALT"), "individual", reference counts "r" and alternative allele counts "a"
#' @param phased A VCF table, i.e. genotypes should start at column 11. Should have phased genotypes e.g. 0|1. Missing values, and equivalent unphased values, are allowed and will be ignored. Column names must correspond to the individuals and be compatible with `ase` and `meta` tables.
#' @param meta A table with an "individual" column and environmental factor "x". x should be normalized to some sensible scaling.
#' @return Tibble with log likelihoods for all tested SNP pairs.
#' @import dplyr
#' @import magrittr
#' @import foreach
#' @importFrom rstan optimizing
#' @export
eagle1S = function(ase,
                   phased,
                   meta,
                   checkpoint_dir = NULL,
                   cisdist=100000,
                   min_samples=10,
                   concShape=1.001,
                   concRate=0.001,
                   min_reads=1000,
                   max_iter = 30) {

  phased %<>% unite( pos_alt, POS, ALT, sep="_", remove=F) %>%
    distinct(pos_alt, .keep_all = TRUE)
  class(phased)="data.frame"

  rownames( phased )=phased$pos_alt

  ase %<>%  mutate( pos_alt=paste(POS, ALT, sep="_"),
                    individual=as.character(individual)) %>%
    filter( pos_alt %in% rownames(phased) ) %>%
    mutate( geno=phased[ cbind(as.character(pos_alt), individual) ] )

  exonic_snps = unique(ase$POS)

  foreach(exonic_snp_pos=exonic_snps,
                 .errorhandling = if (interactive()) "stop" else "remove",
                 .combine = bind_rows) %dopar% {

     # checkpointing per "gene"
     if (!is.null(checkpoint_dir)){
       check_fn=paste0(checkpoint_dir,exonic_snp_pos,".txt.gz")
       if (file.exists(check_fn) & !interactive()) {
         return(read_tsv(check_fn) %>%
                  mutate(exonic_snp_pos=exonic_snp_pos))
       }
     }

     gene_ase = ase %>% filter( POS == exonic_snp_pos )
     if (nrow(gene_ase) < min_samples) return(NULL)

     allelic_count_total=sum(gene_ase$a+gene_ase$r)
     cat("Allelic total count ", allelic_count_total, "\n")
     if (allelic_count_total < min_reads) return(NULL)

     cis_snps=phased %>%
       filter( (exonic_snp_pos-cisdist) < POS,
               (exonic_snp_pos+cisdist) > POS ) %>% .$POS

     # iterate over cis SNPs
     temp_results = foreach(snp_pos=cis_snps,
                            .errorhandling = if (interactive()) "stop" else "remove",
                            .combine = bind_rows) %do% {

       reg_geno = (phased %>% filter( POS == snp_pos ))[,11:ncol(phased)] %>% as.matrix()

       if (nrow(reg_geno) != 1) {
         print("Skipping >biallelic site")
         return(NULL)
       }

       reg_geno=data_frame(individual=colnames(reg_geno),
                           reg_geno=as.character(reg_geno))

       # join the ASE and cisSNP phased genotypes
       ase_temp = gene_ase %>% inner_join(reg_geno, by="individual") %>%
         filter( geno %in% phased_hets, # signal only comes from het exonic SNPs
                 reg_geno %in% phased_types, # require phased cisSNP
                 (r+a) > 0 ) %>% # must have some coverage
         mutate( het_x=ifelse(reg_geno %in% phased_hets, # if cisSNP is het
                              ifelse(geno == reg_geno,1,-1), # is it in phase with the exonicSNP?
                              0) )

       if (nrow(ase_temp) < min_samples) return(NULL)
       num_het_snps=sum(ase_temp$het_x != 0)
       if (num_het_snps < min_samples) return(NULL) # no heterozygous regulatory SNPs

       ase_temp %<>% left_join(meta, by="individual")

       # LM for initialization of GLM
       coverage =  with(ase_temp, a+r)
       y = logit( ase_temp$a/coverage )

       get_init=function(x_mat) {
         beta_init = solve( t(x_mat) %*% x_mat, t(x_mat) %*% y ) # LM
         # now try to get a sensible initialization for the concentration param
         # BB variance is n*p*(1-p)*(conc+n)/(conc+1).
         # np(1-p) is binomial variance.
         # Second term: (conc+1+(n-1))/(conc+1)=1+(n-1)/(conc+1).
         prob = inv_logit( x_mat %*% beta_init )
         fitted_a = coverage * prob
         var_a = (ase_temp$a - fitted_a)^2
         bin_var = coverage * prob * (1-prob)
         #(coverage - 1) / (var_a / bin_var - 1) - 1 # method of moments
         conc_init = mean((coverage + 1) / (var_a / bin_var + 1) - 1) # like adding Gamma(2,2) pseudocounts
         if (conc_init < 1.5 | conc_init > 100.) conc_init=10.
         list(beta=as.numeric(beta_init) %>% as.array(),
              conc=conc_init)
       }

       # are we testing the exonic SNP itself (or a perfectly linked SNP)
       testing_self_flag = length(unique(ase_temp$het_x)) == 1

       # full model with eqtl and gxe
       x_full = if (!testing_self_flag)
         model.matrix( ~ het_x + x:het_x, data=ase_temp ) else
           model.matrix( ~ x, data=ase_temp )

       stan_dat = list(N=nrow(x_full),
                       P=ncol(x_full),
                       x=x_full,
                       ys=ase_temp$a,
                       ns=ase_temp$a + ase_temp$r,
                       concShape=concShape,
                       concRate=concRate)

       fit_full=if (det(t(x_full) %*% x_full) > 0) {
         optimizing(stanmodels$bb,
                  data=stan_dat,
                  init=get_init(x_full),
                  iter=max_iter)
       } else {
         list(value=NA,
              par=rep(NA, 1 + ncol(x_full)))
       }
       # TODO: does this work when testing the exonic SNP itself?
       names(fit_full$par) = c("conc","intercept", if (!testing_self_flag) "b_eqtl", "b_gxe")

       # eqtl but no gxe model
       x_eqtl = if (!testing_self_flag)
         model.matrix( ~ het_x , data=ase_temp ) else
           model.matrix( ~ 1, data=ase_temp )
       stan_dat$x=x_eqtl
       stan_dat$N=nrow(x_eqtl)
       stan_dat$P=ncol(x_eqtl)
       fit_eqtl=optimizing(stanmodels$bb,
                           data=stan_dat,
                           init=get_init(x_eqtl),
                           iter=max_iter)$value

       # null model, no eQTL or GxE
       x_0 = model.matrix( ~ 1, data=ase_temp )
       stan_dat$x=x_0
       stan_dat$N=nrow(x_0)
       stan_dat$P=ncol(x_0)
       fit_0=optimizing(stanmodels$bb,
                        data=stan_dat,
                        init=get_init(x_0),
                        iter=max_iter)$value

       df=ncol(x_full) - ncol(x_eqtl)

       if (F) { # debugging code
         ase_temp %>% mutate(het=reg_geno,
                             coverage=r+a,
                             ar = a/coverage,
                             in_phase=geno == reg_geno) %>%
           ggplot(aes( x,ar,size=coverage,col=factor(het_x))) + geom_point() + ylim(0,1) +
           xlab("Environmental factor") + ylab("Phased allelic ratio")
         pchisq( 2.0*(fit_full$value - fit_eqtl), df = df, lower.tail = F)
         pchisq( 2.0*(fit_eqtl - fit_0), df = 1, lower.tail = F)
       }

       num_het_ind=length(unique(ase_temp %>% filter(het_x != 0) %>% .$individual)) # for performance analysis later

       data_frame(total_count=sum(coverage),
                  num_het_snps=num_het_snps,
                  num_het_ind=num_het_ind,
                  reg_snp_pos=snp_pos,
                  df=df,
                  l0=fit_0,
                  l_geno=fit_eqtl,
                  l_interact=fit_full$value,
                  ) %>% cbind(as_data_frame(as.list(fit_full$par)))
     }

     if (!is.null(checkpoint_dir)){
       print("Saving results")
       checkpoint_file= gzfile( check_fn,"w")
       temp_results %>% write_tsv(checkpoint_file) # write_tsv
       close(checkpoint_file)
     }

     temp_results %>%
       mutate(exonic_snp_pos=exonic_snp_pos)
    }
}
