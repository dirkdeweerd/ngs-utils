ChromosomalFractionPerSample = function(bins)
{
  chr.frac = NULL
  if (length(bins) == 2)
  {
    bins.fwd <- bins[[1]]
    bins.rev <- bins[[2]]
    chr.frac.fwd = rowSums(bins.fwd) / sum(bins.fwd) / 2
    chr.frac.rev = rowSums(bins.rev) / sum(bins.rev) / 2
    chr.frac = rbind(chr.frac, c(chr.frac.fwd, chr.frac.rev))
  }
  else
  {
    bins.all <- bins[[1]]
    chr.frac.all = rowSums(bins.all) / sum(bins.all)
    chr.frac = rbind(chr.frac, chr.frac.all)
  }
  return(chr.frac)
}