library(vcfR)
library(rehh)


# Read in arguments
args = commandArgs(trailingOnly=TRUE)

pop1 = args[1]
pop2 = args[2]
out = args[3]

hh1 = data2haplohh(pop1)
hh2 = data2haplohh(pop2)

ss1 = scan_hh(hh1) #Calc iHH and iES
ss2 = scan_hh(hh2)

ihs = ihh2ihs(ss1)
xp.ehh = ies2xpehh(ss1, ss2)
colnames(ihs)[4] = "ihs.p"
colnames(xp.ehh)[4] = "xpehh.p"

HAP = merge(ihs, xp.ehh, by=c("CHR","POSITION"))

write.table(HAP, out, quote=F)

#Write out the R objects...specifically the hh1 object