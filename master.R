rm(list = ls())
SPATH = paste(getwd(),"/Data_and_Code_GPMMCM",sep="")

# Re-produce Figure 1
source(paste(SPATH, '/Code/fig1.R', sep=''))

# Re-produce Figure 2
source(paste(SPATH, '/Code/fig2.R', sep=''))

# Re-produce Figure 3. Figure 4
source(paste(SPATH, '/Code/fig3.4.R', sep=''))

# Re-produce Figure 5
source(paste(SPATH, '/Code/fig5.R', sep=''))

# Re-produce Figure S1
source(paste(SPATH, '/Code/Fig.S1.R', sep=''))

# Re-produce Figure S2(a)
source(paste(SPATH, '/Code/Fig.S2(a).R', sep=''))

# Re-produce Figure S2(b)
source(paste(SPATH, '/Code/Fig.S2(b).R', sep=''))

# Re-produce Table 1,  Table 2
source(paste(SPATH, '/Code/table1.2.R', sep=''))

# Re-produce Table 3
source(paste(SPATH, '/Code/table3.R', sep=''))

# Re-produce faithful_code.RData
source(paste(SPATH, '/Code/faithful_code.R', sep=''))

# Re-produce faithful_code_01.RData
source(paste(SPATH, '/Code/faithful_code_01.R', sep=''))

# Re-produce hawks_new.RData
source(paste(SPATH, '/Code/fit_hawkdata.R', sep=''))

# Re-produce fit_DMS
source(paste(SPATH, '/Code/fit_DMSdata.R', sep=''))
