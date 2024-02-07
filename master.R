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

# Re-produce Figure 6
source(paste(SPATH, '/Code/fig6.R', sep=''))

# Re-produce Table 2,  Table 3
source(paste(SPATH, '/Code/table2.3.R', sep=''))

# Re-produce Table 4 
source(paste(SPATH, '/Code/table4.R', sep=''))

# Re-produce faithful_code.RData
source(paste(SPATH, '/Code/faithful_code.R', sep=''))

# Re-produce faithful_code_01.RData
source(paste(SPATH, '/Code/faithful_code_01.R', sep=''))

# Re-produce hawks_code.RData
source(paste(SPATH, '/Code/hawks_code.R', sep=''))

# Re-produce sim_VVE_code.RData
source(paste(SPATH, '/Code/sim_VVE_code.R', sep=''))
