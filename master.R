rm(list = ls())
SPATH = paste(getwd(),"/Data_and_Code_GPMMCM",sep="")

# Re-produce Figure 1
source(paste(SPATH, '/Code/fig1.R', sep=''))

# Re-produce Figure 2
source(paste(SPATH, '/Code/fig2.R', sep=''))

# Re-produce Figure 3. Figure 4
source(paste(SPATH, '/Code/fig3.4.R', sep=''))

# Re-produce Figure 5. Figure 6
source(paste(SPATH, '/Code/fig5.6.R', sep=''))

# Re-produce Figure 7
source(paste(SPATH, '/Code/fig7.R', sep=''))

# Re-produce Figure 8
source(paste(SPATH, '/Code/fig8.R', sep=''))


# Re-produce Table 2,Table 3, Table 4
source(paste(SPATH, '/Code/table2.3.4.R', sep=''))

# Re-produce Table 5,  Table 6
source(paste(SPATH, '/Code/table5.6.R', sep=''))

# Re-produce Table 7 
source(paste(SPATH, '/Code/table7.R', sep=''))

# Re-produce faithful_code.Rdata
source(paste(SPATH, '/Code/faithful_code.R', sep=''))

# Re-produce faithful_code_01.Rdata
source(paste(SPATH, '/Code/faithful_code_01.R', sep=''))

# Re-produce vdeq_code.Rdata
source(paste(SPATH, '/Code/vdeq_code.R', sep=''))

# Re-produce hawks_code.Rdata
source(paste(SPATH, '/Code/hawks_code.R', sep=''))

# Re-produce sim_VVE_code.R
source(paste(SPATH, '/Code/sim_VVE_code.R', sep=''))
