*********************************************************************************************************
* 		Important file paths
*********************************************************************************************************
	*Log_files
	global logs_cs /home/log_files................
	*File path data bases CS:
	global bases_cs /home/data_bases..................
	*File path results
	global resultados_cs /home/results............
	
************************************************************************************************************
** Cost-Sharing in Medical Care Can Increase Adult Mortality: Evidence from Colombia
** Aim: impact analyses using a Dinamyc Regression Discontinuity (Cellini-2010) - phase 1 of the analysis

//  task:      	ITT-LLR analyses for the number of outpatient services(OUT)
//  authors:    Giancarlo Buitrago, Javier Amaya-Nieto, Marcos Vera-Hernandez and Grant Miller 
//  date: 		February 20 / 2024
************************************************************************************************************
	
	
	log using "${logs_cs}/log_ITT_LLR_services_OUT.smcl", replace
		
	**** Loading the data(not publicly available on request to the Ministry of Health of Colombia)
	
	use PersonaBasicaID mmw mes OUT using "${bases_cs}/Panel_2011_2018_full_OUT.dta", clear
	
		** Variables meaning
			* mmw: monthly minimum wage
			* mes: month of the observation registered
			* Out: Number of the outpatient services per month for each register

	*******************************************************************************************************	
	*******************		 Centrando el mmw

		gen mmwC=mmw-5

	*******************************************************************************************************	
	*******************		 Setting the dataset as time series

		sort PersonaBasicaID mes 
		tsset PersonaBasicaID  mes
		
** Setting a matrix to collect results from all the lags about ITT results
		matrix index=J(96,10,-99)
		matrix colnames index = ITT CI_low CI_high p Ori_obs Lag Efec_obs Obtimal_BW Mean_left Pro_change
		matrix list index
	
	** Seting the row to collect data
		local irow=0
		
	** Setting a loop to run the regression discontinuity (using the package rdrobust[https://rdpackages.github.io/rdrobust/]) for each lag in the 96 months period. 
	
		foreach lag of numlist 95(1)0 { 

			** generating the lagged mmw to feed the rdrobust. 
				by PersonaBasicaID : gen  mmwC_lag`lag'= l`lag'.mmwC

			**************************************** RD ANALYSIS

				display "*******************  OUT Lag = `lag' **********************"

				rdrobust OUT mmwC_lag`lag', c(0) vce(cluster PersonaBasicaID) masspoints(off)  

		
			** Collecting results
				local ++irow
				
				matrix index[`irow',1]=e(tau_cl)
				matrix index[`irow',2]=e(ci_l_cl)
				matrix index[`irow',3]=e(ci_r_cl)
				matrix index[`irow',4]=e(pv_cl)
				matrix index[`irow',5]=e(N)	
				matrix index[`irow',6]=`lag'
				matrix index[`irow',7]=`e(N_h_l)'+`e(N_h_r)'
				matrix index[`irow',8]=e(h_l)
				matrix index[`irow',9]=e(tau_cl_l)
				matrix index[`irow',10]=(e(tau_cl)/e(tau_cl_l))*100	
											
		** Dropping the lagged mmw variable created at the beggining of the loop.  
			drop mmwC_lag`lag'			
		
		
		** Exporting results 
		putexcel set "${resultados_cs}/ITT_LLR_OUT", replace
		putexcel A1=matrix(index), names		
		
		}

*******************************************************************************************************	
log close

