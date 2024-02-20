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
** Aim: impact analyses using a Dinamyc Regression Discontinuity (Cellini-2010) - phase 2 of the analysis

//  task:      	ITT-LLR analyses for the number of outpatient services(OUT)
//  authors:    Giancarlo Buitrago, Javier Amaya-Nieto, Marcos Vera-Hernandez and Grant Miller 
//  date: 		February 20 / 2024
************************************************************************************************************


** Opening the log file
	log using "${logs_cs}/Bootstrapp_LLR_OUT.smcl", append

	
** Loading the dataset	
	use PersonaBasicaID mes mmw OUT using  "${bases_cs}/Panel_2011_2018_full_OUT.dta", clear 

		** Variables meaning
			* PersonaBasicaID: anonimized ID for each subject in the database
			* mmw: monthly minimum wage
			* mes: month of the observation registered
			* Out: Number of the outpatient services per month for each register



*******************************************************************************************************	
*******************		 Centering the MMW and generating capy variable

	gen mmwC=mmw-5
	gen byte Copay=1 if mmwC>=0 & mmwC!=.
	replace Copay=0 if mmwC<0 & mmwC!=.

	
* Setting a loop to calculate the 500 distribution of each estimator for each lag

	foreach lag of numlist 0(1)95 {

*******************		 Setting data in time series 
	display "*******************before the tsset time: $S_DATE  $S_TIME **********************"
	
	
	sort PersonaBasicaID mes 
	tsset PersonaBasicaID  mes
		
		
	display "*******************after the tsset time: $S_DATE  $S_TIME **********************"
		
	** Generating the lagged mmw
	
		by PersonaBasicaID : gen  mmwC_lag`lag'= l`lag'.mmwC

	**************************************** RD ANALYSIS

		display "******************* Outcome: OUT  Lag = `lag' time: $S_DATE  $S_TIME **********************"

		rdrobust OUT mmwC_lag`lag', c(0) vce(cluster PersonaBasicaID)  masspoints(off)
		
		* Collecting the original estimator
		local bw=e(h_l)
		
		* Setting again the time series to run the bootstrap
		tsset, clear
		
		display "******************* Outcome: OUT Lag = `lag'  BW=`bw' time: $S_DATE  $S_TIME **********************"
		
		bootstrap, reps(500) seed(1) cluster(PersonaBasicaID) idcluster(newID) saving(${resultados_cs}/OUT/Boots_LLR_OUT_lag`lag'.dta, every(100) replace): rd OUT mmwC_lag`lag', z0(0) mbw(100) cluster(PersonaBasicaID)  bwidth(`bw')
		
		display "******************* End of the bootstrap time: $S_DATE  $S_TIME **********************"
		
		** Collecting the results of the bootstrap process 
		matrix boot_or_estim=J(1,3,-99)
		matrix colnames boot_or_estim = lag coeficient boot_SE		
		matrix boot_or_estim[1,1]= `lag'
		matrix boot_or_estim[1,2]= e(b)
		matrix boot_or_estim[1,3]= e(se)
		
		** Exporting the data
		putexcel set "${resultados_cs}/PI/Boots_Ori_Est_OUT_LLR_lag`lag'", replace
		putexcel A1=matrix(boot_or_estim), names
		
		** Dropping the lagged mmw variable created at the beggining of the loop
		drop mmwC_lag`lag' 

}


******************************************************************************************************
log close
		
		
		
		
		
		
