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
** Aim: impact analyses using a Weibull survival analysis to calculate mortality to 8 years

//  task:      	ITT-LLR analyses for mortality
//  authors:    Giancarlo Buitrago, Javier Amaya-Nieto, Marcos Vera-Hernandez and Grant Miller 
//  date: 		February 20 / 2024
************************************************************************************************************

	log using "${logs_cs}/log_ITT_LLR_Mortality.smcl", replace


** Loading the data(not publicly available on request to the Ministry of Health of Colombia)
		use PersonaBasicaID age women eps_public region mes mmw death m_death using "${bases_cs}/MMW_long_survival_2011-2018.dta", clear
		
	** Variables meaning
		* PersonaBasicaID: anonimized ID for each subject in the database
		* Age: age of each subject in the data set
		* women: dichotomous variable indicating if the person is a woman or not
		* eps_public: dichotomous variable indicating if the person was affiliated to a public insurer
		* region: categoric variable indicating the region in which the person is afilliated 
		* mmw: monthly minimum wage
		* mes: month of the observation registered
		* Death: dichotomous variable indicating if the person is death or not
		* m_death: month of the death during the follow up period


** Selecting subjects arround the threshold to guaranty comparable populations and creating useful variables for the analysis


	keep if 0.5>abs(5-mmw)

	gen Copay=1 if mmw>5 & mmw!=.
	replace Copay=0 if mmw<=5
	
	sort PersonaBasicaID mes
	
	gen time= m_death- mes if death==1
	replace time=708-mes if death==0
	
	stset time , failure(death) scale(12)
	
	gen mmw_2=mmw*mmw													
	
	by PersonaBasicaID : gen c=_n
	by PersonaBasicaID: egen max_c=max(c)
			
******************************************************************************************************	
*******************		 setting the matrix to collect results

	matrix index=J(4,6,-99)
	matrix colnames index = HR CI_low CI_high p obs Months 			
	matrix list index
	local irow=0				



	** running the regression analysis including subjects with at least 12 months of follow up
	streg Copay age women eps_public i.region mmw mmw_2 if max_c>=12, vce(cluster PersonaBasicaID) distribution(weibull)
	
	** Making predictions of the survival
	predict double S0_1, csurv
	
	** Collecting results
	local ++irow
	matrix index[`irow',1]=el(r(table),1,1)
	matrix index[`irow',2]=el(r(table),5,1)
	matrix index[`irow',3]=el(r(table),6,1)
	matrix index[`irow',4]=el(r(table),4,1)
	matrix index[`irow',5]=e(N)	
	matrix index[`irow',6]=`dur'

** Identifying the median survival in those who reached a 8 years follwo up
	* Calculating the mortality for each 10 000 persons under the threshold
	sum S0_1 if time==96 & Copay ==0, det
	local t0s = (1-r(p50))*10000

	* Calculating the mortality for each 10 000 persons over the threshold
	sum S0_1 if time==96 & Copay ==1, det
	local t1s= (1-r(p50))*10000

	* Calculating the absolute difference
	local rd=`t1s'-`t0s'

	* Exporting results
	outreg2 using "${resultados_cs}/Mortality_prediction_RD_weibull_", append excel ctitle(2MMW) keep(Copay) eform addtext(BW, "0.5 MMW", Order, "Second") addstat(8-year mortality risk below threshold (per 10000), `t0s', 8-year mortality risk above threshold (per 10000), `t1s', Absolute difference of mortality risk (per 10000), `rd')

	
drop S0_1
	

** Exportando the results matrix to excel. 
	
	putexcel set "${resultados_cs}/Survival_ID_weibull_coefplot", replace
	putexcel A1=matrix(index), names		

	import excel "${resultados_cs}/Survival_ID_weibull_coefplot", sheet("Sheet1") firstrow clear

******************************************************************************************************************************************************
******************************* Making the coef plot


	import excel "${resultados_cs}/Survival_ID_weibull_coefplot", sheet("Sheet1") firstrow clear

	** Graph

	twoway rcap CI_high CI_low  Months, lcolor(black) || scatter HR Months,	graphregion(color(white)) ///
	mlabcolor(black) msymbol(S) mcolor(black) msize(vsmall) mlabel(HR) ///
	xlabel(, labsize(small)) xlabel(,  valuelabel) xlabel(, angle(0)) xtitle("Months of exposure", size(small)) /// 
	yline(1, lpattern(dash) lcolor(black) lwidth(medium))  ///  
	ytitle(Hazard Ratio, size(small))  ///
	yscale(range(0.9 1.4)) ylabel(0.9(0.1)1.4, labsize(small) format(%9.1f))  ///
	legend(label(1 "95% CI")) legend(size(small))

	graph save "${resultados_cs}/Survival_ID_weibull_coefplot.gph", replace 
	graph export "${resultados_cs}/Survival_ID_weibull_coefplot.png", as(png) name("Graph") replace


*********************************
log close










